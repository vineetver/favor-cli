use arrow::array::{
    Array, BooleanArray, Float32Array, Float64Array, Int32Array, Int64Array, StringArray,
    UInt32Array, UInt64Array,
};
use arrow::compute::kernels::cmp;
use arrow::datatypes::DataType;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::{ArrowPredicateFn, ArrowReaderMetadata, RowFilter};
use parquet::arrow::ProjectionMask;

use crate::error::CohortError;

use super::parquet_scroller::RowFilterFactory;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Eq,
    Neq,
    Lt,
    Lte,
    Gt,
    Gte,
    Contains,
}

#[derive(Debug, Clone)]
pub enum Literal {
    Int(i64),
    Float(f64),
    Str(String),
}

#[derive(Debug, Clone)]
pub struct Clause {
    pub column: String,
    pub op: Op,
    pub literal: Literal,
}

#[derive(Debug, Clone)]
pub struct CompiledFilter {
    text: String,
    clauses: Vec<Clause>,
}

impl CompiledFilter {
    pub fn parse(input: &str) -> Result<Self, String> {
        let text = input.trim().to_string();
        if text.is_empty() {
            return Err("filter is empty".into());
        }
        let mut p = Parser::new(&text);
        let clauses = p.parse_filter()?;
        Ok(Self { text, clauses })
    }

    pub fn clauses(&self) -> &[Clause] {
        &self.clauses
    }
}

struct Parser<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Self {
        Self { src, pos: 0 }
    }

    fn parse_filter(&mut self) -> Result<Vec<Clause>, String> {
        let mut out = Vec::new();
        out.push(self.parse_clause()?);
        loop {
            self.skip_ws();
            if self.eof() {
                break;
            }
            if !self.consume_keyword("AND") {
                return Err(format!("expected AND at offset {}", self.pos));
            }
            out.push(self.parse_clause()?);
        }
        Ok(out)
    }

    fn parse_clause(&mut self) -> Result<Clause, String> {
        self.skip_ws();
        let column = self.parse_ident()?;
        self.skip_ws();
        let op = self.parse_op()?;
        self.skip_ws();
        let literal = self.parse_literal()?;
        Ok(Clause { column, op, literal })
    }

    fn parse_ident(&mut self) -> Result<String, String> {
        let start = self.pos;
        let bytes = self.src.as_bytes();
        if start >= bytes.len() || !is_ident_start(bytes[start]) {
            return Err(format!("expected column name at offset {start}"));
        }
        while self.pos < bytes.len() && is_ident_cont(bytes[self.pos]) {
            self.pos += 1;
        }
        Ok(self.src[start..self.pos].to_string())
    }

    fn parse_op(&mut self) -> Result<Op, String> {
        let bytes = self.src.as_bytes();
        let rest = &self.src[self.pos..];
        if rest.starts_with("!=") {
            self.pos += 2;
            return Ok(Op::Neq);
        }
        if rest.starts_with("<=") {
            self.pos += 2;
            return Ok(Op::Lte);
        }
        if rest.starts_with(">=") {
            self.pos += 2;
            return Ok(Op::Gte);
        }
        if self.consume_keyword("contains") {
            return Ok(Op::Contains);
        }
        if self.pos < bytes.len() {
            let c = bytes[self.pos];
            if c == b'=' {
                self.pos += 1;
                return Ok(Op::Eq);
            }
            if c == b'<' {
                self.pos += 1;
                return Ok(Op::Lt);
            }
            if c == b'>' {
                self.pos += 1;
                return Ok(Op::Gt);
            }
        }
        Err(format!("expected operator at offset {}", self.pos))
    }

    fn parse_literal(&mut self) -> Result<Literal, String> {
        let bytes = self.src.as_bytes();
        if self.pos >= bytes.len() {
            return Err("expected literal".into());
        }
        let c = bytes[self.pos];
        if c == b'"' || c == b'\'' {
            let quote = c;
            self.pos += 1;
            let start = self.pos;
            while self.pos < bytes.len() && bytes[self.pos] != quote {
                self.pos += 1;
            }
            if self.pos >= bytes.len() {
                return Err("unterminated string literal".into());
            }
            let s = self.src[start..self.pos].to_string();
            self.pos += 1;
            return Ok(Literal::Str(s));
        }
        let start = self.pos;
        while self.pos < bytes.len() {
            let b = bytes[self.pos];
            if b.is_ascii_whitespace() {
                break;
            }
            self.pos += 1;
        }
        if start == self.pos {
            return Err("expected literal".into());
        }
        let token = &self.src[start..self.pos];
        if let Ok(i) = token.parse::<i64>() {
            return Ok(Literal::Int(i));
        }
        if let Ok(f) = token.parse::<f64>() {
            return Ok(Literal::Float(f));
        }
        Ok(Literal::Str(token.to_string()))
    }

    fn skip_ws(&mut self) {
        let bytes = self.src.as_bytes();
        while self.pos < bytes.len() && bytes[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn eof(&self) -> bool {
        self.pos >= self.src.len()
    }

    fn consume_keyword(&mut self, kw: &str) -> bool {
        let rest = &self.src[self.pos..];
        if rest.len() < kw.len() {
            return false;
        }
        let head = &rest[..kw.len()];
        if !head.eq_ignore_ascii_case(kw) {
            return false;
        }
        let after = rest.as_bytes().get(kw.len()).copied();
        if let Some(b) = after {
            if is_ident_cont(b) {
                return false;
            }
        }
        self.pos += kw.len();
        true
    }
}

fn is_ident_start(b: u8) -> bool {
    b.is_ascii_alphabetic() || b == b'_'
}

fn is_ident_cont(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

impl RowFilterFactory for CompiledFilter {
    fn describe(&self) -> &str {
        &self.text
    }

    fn build(&self, meta: &ArrowReaderMetadata) -> Result<RowFilter, CohortError> {
        let schema = meta.schema();
        let mut field_indices: Vec<usize> = Vec::new();
        for clause in &self.clauses {
            let idx = schema.index_of(&clause.column).map_err(|_| {
                CohortError::Input(format!("filter column not found: {}", clause.column))
            })?;
            if !field_indices.contains(&idx) {
                field_indices.push(idx);
            }
            check_clause_against_field(&clause, schema.field(idx))?;
        }
        let mask =
            ProjectionMask::leaves(meta.parquet_schema(), field_indices.iter().copied());
        let clauses = self.clauses.clone();
        let predicate = ArrowPredicateFn::new(mask, move |batch: RecordBatch| {
            apply_clauses(&clauses, &batch)
        });
        Ok(RowFilter::new(vec![Box::new(predicate)]))
    }
}

fn check_clause_against_field(
    clause: &Clause,
    field: &arrow::datatypes::Field,
) -> Result<(), CohortError> {
    let dt = field.data_type();
    let is_numeric = matches!(
        dt,
        DataType::Int8
            | DataType::Int16
            | DataType::Int32
            | DataType::Int64
            | DataType::UInt8
            | DataType::UInt16
            | DataType::UInt32
            | DataType::UInt64
            | DataType::Float32
            | DataType::Float64
    );
    let is_string = matches!(
        dt,
        DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View
    );

    match (clause.op, &clause.literal) {
        (Op::Contains, _) if !is_string => Err(CohortError::Input(format!(
            "contains requires a string column: {}",
            clause.column
        ))),
        (_, Literal::Int(_) | Literal::Float(_)) if !is_numeric => Err(CohortError::Input(
            format!("numeric literal against non-numeric column: {}", clause.column),
        )),
        (_, Literal::Str(s)) if is_numeric => {
            if s.parse::<f64>().is_err() {
                Err(CohortError::Input(format!(
                    "non-numeric literal against numeric column: {} = {s}",
                    clause.column
                )))
            } else {
                Ok(())
            }
        }
        _ => Ok(()),
    }
}

fn apply_clauses(
    clauses: &[Clause],
    batch: &RecordBatch,
) -> Result<BooleanArray, arrow::error::ArrowError> {
    let mut combined: Option<BooleanArray> = None;
    for clause in clauses {
        let idx = batch
            .schema()
            .index_of(&clause.column)
            .map_err(|e| arrow::error::ArrowError::ComputeError(format!("{e}")))?;
        let column = batch.column(idx);
        let mask = eval_clause(clause, column.as_ref())?;
        combined = Some(match combined.take() {
            None => mask,
            Some(prev) => arrow::compute::and(&prev, &mask)?,
        });
    }
    Ok(combined.unwrap_or_else(|| BooleanArray::from(vec![true; batch.num_rows()])))
}

fn eval_clause(
    clause: &Clause,
    array: &dyn Array,
) -> Result<BooleanArray, arrow::error::ArrowError> {
    use arrow::error::ArrowError;
    let dt = array.data_type().clone();
    if matches!(clause.op, Op::Contains) {
        let needle = literal_to_string(&clause.literal);
        return Ok(string_contains(array, &needle)?);
    }
    let lit_f = literal_to_f64(&clause.literal);
    let lit_str = literal_to_string(&clause.literal);

    macro_rules! cmp_numeric {
        ($arr_ty:ty, $native:ty) => {{
            let arr = array
                .as_any()
                .downcast_ref::<$arr_ty>()
                .ok_or_else(|| ArrowError::ComputeError("downcast".into()))?;
            let lit = lit_f as $native;
            let scalar = <$arr_ty>::from(vec![lit; arr.len()]);
            apply_op(clause.op, arr, &scalar)?
        }};
    }

    let result = match dt {
        DataType::Int32 => cmp_numeric!(Int32Array, i32),
        DataType::Int64 => cmp_numeric!(Int64Array, i64),
        DataType::UInt32 => cmp_numeric!(UInt32Array, u32),
        DataType::UInt64 => cmp_numeric!(UInt64Array, u64),
        DataType::Float32 => cmp_numeric!(Float32Array, f32),
        DataType::Float64 => cmp_numeric!(Float64Array, f64),
        DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View => {
            let arr = array
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| ArrowError::ComputeError("string downcast".into()))?;
            let scalar = StringArray::from(vec![lit_str.clone(); arr.len()]);
            apply_op(clause.op, arr, &scalar)?
        }
        DataType::Boolean => {
            let want_true = matches!(lit_str.to_ascii_lowercase().as_str(), "true" | "1");
            let arr = array
                .as_any()
                .downcast_ref::<BooleanArray>()
                .ok_or_else(|| ArrowError::ComputeError("bool downcast".into()))?;
            let scalar = BooleanArray::from(vec![want_true; arr.len()]);
            apply_op(clause.op, arr, &scalar)?
        }
        other => {
            return Err(ArrowError::ComputeError(format!(
                "unsupported column type for filter: {other:?}"
            )))
        }
    };
    Ok(result)
}

fn apply_op(
    op: Op,
    lhs: &dyn arrow::array::Datum,
    rhs: &dyn arrow::array::Datum,
) -> Result<BooleanArray, arrow::error::ArrowError> {
    match op {
        Op::Eq => cmp::eq(lhs, rhs),
        Op::Neq => cmp::neq(lhs, rhs),
        Op::Lt => cmp::lt(lhs, rhs),
        Op::Lte => cmp::lt_eq(lhs, rhs),
        Op::Gt => cmp::gt(lhs, rhs),
        Op::Gte => cmp::gt_eq(lhs, rhs),
        Op::Contains => unreachable!("contains handled before apply_op"),
    }
}

fn string_contains(
    array: &dyn Array,
    needle: &str,
) -> Result<BooleanArray, arrow::error::ArrowError> {
    use arrow::error::ArrowError;
    let arr = array
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| ArrowError::ComputeError("string downcast for contains".into()))?;
    let mut out = Vec::with_capacity(arr.len());
    for i in 0..arr.len() {
        if arr.is_null(i) {
            out.push(false);
        } else {
            out.push(arr.value(i).contains(needle));
        }
    }
    Ok(BooleanArray::from(out))
}

fn literal_to_f64(lit: &Literal) -> f64 {
    match lit {
        Literal::Int(i) => *i as f64,
        Literal::Float(f) => *f,
        Literal::Str(s) => s.parse::<f64>().unwrap_or(f64::NAN),
    }
}

fn literal_to_string(lit: &Literal) -> String {
    match lit {
        Literal::Int(i) => i.to_string(),
        Literal::Float(f) => f.to_string(),
        Literal::Str(s) => s.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::datatypes::{Field, Schema};
    use std::sync::Arc;

    fn parse_ok(s: &str) -> CompiledFilter {
        CompiledFilter::parse(s).expect("parse should succeed")
    }

    #[test]
    fn parse_eq() {
        let f = parse_ok("a = 1");
        assert_eq!(f.clauses[0].op, Op::Eq);
    }

    #[test]
    fn parse_neq() {
        assert_eq!(parse_ok("a != 1").clauses[0].op, Op::Neq);
    }

    #[test]
    fn parse_lt() {
        assert_eq!(parse_ok("a < 1").clauses[0].op, Op::Lt);
    }

    #[test]
    fn parse_lte() {
        assert_eq!(parse_ok("a <= 1").clauses[0].op, Op::Lte);
    }

    #[test]
    fn parse_gt() {
        assert_eq!(parse_ok("a > 1").clauses[0].op, Op::Gt);
    }

    #[test]
    fn parse_gte() {
        assert_eq!(parse_ok("a >= 1").clauses[0].op, Op::Gte);
    }

    #[test]
    fn parse_contains_case_insensitive() {
        assert_eq!(parse_ok("name CONTAINS missense").clauses[0].op, Op::Contains);
        assert_eq!(parse_ok("name contains missense").clauses[0].op, Op::Contains);
    }

    #[test]
    fn parse_three_clause_and() {
        let f = parse_ok("a > 1 AND b <= 2 AND c contains x");
        assert_eq!(f.clauses.len(), 3);
        assert_eq!(f.clauses[0].column, "a");
        assert_eq!(f.clauses[1].column, "b");
        assert_eq!(f.clauses[2].column, "c");
    }

    #[test]
    fn parse_quoted_string_literal() {
        let f = parse_ok("name = \"hello world\"");
        match &f.clauses[0].literal {
            Literal::Str(s) => assert_eq!(s, "hello world"),
            other => panic!("expected str literal, got {other:?}"),
        }
    }

    #[test]
    fn parse_error_missing_column() {
        assert!(CompiledFilter::parse("= 1").is_err());
    }

    #[test]
    fn parse_error_missing_op() {
        assert!(CompiledFilter::parse("a 1").is_err());
    }

    #[test]
    fn parse_error_unterminated_string() {
        assert!(CompiledFilter::parse("a = \"hello").is_err());
    }

    #[test]
    fn parse_error_empty() {
        assert!(CompiledFilter::parse("   ").is_err());
    }

    #[test]
    fn check_clause_rejects_string_in_numeric_col() {
        let field = Field::new("cadd", DataType::Float64, false);
        let clause = Clause {
            column: "cadd".into(),
            op: Op::Gt,
            literal: Literal::Str("not a number".into()),
        };
        assert!(check_clause_against_field(&clause, &field).is_err());
    }

    #[test]
    fn check_clause_accepts_string_with_numeric_text_in_numeric_col() {
        let field = Field::new("cadd", DataType::Float64, false);
        let clause = Clause {
            column: "cadd".into(),
            op: Op::Gt,
            literal: Literal::Str("25.0".into()),
        };
        assert!(check_clause_against_field(&clause, &field).is_ok());
    }

    #[test]
    fn check_clause_rejects_contains_on_numeric() {
        let field = Field::new("cadd", DataType::Float64, false);
        let clause = Clause {
            column: "cadd".into(),
            op: Op::Contains,
            literal: Literal::Str("25".into()),
        };
        assert!(check_clause_against_field(&clause, &field).is_err());
    }

    #[test]
    fn end_to_end_predicate_eval() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("cadd", DataType::Float64, false),
            Field::new("consequence", DataType::Utf8, false),
        ]));
        let cadd: Arc<dyn Array> =
            Arc::new(Float64Array::from(vec![10.0_f64, 25.5, 30.0, 5.0, 22.0]));
        let csq: Arc<dyn Array> = Arc::new(StringArray::from(vec![
            "synonymous",
            "missense_variant",
            "stop_gained",
            "missense_variant",
            "missense_variant",
        ]));
        let batch = RecordBatch::try_new(schema.clone(), vec![cadd, csq]).unwrap();

        let f = parse_ok("cadd > 20 AND consequence contains missense");
        let mask = apply_clauses(f.clauses(), &batch).unwrap();
        let v: Vec<bool> = (0..mask.len()).map(|i| mask.value(i)).collect();
        assert_eq!(v, vec![false, true, false, false, true]);
    }
}
