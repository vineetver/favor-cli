use crossterm::event::{KeyCode, KeyModifiers};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ActionScope {
    Global,
    Workspace,
    Setup,
    Transform,
    FilePicker,
    Run,
    Help,
    Palette,
    Variant,
}

impl ActionScope {
    pub fn title(&self) -> &'static str {
        match self {
            Self::Global => "Global",
            Self::Workspace => "Workspace",
            Self::Setup => "Setup",
            Self::Transform => "Transform",
            Self::FilePicker => "File picker",
            Self::Run => "Run",
            Self::Help => "Help",
            Self::Palette => "Palette",
            Self::Variant => "Variant browser",
        }
    }

    pub fn ordered() -> &'static [ActionScope] {
        &[
            ActionScope::Global,
            ActionScope::Workspace,
            ActionScope::Setup,
            ActionScope::Transform,
            ActionScope::FilePicker,
            ActionScope::Run,
            ActionScope::Help,
            ActionScope::Palette,
            ActionScope::Variant,
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Action {
    Quit,
    OpenPalette,
    OpenHelp,
    ClosePalette,
    ClosePaletteAndRun,

    WorkspaceUp,
    WorkspaceDown,
    WorkspaceRescan,
    WorkspaceOpenSetup,
    WorkspaceOpenFocused,

    TransformNextField,
    TransformPrevField,
    TransformActivate,
    TransformToggleBool,
    TransformClearField,
    TransformCancel,

    PickerUp,
    PickerDown,
    PickerInto,
    PickerParent,
    PickerSelect,
    PickerCreate,
    PickerTypePath,
    PickerCancel,

    RunCancelRequest,
    RunReturn,

    HelpClose,

    SetupCancel,
    SetupConfirm,
    SetupNext,
    SetupPrev,
    SetupToggle,
    SetupToggleAll,

    VariantScrollRowUp,
    VariantScrollRowDown,
    VariantPrevRowGroup,
    VariantNextRowGroup,
    VariantJumpStart,
    VariantJumpEnd,
    VariantOpenFilter,
    VariantFilterSubmit,
    VariantFilterClear,
    VariantOpenColumnPicker,
    VariantColumnToggle,
    VariantColumnPickerClose,
    VariantOpenCarrierView,
    VariantCloseCarrierView,
    VariantClose,
}

pub type KeyBinding = (KeyCode, KeyModifiers);

const ACTIONS_ALL: &[Action] = &[
    Action::Quit,
    Action::OpenPalette,
    Action::OpenHelp,
    Action::ClosePalette,
    Action::ClosePaletteAndRun,
    Action::WorkspaceUp,
    Action::WorkspaceDown,
    Action::WorkspaceRescan,
    Action::WorkspaceOpenSetup,
    Action::WorkspaceOpenFocused,
    Action::SetupPrev,
    Action::SetupNext,
    Action::SetupConfirm,
    Action::SetupCancel,
    Action::SetupToggle,
    Action::SetupToggleAll,
    Action::TransformPrevField,
    Action::TransformNextField,
    Action::TransformActivate,
    Action::TransformToggleBool,
    Action::TransformClearField,
    Action::TransformCancel,
    Action::PickerUp,
    Action::PickerDown,
    Action::PickerInto,
    Action::PickerParent,
    Action::PickerSelect,
    Action::PickerCreate,
    Action::PickerTypePath,
    Action::PickerCancel,
    Action::RunCancelRequest,
    Action::RunReturn,
    Action::HelpClose,
    Action::VariantScrollRowUp,
    Action::VariantScrollRowDown,
    Action::VariantPrevRowGroup,
    Action::VariantNextRowGroup,
    Action::VariantJumpStart,
    Action::VariantJumpEnd,
    Action::VariantOpenFilter,
    Action::VariantFilterSubmit,
    Action::VariantFilterClear,
    Action::VariantOpenColumnPicker,
    Action::VariantColumnToggle,
    Action::VariantColumnPickerClose,
    Action::VariantOpenCarrierView,
    Action::VariantCloseCarrierView,
    Action::VariantClose,
];

impl Action {
    pub fn scope(&self) -> ActionScope {
        match self {
            Self::Quit
            | Self::OpenPalette
            | Self::OpenHelp
            | Self::ClosePalette
            | Self::ClosePaletteAndRun => ActionScope::Global,

            Self::WorkspaceUp
            | Self::WorkspaceDown
            | Self::WorkspaceRescan
            | Self::WorkspaceOpenSetup
            | Self::WorkspaceOpenFocused => ActionScope::Workspace,

            Self::SetupCancel
            | Self::SetupConfirm
            | Self::SetupNext
            | Self::SetupPrev
            | Self::SetupToggle
            | Self::SetupToggleAll => ActionScope::Setup,

            Self::TransformNextField
            | Self::TransformPrevField
            | Self::TransformActivate
            | Self::TransformToggleBool
            | Self::TransformClearField
            | Self::TransformCancel => ActionScope::Transform,

            Self::PickerUp
            | Self::PickerDown
            | Self::PickerInto
            | Self::PickerParent
            | Self::PickerSelect
            | Self::PickerCreate
            | Self::PickerTypePath
            | Self::PickerCancel => ActionScope::FilePicker,

            Self::RunCancelRequest | Self::RunReturn => ActionScope::Run,

            Self::HelpClose => ActionScope::Help,

            Self::VariantScrollRowUp
            | Self::VariantScrollRowDown
            | Self::VariantPrevRowGroup
            | Self::VariantNextRowGroup
            | Self::VariantJumpStart
            | Self::VariantJumpEnd
            | Self::VariantOpenFilter
            | Self::VariantFilterSubmit
            | Self::VariantFilterClear
            | Self::VariantOpenColumnPicker
            | Self::VariantColumnToggle
            | Self::VariantColumnPickerClose
            | Self::VariantOpenCarrierView
            | Self::VariantCloseCarrierView
            | Self::VariantClose => ActionScope::Variant,
        }
    }

    pub fn title(&self) -> &'static str {
        match self {
            Self::Quit => "Quit",
            Self::OpenPalette => "Command palette",
            Self::OpenHelp => "Help",
            Self::ClosePalette => "Close palette",
            Self::ClosePaletteAndRun => "Run selected",

            Self::WorkspaceUp => "Move focus up",
            Self::WorkspaceDown => "Move focus down",
            Self::WorkspaceRescan => "Rescan workspace",
            Self::WorkspaceOpenSetup => "Open setup",
            Self::WorkspaceOpenFocused => "Open focused artifact",

            Self::TransformNextField => "Next field",
            Self::TransformPrevField => "Previous field",
            Self::TransformActivate => "Activate field",
            Self::TransformToggleBool => "Toggle bool",
            Self::TransformClearField => "Clear field",
            Self::TransformCancel => "Back",

            Self::PickerUp => "Picker up",
            Self::PickerDown => "Picker down",
            Self::PickerInto => "Enter directory",
            Self::PickerParent => "Parent directory",
            Self::PickerSelect => "Select current directory",
            Self::PickerCreate => "Create directory",
            Self::PickerTypePath => "Type path",
            Self::PickerCancel => "Cancel picker",

            Self::RunCancelRequest => "Cancel run",
            Self::RunReturn => "Return",

            Self::HelpClose => "Close help",

            Self::SetupCancel => "Cancel setup",
            Self::SetupConfirm => "Confirm",
            Self::SetupNext => "Next",
            Self::SetupPrev => "Previous",
            Self::SetupToggle => "Toggle",
            Self::SetupToggleAll => "Toggle all",

            Self::VariantScrollRowUp => "Scroll up",
            Self::VariantScrollRowDown => "Scroll down",
            Self::VariantPrevRowGroup => "Previous row group",
            Self::VariantNextRowGroup => "Next row group",
            Self::VariantJumpStart => "Jump to start",
            Self::VariantJumpEnd => "Jump to end",
            Self::VariantOpenFilter => "Filter",
            Self::VariantFilterSubmit => "Apply filter",
            Self::VariantFilterClear => "Clear filter",
            Self::VariantOpenColumnPicker => "Columns",
            Self::VariantColumnToggle => "Toggle column",
            Self::VariantColumnPickerClose => "Close column picker",
            Self::VariantOpenCarrierView => "Show carriers",
            Self::VariantCloseCarrierView => "Close carrier panel",
            Self::VariantClose => "Close browser",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::Quit => "exit cohort",
            Self::OpenPalette => "open the command palette",
            Self::OpenHelp => "show keys for the current screen",
            Self::ClosePalette => "dismiss the palette",
            Self::ClosePaletteAndRun => "run the highlighted action",

            Self::WorkspaceUp => "select the previous artifact",
            Self::WorkspaceDown => "select the next artifact",
            Self::WorkspaceRescan => "rescan roots for artifacts",
            Self::WorkspaceOpenSetup => "configure tier, root, packs, environment",
            Self::WorkspaceOpenFocused => "open a transform for the focused artifact",

            Self::TransformNextField => "move to the next form field",
            Self::TransformPrevField => "move to the previous form field",
            Self::TransformActivate => "edit, cycle, or run the focused field",
            Self::TransformToggleBool => "toggle a bool field",
            Self::TransformClearField => "clear the focused field",
            Self::TransformCancel => "leave the form",

            Self::PickerUp => "move selection up in the directory list",
            Self::PickerDown => "move selection down in the directory list",
            Self::PickerInto => "descend into the highlighted directory",
            Self::PickerParent => "go to the parent directory",
            Self::PickerSelect => "select the current directory",
            Self::PickerCreate => "create the current path and select it",
            Self::PickerTypePath => "type a path directly",
            Self::PickerCancel => "close the picker",

            Self::RunCancelRequest => "request cancellation (advisory)",
            Self::RunReturn => "return to the previous screen",

            Self::HelpClose => "close help",

            Self::SetupCancel => "exit the setup wizard",
            Self::SetupConfirm => "confirm the current selection",
            Self::SetupNext => "move to the next option",
            Self::SetupPrev => "move to the previous option",
            Self::SetupToggle => "toggle the highlighted option",
            Self::SetupToggleAll => "toggle every option in the list",

            Self::VariantScrollRowUp => "move row focus up within the row group",
            Self::VariantScrollRowDown => "move row focus down within the row group",
            Self::VariantPrevRowGroup => "jump to the previous parquet row group",
            Self::VariantNextRowGroup => "jump to the next parquet row group",
            Self::VariantJumpStart => "jump to the first row group",
            Self::VariantJumpEnd => "jump to the last row group",
            Self::VariantOpenFilter => "open the filter bar",
            Self::VariantFilterSubmit => "compile and apply the filter expression",
            Self::VariantFilterClear => "drop the active filter",
            Self::VariantOpenColumnPicker => "show or hide individual columns",
            Self::VariantColumnToggle => "toggle the focused column visibility",
            Self::VariantColumnPickerClose => "close the column picker overlay",
            Self::VariantOpenCarrierView => "show carrier samples for the focused variant",
            Self::VariantCloseCarrierView => "close the carrier panel",
            Self::VariantClose => "leave the variant browser",
        }
    }

    pub fn all() -> &'static [Action] {
        ACTIONS_ALL
    }

    pub fn default_key(&self) -> Option<KeyBinding> {
        let none = KeyModifiers::NONE;
        match self {
            Self::Quit => Some((KeyCode::Char('q'), none)),
            Self::OpenPalette => Some((KeyCode::Char('p'), KeyModifiers::CONTROL)),
            Self::OpenHelp => Some((KeyCode::Char('?'), none)),
            Self::ClosePalette => Some((KeyCode::Esc, none)),
            Self::ClosePaletteAndRun => Some((KeyCode::Enter, none)),

            Self::WorkspaceUp => Some((KeyCode::Char('k'), none)),
            Self::WorkspaceDown => Some((KeyCode::Char('j'), none)),
            Self::WorkspaceRescan => Some((KeyCode::Char('r'), none)),
            Self::WorkspaceOpenSetup => Some((KeyCode::Char('s'), none)),
            Self::WorkspaceOpenFocused => Some((KeyCode::Enter, none)),

            Self::TransformNextField => Some((KeyCode::Tab, none)),
            Self::TransformPrevField => Some((KeyCode::BackTab, KeyModifiers::SHIFT)),
            Self::TransformActivate => Some((KeyCode::Enter, none)),
            Self::TransformToggleBool => Some((KeyCode::Char(' '), none)),
            Self::TransformClearField => Some((KeyCode::Backspace, none)),
            Self::TransformCancel => Some((KeyCode::Esc, none)),

            Self::PickerUp => Some((KeyCode::Char('k'), none)),
            Self::PickerDown => Some((KeyCode::Char('j'), none)),
            Self::PickerInto => Some((KeyCode::Enter, none)),
            Self::PickerParent => Some((KeyCode::Left, none)),
            Self::PickerSelect => Some((KeyCode::Char(' '), none)),
            Self::PickerCreate => Some((KeyCode::Char('c'), none)),
            Self::PickerTypePath => Some((KeyCode::Char('/'), none)),
            Self::PickerCancel => Some((KeyCode::Esc, none)),

            Self::RunCancelRequest => Some((KeyCode::Char('c'), none)),
            Self::RunReturn => Some((KeyCode::Enter, none)),

            Self::HelpClose => Some((KeyCode::Esc, none)),

            Self::SetupCancel => Some((KeyCode::Esc, none)),
            Self::SetupConfirm => Some((KeyCode::Enter, none)),
            Self::SetupNext => Some((KeyCode::Down, none)),
            Self::SetupPrev => Some((KeyCode::Up, none)),
            Self::SetupToggle => Some((KeyCode::Char(' '), none)),
            Self::SetupToggleAll => Some((KeyCode::Char('a'), none)),

            Self::VariantScrollRowUp => Some((KeyCode::Char('k'), none)),
            Self::VariantScrollRowDown => Some((KeyCode::Char('j'), none)),
            Self::VariantPrevRowGroup => Some((KeyCode::Char('{'), none)),
            Self::VariantNextRowGroup => Some((KeyCode::Char('}'), none)),
            Self::VariantJumpStart => Some((KeyCode::Char('g'), none)),
            Self::VariantJumpEnd => Some((KeyCode::Char('G'), KeyModifiers::SHIFT)),
            Self::VariantOpenFilter => Some((KeyCode::Char('/'), none)),
            Self::VariantFilterSubmit => Some((KeyCode::Enter, none)),
            Self::VariantFilterClear => Some((KeyCode::Char('x'), none)),
            Self::VariantOpenColumnPicker => Some((KeyCode::Char('!'), none)),
            Self::VariantColumnToggle => Some((KeyCode::Char(' '), none)),
            Self::VariantColumnPickerClose => Some((KeyCode::Esc, none)),
            Self::VariantOpenCarrierView => Some((KeyCode::Char('c'), none)),
            Self::VariantCloseCarrierView => Some((KeyCode::Esc, none)),
            Self::VariantClose => Some((KeyCode::Char('q'), none)),
        }
    }
}

pub fn format_binding(code: KeyCode, mods: KeyModifiers) -> String {
    let mut parts: Vec<String> = Vec::new();
    if mods.contains(KeyModifiers::CONTROL) {
        parts.push("Ctrl".to_string());
    }
    if mods.contains(KeyModifiers::ALT) {
        parts.push("Alt".to_string());
    }
    if mods.contains(KeyModifiers::SHIFT) && !matches!(code, KeyCode::BackTab) {
        parts.push("Shift".to_string());
    }
    let key = match code {
        KeyCode::Char(' ') => "Space".to_string(),
        KeyCode::Char(c) => c.to_string(),
        KeyCode::Enter => "Enter".to_string(),
        KeyCode::Esc => "Esc".to_string(),
        KeyCode::Tab => "Tab".to_string(),
        KeyCode::BackTab => "Shift-Tab".to_string(),
        KeyCode::Backspace => "Backspace".to_string(),
        KeyCode::Up => "Up".to_string(),
        KeyCode::Down => "Down".to_string(),
        KeyCode::Left => "Left".to_string(),
        KeyCode::Right => "Right".to_string(),
        KeyCode::Home => "Home".to_string(),
        KeyCode::End => "End".to_string(),
        KeyCode::PageUp => "PageUp".to_string(),
        KeyCode::PageDown => "PageDown".to_string(),
        KeyCode::Delete => "Delete".to_string(),
        KeyCode::Insert => "Insert".to_string(),
        other => format!("{other:?}"),
    };
    parts.push(key);
    parts.join("-")
}

pub struct KeyMap {
    entries: Vec<(KeyCode, KeyModifiers, Action)>,
}

impl KeyMap {
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    pub fn bind(mut self, code: KeyCode, mods: KeyModifiers, action: Action) -> Self {
        self.entries.push((code, mods, action));
        self
    }

    pub fn lookup(&self, code: KeyCode, mods: KeyModifiers) -> Option<Action> {
        for (c, m, a) in &self.entries {
            if *c == code && match_mods(*m, mods) {
                return Some(*a);
            }
        }
        None
    }
}

impl Default for KeyMap {
    fn default() -> Self {
        Self::new()
    }
}

fn match_mods(want: KeyModifiers, got: KeyModifiers) -> bool {
    let mask = KeyModifiers::CONTROL | KeyModifiers::ALT;
    (want & mask) == (got & mask)
}
