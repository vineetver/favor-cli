#!/usr/bin/env sh
set -e

REPO="vineetver/cohort-cli"
INSTALL_DIR="${COHORT_INSTALL_DIR:-$HOME/.local/bin}"

main() {
    need_cmd curl
    need_cmd tar

    OS=$(uname -s | tr '[:upper:]' '[:lower:]')
    ARCH=$(uname -m)

    case "$OS" in
        linux)
            case "$ARCH" in
                x86_64|amd64) ARTIFACT="cohort-x86_64-linux" ;;
                *) err "Unsupported architecture: $ARCH. Build from source: cargo install --git https://github.com/$REPO" ;;
            esac
            ;;
        darwin)
            case "$ARCH" in
                x86_64)        ARTIFACT="cohort-x86_64-macos" ;;
                arm64|aarch64) ARTIFACT="cohort-aarch64-macos" ;;
                *) err "Unsupported architecture: $ARCH" ;;
            esac
            ;;
        *) err "Unsupported OS: $OS. Build from source: cargo install --git https://github.com/$REPO" ;;
    esac

    LATEST=$(curl -fsSL "https://api.github.com/repos/$REPO/releases/latest" | grep '"tag_name"' | head -1 | sed 's/.*"tag_name": *"//;s/".*//')
    if [ -z "$LATEST" ]; then
        err "Could not determine latest release. Check https://github.com/$REPO/releases"
    fi

    URL="https://github.com/$REPO/releases/download/$LATEST/$ARTIFACT.tar.gz"

    echo "Installing cohort $LATEST ($ARTIFACT)..."
    echo "  From: $URL"
    echo "  To:   $INSTALL_DIR/cohort"

    TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    curl -fsSL "$URL" -o "$TMPDIR/$ARTIFACT.tar.gz"
    tar xzf "$TMPDIR/$ARTIFACT.tar.gz" -C "$TMPDIR"

    mkdir -p "$INSTALL_DIR"
    mv "$TMPDIR/cohort" "$INSTALL_DIR/cohort"
    chmod +x "$INSTALL_DIR/cohort"

    # Add to PATH if not already there
    case ":$PATH:" in
        *":$INSTALL_DIR:"*) ;;
        *)
            SHELL_NAME=$(basename "$SHELL" 2>/dev/null || echo "bash")
            case "$SHELL_NAME" in
                zsh)  RC="$HOME/.zshrc" ;;
                fish) RC="$HOME/.config/fish/config.fish" ;;
                *)    RC="$HOME/.bashrc" ;;
            esac

            if [ -f "$RC" ] && grep -q "$INSTALL_DIR" "$RC" 2>/dev/null; then
                : # already in rc file, just not in current session
            else
                echo "export PATH=\"$INSTALL_DIR:\$PATH\"" >> "$RC"
                echo "  Added $INSTALL_DIR to PATH in $RC"
            fi
            export PATH="$INSTALL_DIR:$PATH"
            ;;
    esac

    echo ""
    echo "Installed cohort $LATEST"
    echo "Run 'cohort setup' to configure."
    echo ""
    echo "To uninstall: cohort uninstall (or just rm $INSTALL_DIR/cohort)"
}

need_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        err "need '$1' (not found in PATH)"
    fi
}

err() {
    echo "error: $1" >&2
    exit 1
}

main
