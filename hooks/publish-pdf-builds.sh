#!/usr/bin/env bash
# Build mumu_cancellation.pdf and rewrite the orphan branch `pdf-builds`
# to a single-commit tip containing only that PDF at the branch root.
# Force-pushes to remotes `github` and `origin` (via proxychains4 -q).
#
# Usage:
#   ./hooks/publish-pdf-builds.sh           # make + publish
#   ./hooks/publish-pdf-builds.sh --no-build  # publish existing PDF
#   ./hooks/publish-pdf-builds.sh --dry-run   # build (unless --no-build) and
#                                             # assemble orphan commit locally;
#                                             # do not push
#
# Safe for post-commit: failures print a warning and exit 0.

set -u

ROOT="$(git -C "${BASH_SOURCE[0]%/*}/.." rev-parse --show-toplevel 2>/dev/null)" || exit 0
cd "$ROOT" || exit 0

REPORT_DIR="$ROOT/06-mc-multiperipheral/report"
PDF_SRC="$REPORT_DIR/mumu_cancellation.pdf"
PDF_NAME="mumu_cancellation.pdf"
BRANCH="pdf-builds"
REMOTES=(github origin)

NO_BUILD=0
DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --no-build) NO_BUILD=1 ;;
    --dry-run)  DRY_RUN=1 ;;
    -h|--help)
      sed -n '2,16p' "$0" | sed 's/^# \?//'
      exit 0
      ;;
    *)
      echo "publish-pdf-builds: unknown option: $arg" >&2
      exit 0
      ;;
  esac
done

if [[ -d "$(git rev-parse --git-path rebase-merge 2>/dev/null)" ]] || \
   [[ -d "$(git rev-parse --git-path rebase-apply 2>/dev/null)" ]]; then
  echo "publish-pdf-builds: rebase in progress — skipping." >&2
  exit 0
fi

if [[ "$NO_BUILD" -eq 0 ]]; then
  echo "publish-pdf-builds: building $PDF_NAME ..."
  if ! make -C "$REPORT_DIR" pdf; then
    echo "publish-pdf-builds: ERROR: make pdf failed — not publishing." >&2
    exit 0
  fi
fi

if [[ ! -s "$PDF_SRC" ]]; then
  echo "publish-pdf-builds: ERROR: missing or empty $PDF_SRC" >&2
  exit 0
fi

SHORT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
COMMIT_MSG="$(git log -1 --pretty=format:'%s' 2>/dev/null || echo '(no commit)')"
PAGES="$(pdfinfo "$PDF_SRC" 2>/dev/null | awk '/^Pages:/{print $2}')"
PAGES="${PAGES:-?}"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

STAGE="$TMPDIR/stage"
mkdir -p "$STAGE"
cp "$PDF_SRC" "$STAGE/$PDF_NAME"

# Local orphan commit via commit-tree (does not touch the main worktree).
blob="$(git hash-object -w "$STAGE/$PDF_NAME")" || exit 0
tree="$(printf '100644 blob %s\t%s\n' "$blob" "$PDF_NAME" | git mktree)" || exit 0
msg="pdf: ${SHORT_SHA} — ${COMMIT_MSG} [${PDF_NAME}]"
orphan="$(git commit-tree "$tree" -m "$msg")" || exit 0

echo "publish-pdf-builds: orphan tip $(git rev-parse --short "$orphan") ($PDF_NAME, ${PAGES} pages) from main@${SHORT_SHA}"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "publish-pdf-builds: dry-run — not updating refs or remotes."
  echo "publish-pdf-builds: would: git update-ref refs/heads/${BRANCH} ${orphan}"
  for r in "${REMOTES[@]}"; do
    if git remote get-url "$r" >/dev/null 2>&1; then
      echo "publish-pdf-builds: would: proxychains4 -q git push --force $r ${orphan}:refs/heads/${BRANCH}"
    else
      echo "publish-pdf-builds: dry-run: remote '$r' missing — would skip"
    fi
  done
  exit 0
fi

git update-ref "refs/heads/${BRANCH}" "$orphan" || {
  echo "publish-pdf-builds: ERROR: update-ref ${BRANCH} failed" >&2
  exit 0
}

push_one() {
  local remote="$1"
  if ! git remote get-url "$remote" >/dev/null 2>&1; then
    echo "publish-pdf-builds: remote '$remote' not configured — skip." >&2
    return 0
  fi
  echo "publish-pdf-builds: force-pushing ${BRANCH} -> ${remote} ..."
  if ! proxychains4 -q git push --force "$remote" "refs/heads/${BRANCH}:refs/heads/${BRANCH}"; then
    echo "publish-pdf-builds: ERROR: push to ${remote} failed" >&2
    return 0
  fi
  echo "publish-pdf-builds: ${remote}/${BRANCH} updated."
}

for r in "${REMOTES[@]}"; do
  push_one "$r"
done

exit 0
