# Git hooks

## `post-commit` → orphan `pdf-builds`

When a commit on `main` touches `06-mc-multiperipheral/report/`, the hook
rebuilds `mumu_cancellation.pdf` and force-pushes a single-commit orphan
branch `pdf-builds` (PDF at branch root) to remotes `github` and `origin`.

Raw URLs (after the branch exists):

- https://github.com/mkhellat/nikhef-form-lab/raw/pdf-builds/mumu_cancellation.pdf
- https://codeberg.org/mkhellat/nikhef-form-lab/raw/branch/pdf-builds/mumu_cancellation.pdf

### Install (per clone)

```bash
cp hooks/post-commit .git/hooks/post-commit
chmod +x .git/hooks/post-commit hooks/publish-pdf-builds.sh
```

Git does not auto-install hooks from the working tree; this copy step is
required after each fresh clone.

### Manual publish / dry-run

```bash
./hooks/publish-pdf-builds.sh --dry-run   # build + show what would push
./hooks/publish-pdf-builds.sh             # build + update ref + push both remotes
./hooks/publish-pdf-builds.sh --no-build  # reuse existing PDF
make -C 06-mc-multiperipheral/report pdf-builds
```

Pushes use `proxychains4 -q` (same network constraint as ordinary remotes).
Failures in the hook print a warning and exit 0 — they do not undo the
source commit.
