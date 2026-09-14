#!/bin/sh

set -eu

development_root="/Users/zhikuan_quan/Documents/ALRD"
github_root="/Users/zhikuan_quan/Documents/GitHub/ALRD"
git_bin="/Users/zhikuan_quan/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/fallback/git"

if [ ! -d "$development_root/.git" ]; then
  echo "Development repository not found: $development_root" >&2
  exit 1
fi

if [ ! -d "$github_root/.git" ]; then
  echo "GitHub-connected repository not found: $github_root" >&2
  exit 1
fi

echo "Development workspace status"
"$git_bin" -C "$development_root" status --short --branch

echo
echo "GitHub checkout status"
"$git_bin" -C "$github_root" status --short --branch

echo
echo "Content differences"
diff -qr \
  --exclude=.git \
  --exclude=.DS_Store \
  --exclude=.Rhistory \
  --exclude=tmp \
  --exclude=output \
  "$development_root" "$github_root" || true
