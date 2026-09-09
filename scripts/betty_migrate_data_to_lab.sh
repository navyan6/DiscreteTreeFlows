#!/bin/bash
# Move large data/ subtrees from Betty home → lab storage and symlink back.
# Run when home quota is exceeded; jobs still read data/* via symlinks.
set -euo pipefail
LAB="${LAB_DATA:-/vast/projects/pranam/lab/nnori/DiscreteTreeFlows_data}"
REPO="${ROOT:-$HOME/DiscreteTreeFlows}/data"
mkdir -p "$LAB"

migrate() {
  local name="$1"
  local src="$REPO/$name"
  local dst="$LAB/$name"
  if [[ -L "$src" ]]; then echo "skip (symlink) $name"; return 0; fi
  if [[ ! -d "$src" ]]; then echo "skip (missing) $name"; return 0; fi
  if [[ -e "$dst" ]]; then echo "skip (lab exists) $dst"; return 0; fi
  echo "mv $src -> $dst"
  mv "$src" "$dst"
  ln -s "$dst" "$src"
  echo "linked $src -> $dst"
}

for d in "$@"; do migrate "$d"; done

echo "Home data usage: $(du -sh "$REPO" 2>/dev/null | cut -f1)"
touch "${ROOT:-$HOME/DiscreteTreeFlows}/.write_test" && rm "${ROOT:-$HOME/DiscreteTreeFlows}/.write_test" \
  && echo "write test OK" || echo "WARNING: home still read-only"
