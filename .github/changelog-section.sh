#!/bin/sh
# Print the CHANGELOG.md section for one version: the lines between its
# "## X.Y.Z ..." heading and the next "## " heading, with trailing blank lines
# dropped. Exits 1 if the version has no section, or the section is empty.
#
# Usage: changelog-section.sh X.Y.Z [CHANGELOG.md]
set -eu
VERSION=$1
FILE=${2:-CHANGELOG.md}

awk -v v="$VERSION" '
  /^## / { if (found) exit; if ($2 == v) { found = 1; next } }
  found  { lines[++n] = $0; if ($0 ~ /[^[:space:]]/) last = n }
  END    {
    if (!found || !last) exit 1
    for (i = 1; i <= last; i++) print lines[i]
  }
' "$FILE" || {
  echo "CHANGELOG: no non-empty '## $VERSION' section in $FILE" >&2
  exit 1
}
