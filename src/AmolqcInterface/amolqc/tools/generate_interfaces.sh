#!/bin/bash

set -e

dir=$(dirname "$0")
f2py="$dir/f2py/f2py"

if [[ ! -f "$f2py" ]]; then
  echo "f2py not found. Expected in $f2py" >&2
  exit 1
fi

if [[ $# = 0 ]]; then
  echo "Need output module name" >&2
  exit 1
fi
outmod="$1"
out="$outmod.f90"
shift

if [[ $# = 0 ]]; then
  set .
fi

> "$out"
tmp=$(mktemp)

echo "! Interfaces automatically generated with generate_interfaces.sh" >> "$out"

for target; do
  while IFS= read -d $'\0' file; do
    echo -n "Checking $file... "
    module=$(grep -m1 -iR "^\s*module" "$file" || true) # TODO remove
    if [[ -z "$module" ]]; then
      echo "no module found"
      continue
    fi

    #echo "! file $file" >> "$out" # XXX ausgeben wenn module <name> entfernt wird
    python "$f2py" --quiet --no-lower --overwrite-signature -h "$tmp" "$file" > /dev/null

    #grep -v "\(end \)\?python module" "$tmp" >> "$out"
    cat "$tmp" >> "$out"
    echo "done"

  done < <(find "$target" -iname "*.f" -or -iname "*.f90" -print0)
done
