#!/bin/bash
tmpfile="$1.clean_warnings.bash.tmp"
grep -v '^ *"WARNING:pystan.*",' "$1" > "$tmpfile"
mv "$tmpfile" "$1"
rm -f "$tmpfile"
