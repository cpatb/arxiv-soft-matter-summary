#!/bin/bash
# Copy the Windows tree into ~/pj, strip CRLF, source env, run the given command.
# Windows->WSL line endings: Python's write on Windows emits CRLF, and a shell script
# with CRLF dies in bash on `do\r`. Strip on every sync.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj
cp -r "$SRC"/* ~/pj/ 2>/dev/null
find ~/pj -name "*.py" -o -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
exec "$@"
