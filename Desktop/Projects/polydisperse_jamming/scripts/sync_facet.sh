#!/bin/bash
# copy facet2 raw results (jsonl, census npz, logs) from the WSL run tree to the repo. Raw results only.
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp ~/pj/out/lsdem/ls_facet2*.jsonl ~/pj/out/lsdem/contacts_facet2*.npz ~/pj/out/lsdem/facet2*.log "$SRC/out/lsdem/" 2>/dev/null
wc -l "$SRC/out/lsdem/ls_facet2.jsonl" 2>/dev/null
