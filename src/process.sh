#!/bin/bash
FILES=/pdb/*
for f in *.pdb; do
    # do some stuff here with "$f"
    # remember to quote it or spaces may misbehave
    echo "$f"
    ./test.exe "$f"
done
