#!/bin/bash
FILES=/pdb/*
for f in *.pdb; do
    echo "$f"
    ./test.exe "$f"
done
