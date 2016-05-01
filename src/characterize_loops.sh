#!/bin/bash
FILES=/pdb/*
for f in *.pdb; do
    echo "$f"
    ./characterize.exe "$f" separation_analysis_good.csv 3
done
