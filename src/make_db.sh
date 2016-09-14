#!/bin/bash
# Constructs master list of pdb files. Place in directory of pdb files, as well
# as pdblist_creation.exe

for f in *.pdb
do
    echo $f;
    ./pdblist_creation.exe "$f" pdblist.dat;
done
./database_creation.exe pdblist.dat looplist.dat grid.dat
