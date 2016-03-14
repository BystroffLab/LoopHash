#ifndef loop_generation_h
#define loop_generation_h

#include "protein.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstring>
#include <utility>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

typedef std::vector<std::vector<std::vector<int> > > vec_3D;

/*
  Generates loops given a random access file of the following format:

  Every residue record is 70 bytes:
  4 letter pdb code (4 bytes) - 1 letter aa code (1 byte) - int residue # (4 bytes) - char chain id (1 byte) - 10 bytes total
  Coordinates for N, CA, C, O, CB 12 bytes per atom (4 bytes per float) - 60 total

  Writes to a second random access file. Each record is X bytes long, containing
  numerical pointers to the first and last record numbers of the loop. Each
  entry in a record should be 8 bytes (2 ints)

  ------------------------------------------------------------------------------
  countLoops

  Loops should only have 1 parent chain, unless the 2 chains are contiguous and
  have CAs than 4 angstroms apart.

  First pass - 3D vector of int pairs with indexes using int multiples of CA-CA,CB-CB, and
  loop length. pair.first = number of total loops counted | pair.second = number of loops actually inserted


  writeLoops

  Takes random access file, 3D array of int triples. Contructs random access file by looping
  through 3D array, getting record location by keeping track of all previous records generated
  and multiplying loops x residue length x 8 bytes. Keeps track of how filled the record is



*/

void test_raf_in (char* filename);

void countLoops(char* filename, int length, vec_3D& grid);

void writeLoops(char* pdbselect, char* loopfile, char* gridfile, const vec_3D& grid);

void readLoops(char* loopfile);

void readGrid(char* gridfile);

void pruneDB(char* pdbselect, char* loopfile, char* gridfile);

void db_query(float dCA, float dCB, int loop_length, char* pdbselect, char* loopfile, char* gridfile);



#endif
