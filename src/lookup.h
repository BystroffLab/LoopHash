//Lookup object header file
#ifndef lookup_h
#define lookup_h

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>
#include <functional>
#include <list>
#include <string>
#include <vector>
#include "superimposer.h"
#include "protein.h"

//Loop struct for more convenience
struct Loop{
  std::vector<std::vector<float> > coordinates;
  std::vector<char> sequence;
  float rmsd;
};

class Lookup{
public:

  //Constructors
  Lookup();
  Lookup(Protein* protein, int start, int end);

  //Getters
  int size() const{ return results.size(); }
  std::list<Loop> getResults() const{ return results; }
  std::vector<std::vector<float> > getOriginal() const{ return original_loop; }

  //Setters
  void setMin(int x){ min_results = x; }
  void setMax(int x){ max_results = x; }
  void setCutoff(float x){ rmsd_cutoff = x; }
  void setDB(char* pdb, char* loops, char* grid);
  void setSequence(std::string s, float identity);
  void setRange(int min_length, int max_length);

  //Miscellaneous
  void run();

private:
  //Lookup results, original loop, database files
  std::list<Loop> results;
  Protein* scaffold;
  std::vector<std::vector<float> > original_loop;
  std::vector<char*> database_files; // [0]=pdb select, [1]=loop db, [2]=grid

  //Various parameters for lookup
  int scaffold_start;
  int scaffold_end;
  int min_results;                          //Minimum number of loops to try
  int max_results;                          //Max number of loops to try
  int length_range[2];              //Range of lengths to try
  float rmsd_cutoff;                //Highest RMSD to accept
  bool filter;                      //Are we filtering by sequence?
  std::string sequence_filter;      //Sequence filter
  float sequence_identity_cutoff;   //Minimum identity

  //Run helper
  void runHelper(float CA_CA, float CB_CB, int loop_length);

};

#endif
