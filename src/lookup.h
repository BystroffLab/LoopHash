//Lookup object header file
#ifndef lookup_h
#define lookup_h

#include "superimposer.h"
#include "protein.h"
#include "complex.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>
#include <functional>
#include <list>
#include <string>
#include <string.h>
#include <vector>



struct Loop{
  std::vector<std::vector<float> > coordinates;
  std::vector<char> sequence;
  float rmsd;
};



class Lookup{
public:

  // Constructor
  Lookup(char* input_file);
  void parse(char* input_file);

  // Getters
  int size() const{ return results.size(); }
  std::list<Loop> getResults() const{ return results; }
  std::vector<std::vector<float> > getOriginal() const{ return original_loop; }

  // Setters
  void setMin(int x);
  void setMax(int x){ max_results = x; }
  void setCutoff(float x){ duplicate_threshold = x; }
  void setCollisions(float x){ collision_cutoff = x*x; logmsg("Set collision cutoff to " + std::to_string(collision_cutoff) + "\n");}
  void setSequence(std::string s, float identity);
  void setRange(int min_length, int max_length);

  // Output
  void writeLog();
  void logmsg(std::string msg){ logdump.push_back(msg); }
  void iRosettaOutput();

  // Lookup job
  void run();


private:

  // Containers and database files
  std::list<Loop> results;
  std::list<Loop> results_buffer;
  std::vector<Protein> scaffold;
  Complex complex;
  std::vector<std::vector<float> > original_loop;
  std::vector<std::vector<float> > original_loop_anchors;
  std::vector<char*> database_files; // [0]=pdb select, [1]=loop db, [2]=grid
  char* logfile;
  std::vector<std::string> logdump;

  // Parameters
  int scaffold_start;               // N term anchor
  int scaffold_end;                 // C term anchor
  int min_results;                  // Minimum number of loops return
  int max_results;                  // Max number of loops to return
  int length_range[2];              // Range of loop lengths
  float rmsd_cutoff;                // Highest anchor RMSD to accept
  float collision_cutoff;           // Value below which is considered a steric collision (A^2)
  float duplicate_threshold;        // RMSD threshold for duplicates
  bool filter;                      // Are we filtering by sequence?
  std::string sequence_filter;      // Sequence filter
  float sequence_identity_cutoff;   // Minimum identity
  int symmetry;                     // How many monomers in complex? (NOT IMPLEMENTED)
  bool preserve_sequence;

  // Statistics to report at the end of the search
  unsigned int database_hits;
  unsigned int scaffold_colliding_loops;
  unsigned int complex_colliding_loops;
  unsigned int loop_colliding_loops;
  unsigned int redundant_loops;
  unsigned int bad_fits;

  // Private member functions
  void runHelper(float CA_CA, float CB_CB, int loop_length);
  bool isDuplicate(const Loop &candidate);
  bool isCollision(const Loop &loop1, const Loop &loop2);
  void updateBuffer();
  void cleanDuplicates();
  void cleanBadFits(std::list<Loop>& results);
  void cleanCollisions(std::list<Loop>& results);
  std::vector<std::vector<float> > collectAnchors(const std::vector<std::vector<float> > &loop);
  void superimposeUsingAnchors(Loop &database_loop, const std::vector<std::vector<float> > &original_loop_anchors);
  std::vector<std::vector<Loop> > superimposeForSymmetry();


};

#endif
