//---Protein header file---//

#ifndef protein_h
#define protein_h

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "superimposer.h"

class Protein {
public:
  //Constructor
  Protein(const char* filename);

  //Accessors
  std::vector<std::vector<float> > getCoordinates() const{ return backbone_coordinates; }
  std::string getIdentifier() const{ return identifier; }

  //Output
  bool RAF_out (char* filename);

  typedef std::vector< std::vector<float> > matrix;
  typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;

private:
  //Representation
  std::string identifier;
  int length;
  std::vector<std::vector<float> > backbone_coordinates;
  std::vector<char> residue_type;
  std::vector<int> atom_residue_numbering;
  std::vector<char> chain;

};

//Distance calculators
float ca_ca_dist(std::vector< std::vector<float> > loop);
float cb_cb_dist(std::vector< std::vector<float> > loop);


//Statistics
float RMSD(const std::vector< std::vector<float> >& loop1, const std::vector< std::vector<float> >& loop2);
float standard_deviation(const std::vector<float>& values, float mean);


//Output
void PDB_out (const std::vector< std::vector<float> >& loop, char* filename );


#endif
