//---Protein header file---//

#ifndef protein_h
#define protein_h

#include "superimposer.h"
#include <exception>
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

class Protein {
public:
  //Constructor
  Protein();                      //NOTE: This is unsafe and just used when protein class is a member of another
  Protein(const char* filename);
  Protein(const std::string &new_identifier,
          const std::vector<std::vector<float> > &new_backbone_coordinates,
          const std::vector<char> &new_residue_type,
          const std::vector<int> &new_atom_residue_numbering,
          const std::vector<char> &new_chain
        );


  //Getters
  std::vector<std::vector<float> > getCoordinates() const{ return backbone_coordinates; }
  std::vector<std::vector<float> > getLoop(int start, int end) const;
  int size() const{ return backbone_coordinates.size() / 5; }
  std::string getIdentifier() const{ return identifier; }

  bool RAF_out (char* filename);
  bool is_collision (const std::vector<std::vector<float> >& insertion, int start, int end) const;
  const std::vector<int> findChains();
  std::vector<Protein> splitChains();


  //Some typedefs. ***USE BETTER NAMES FOR TYPEDEFS***
  typedef std::vector< std::vector<float> > matrix;
  typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;

private:
  //Representation
  std::string identifier;
  std::vector<std::vector<float> > backbone_coordinates;
  std::vector<char> residue_type;
  std::vector<int> atom_residue_numbering;
  std::vector<char> chain;

  std::map< std::string, char > amino_codes;
  std::vector< std::vector<float> > alanine;

  std::vector<std::string> collectAtoms(const char* filename);
  void addBetaCarbon(std::vector<std::vector<float> > &residue);
  void parseBackbone(const std::vector<std::string> &cleaned_file);

};



//Distance calculators
float ca_ca_dist(const std::vector< std::vector<float> >& loop);
float cb_cb_dist(const std::vector< std::vector<float> >& loop);
float atom_dist(const std::vector<float>& atom1, const std::vector<float>& atom2);
float atom_dist_fast(const std::vector<float>& atom1, const std::vector<float>& atom2);

//Statistics
float RMSD(const std::vector< std::vector<float> >& loop1, const std::vector< std::vector<float> >& loop2);
float standard_deviation(const std::vector<float>& values, float mean);


//Output
void PDB_out (const std::vector< std::vector<float> >& loop, const std::vector<char> residues, char* filename );
void PDB_out (const std::vector< std::vector<float> >& loop, char* filename );


#endif
