#include "complex.h"



Complex::Complex(){
  return;
}



Complex::Complex(char* input_pdb)
{
  parse(input_pdb);
  return;
}


/*
    Adds a molecule's coordinates to the object
*/
void Complex::addMolecule(char* input_pdb)
{
  parse(input_pdb);
  return;
}



/*
    Checks to see if a loop collides with the coordinates in the complex
*/
bool Complex::isCollision(const std::vector< std::vector<float> > &loop)
{
  if (coordinates.size() == 0){
    return false;
  }

  for (unsigned int i = 0; i < loop.size(); ++i){
    for (unsigned int j = 0; j < coordinates.size(); ++j){
      if (atom_dist_fast(loop[i], coordinates[j]) < 16.0){
        return true;
      }
    }
  }

  return false;
}



/*
    Extracts atom coordinates from a pdb file
*/
void Complex::parse(char* input_pdb){
  std::ifstream in(input_pdb);

  if (!in){
    std::string err(input_pdb);
    throw std::runtime_error("Couldn't parse Complex file " + err + "\n");
    return;
  }

  // Skip to atoms
  std::string line = "";
  std::getline(in, line);
  while (line.substr(0,6) != "ATOM  "){
    std::getline(in, line);
  }

  // Pull first atom
  std::vector<float> xyz;
  xyz.push_back(atof( line.substr(30,8).c_str()));
  xyz.push_back(atof( line.substr(38,8).c_str()));
  xyz.push_back(atof( line.substr(46,8).c_str()));
  coordinates.push_back(xyz);
  xyz.clear();

  // Pull the rest of the atoms
  while (!in.eof()){
    std::getline(in, line);
    if (line.substr(0,6) == "ATOM  " || line.substr(0,6) == "HETATM"){
      xyz.push_back(atof( line.substr(30,8).c_str()));
      xyz.push_back(atof( line.substr(38,8).c_str()));
      xyz.push_back(atof( line.substr(46,8).c_str()));
      coordinates.push_back(xyz);
      xyz.clear();
    }
  }

  return;
}



/*
    Debug
*/
void Complex::printComplex()
{
  std::cout << "Coordinates: " << std::endl;
  for (unsigned int i = 0; i < coordinates.size(); ++i){
    std::cout << "x: " << coordinates[i][0] << " y: " << coordinates[i][1] << " z: " << coordinates[i][2] << std::endl;
  }
  return;
}
