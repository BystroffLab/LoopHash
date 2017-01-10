#ifndef COMPLEX_H
#define COMPLEX_H

#include "protein.h"
#include <fstream>
#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>


/*
  Simple class to store the coordinates of non-scaffold macromolecule
  atoms in a complex.
*/
class Complex{
public:
  Complex();
  Complex(char* input_pdb);

  unsigned int size(){ return coordinates.size(); }
  void addMolecule(char* input_pdb);
  bool isCollision(const std::vector<std::vector<float> > &loop);
  void printComplex();

private:
  std::vector<std::vector<float> > coordinates;
  void parse(char* input_pdb);
};

#endif
