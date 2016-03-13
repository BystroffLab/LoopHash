#include "protein.h"

//------------------Constructor------------------//
Protein::Protein(const char* filename){
  //Initialize hardcoded alanine for adding beta carbon to glycine
  //Taken from: https://www.nyu.edu/pages/mathmol/library/life/life1.html
  std::vector< std::vector<float> > alanine;
  alanine.resize( 5, std::vector<float>(3, 0.0) );
  alanine[0][0] = 0.039; alanine[0][1] = -0.028; alanine[0][2] = 0.000;  //N
  alanine[1][0] = 1.499; alanine[1][1] = -0.043; alanine[1][2] = 0.000;  //CA
  alanine[2][0] = 2.055; alanine[2][1] =  1.361; alanine[2][2] = 0.000;  //C
  alanine[3][0] = 1.321; alanine[3][1] =  2.356; alanine[3][2] = 0.011;  //O
  alanine[4][0] = 1.956; alanine[4][1] = -0.866; alanine[4][2] = -1.217; //CB

  //Initialize map of 3-letter amino acid codes to 1 letter codes
  std::map< std::string, char > amino_codes;
  amino_codes["ALA"] = 'A';   amino_codes["CYS"] = 'C';
  amino_codes["ASP"] = 'D';   amino_codes["GLU"] = 'E';
  amino_codes["PHE"] = 'F';   amino_codes["GLY"] = 'G';
  amino_codes["HIS"] = 'H';   amino_codes["ILE"] = 'I';
  amino_codes["LYS"] = 'K';   amino_codes["LEU"] = 'L';
  amino_codes["MET"] = 'M';   amino_codes["ASN"] = 'N';
  amino_codes["PRO"] = 'P';   amino_codes["GLN"] = 'Q';
  amino_codes["ARG"] = 'R';   amino_codes["SER"] = 'S';
  amino_codes["THR"] = 'T';   amino_codes["VAL"] = 'V';
  amino_codes["TRP"] = 'W';   amino_codes["TYR"] = 'Y';



  //Try to open the PDB file for reading
  std::ifstream pdb_str(filename);
  if (!pdb_str.good()) {
    std::cerr << "Couldn't open " << filename << std::endl;
    throw 0;
  }

  //Check to make sure file isn't empty
  if ( pdb_str.peek() == std::ifstream::traits_type::eof() ){
    std::cerr << "File " << filename << " is empty!" << std::endl;
    throw 0;
  }

  //EACH LINE IS 80 CHARACTERS
  //Get 4 letter identifier
  std::string line;
  std::string token(6 , ' ');
  std::string AA(3 , ' ');
  std::string chainID;
  int residue_number;
  std::vector<float> xyz(3);


  //In the PDB Select data given, records start with ATOM, identifier is just the file name minus .pdb
  /*
  std::getline(pdb_str, line);
  identifier = line.substr( line.size() - 4 , 4 );
  */
  std::string tmp_identifier = filename;
  identifier = tmp_identifier.substr(0, 4);
  std::cout << "IDENTIFIER: " << identifier << " Length: " << identifier.size() << std::endl;

  //Skip ahead to the ATOM records (if necessary)
  while (token != "ATOM  "){
    std::getline(pdb_str, line);
    token = line.substr(0,6);   //Record identifiers are 6 characters long
  }

  //If the first atom is CA, assume this file is only CAs
  std::string ca_test = line.substr(12, 3);
  if (ca_test == " CA"){
    std::cout << "this is getting tripped" << std::endl << std::endl;
    throw 0;
  }

  while ( !pdb_str.eof() ){


    if (token == "CONECT"){ break; }
    //We don't care about HETATMs or TERs, skip to next line
    else if (token == "HETATM" || token == "TER   "){
      std::getline(pdb_str, line);
      token = line.substr(0,6);
    }

    else if (token == "ATOM  "){
      //std::cout << line << std::endl;

      //Some terminal atoms end in OXT, ignore these
      std::string oxt_test = line.substr(12, 3);
      if (oxt_test == "OXT"){ break; }

      //We want coordinates for backbone atoms N, CA, C, O, CB
      AA = line.substr(17,3);
      residue_number = atoi( line.substr(23,3).c_str() );
      chainID = line.substr(21,1);


      //Glycine doesn't have a beta carbon, superpose an alanine onto it to get CB coordinates
      if (AA == "GLY"){

        //Create temporary vector to hold glycine coordinates
        std::vector< std::vector<float> > tmp_gly;
        for (int i = 0; i < 4; ++i){
          xyz[0] =  atof( line.substr(30,8).c_str() );
          xyz[1] =  atof( line.substr(38,8).c_str() );
          xyz[2] =  atof( line.substr(46,8).c_str() );
          tmp_gly.push_back(xyz);
          //Get next line
          std::getline(pdb_str, line);
          if ( pdb_str.eof() ){ break; }

          //Some terminal atoms end in OXT, ignore these
          std::string oxt_test = line.substr(12, 3);
          if (oxt_test == "OXT"){ break; }

          token = line.substr(0,6);
        }


        /*Superpose (N, CA, C, O) of Alanine onto Glycine to get rotation matrix,
        translation vector. Apply these to Ala-CB and use coordinates as Gly-CB.
        */
        matrix tmp_ala(alanine.begin(), alanine.end() - 1);
        std::vector<float> new_CB = alanine[4];
        return_val transformation;

        if (tmp_gly.size() == 4){
          try{
            transformation = superimposer(tmp_gly, tmp_ala, 4);
            superimposer_move(new_CB, transformation.first, transformation.second);
            tmp_gly.push_back(new_CB);
          }
          catch(int){
            std::cerr << "Superimposer failed at resiude " << residue_number << std::endl;
            throw 0;
          }
        }


        //Add glycine to vector of peptide backbone_coordinates if all atoms are there
        if (tmp_gly.size() == 5){
          for (int i = 0; i < 5; ++i){
            residue_type.push_back( amino_codes.find(AA)->second );
            chain.push_back( chainID[0] );
            atom_residue_numbering.push_back(residue_number);
            backbone_coordinates.push_back(tmp_gly[i]);
          }
        }


      } //End of residue

      else{
        //Grab this line and the next four atoms' xyz coordinates, building point vector and pushing info into vectors
        for (int i = 0; i < 5; ++i){

          xyz[0] =  atof( line.substr(30,8).c_str() );
          xyz[1] =  atof( line.substr(38,8).c_str() );
          xyz[2] =  atof( line.substr(46,8).c_str() );

          residue_type.push_back( amino_codes.find(AA)->second );
          chainID = line.substr(21,1);
          chain.push_back( chainID[0] );
          atom_residue_numbering.push_back(residue_number);
          backbone_coordinates.push_back(xyz);

          //Get next line
          std::getline(pdb_str, line);
          if ( pdb_str.eof() ){ break; }
          //std::cout << line << std::endl;
          token = line.substr(0,6);

        }

        //Skip past the rest of the associated atoms in the residue
        if ( pdb_str.eof() ){ break; }
        int tmp_residue_number = atoi( line.substr(23,3).c_str() );
        while (tmp_residue_number == residue_number && (token != "HETATM" && token != "TER   ") ){
          std::getline(pdb_str, line);
          if ( pdb_str.eof() ){ break; }
          token = line.substr(0,6);
          tmp_residue_number = atoi( line.substr(23,3).c_str() );
        }
      } //End of Residue
    } //End of chain
  } //End of ATOM record chain

  pdb_str.close();

  /*Sometimes there's an extra atom or two from a residue not fully included in the chain.
    Pop these off the back of the vector until only complete residues remain
  */
  while (backbone_coordinates.size() % 5 != 0){
    backbone_coordinates.pop_back();
    atom_residue_numbering.pop_back();
    chain.pop_back();
    residue_type.pop_back();
  }

}//End Protein constructor




//------------------Loop Generation------------------//

void Protein::getLoops(int length, std::vector< std::vector< std::vector<float> > >& loop_list){
  //Use a sliding window of the length given to generate loops and push into vector
  if (length == 0 || length == 1){
    std::cerr << "ERROR: Loop of length less than 2 residues requested." << std::endl;
    return;
  }
  if (backbone_coordinates.size() == 0){
    std::cerr << "ERROR: No backbone coordinates." << std::endl;
    return;
  }

  for(unsigned int i = 0; i < (backbone_coordinates.size() / 5) - length -1 ; ++i){
      std::vector< std::vector<float> > residues(backbone_coordinates.begin() + (i*5), backbone_coordinates.begin() + (i*5) + (length * 5));
      loop_list.push_back(residues);
  }

  return;
}//End getLoops




//------------------Distance Calculators------------------//

float ca_ca_dist(std::vector< std::vector<float> > loop){
  if (loop.size() == 0){
    std::cerr << "ERROR: No loop to measure. " << std::endl;
    return 0;
  }
  //CA1 = loop[1];
  //CA2 = loop[loop.size() - 5];
  float x = pow( (loop[1][0] - loop[loop.size() - 4][0]) , 2);
  float y = pow( (loop[1][1] - loop[loop.size() - 4][1]) , 2);
  float z = pow( (loop[1][2] - loop[loop.size() - 4][2]) , 2);

  return sqrt( x + y + z);
}

float cb_cb_dist(std::vector< std::vector<float> > loop){
  if (loop.size() == 0){
    std::cerr << "ERROR: No loop to measure. " << std::endl;
    return 0;
  }
  //CB1 = loop[4];
  //CB2 = loop[loop.size() - 1];
  float x =  pow( (loop[4][0] - loop[loop.size() - 1][0]) , 2);
  float y =  pow( (loop[4][1] - loop[loop.size() - 1][1]) , 2);
  float z =  pow( (loop[4][2] - loop[loop.size() - 1][2]) , 2);

  return sqrt( x + y + z);
}




//------------------Statistics------------------//

float RMSD(const std::vector< std::vector<float> >& loop1, const std::vector< std::vector<float> >& loop2){
  //Loops (or molecules) need to be of the same size
  if (loop1.size() != loop2.size()){
    std::cerr << "Loops not of the same size!" << std::endl;
    throw 0;
  }

  //Make working copies, Initialize sum and n
  std::vector< std::vector<float> > loop1_copy = loop1;
  std::vector< std::vector<float> > loop2_copy = loop2;
  unsigned int n = loop1.size();
  float sum;

  //Center molecules
  std::vector<float> loop1_center(3, 0.0);
  std::vector<float> loop2_center(3, 0.0);

  for (unsigned int i = 0; i < n; ++i){
    for (int j = 0; j < 3; ++j){
      loop1_center[j] = loop1_center[j] + loop1[i][j];
      loop2_center[j] = loop2_center[j] + loop2[i][j];
    }
  }

  for (int i = 0; i < 3; ++i){
    loop1_center[i] = loop1_center[i] / float( loop1.size() );
    loop2_center[i] = loop2_center[i] / float( loop1.size() );
  }

  for (unsigned int i = 0; i < n; ++i){
    for (int j = 0; j < 3; ++j){
      loop1_copy[i][j] = loop1[i][j] - loop1_center[j];
      loop2_copy[i][j] = loop2[i][j] - loop2_center[j];
    }
  }

  //Loop through atoms
  for (unsigned int i = 0; i < n; ++i){
    //Loop through xyz coordinates of each atom, calculate difference, square, and add so sum
    for (int j = 0; j < 3; ++j){
      sum += pow( loop2_copy[i][j] - loop1_copy[i][j] , 2);
    }
  }

  return sqrt( sum / n );
}


float standard_deviation(const std::vector<float>& values, float mean){
  float sum = 0;

  for (unsigned int i = 0; i < values.size(); ++i){
    sum = sum + pow(mean - values[i] , 2);
  }

  return sqrt( sum / values.size() );
}



//------------------Superimposer------------------//
std::pair< std::vector<std::vector<float> > , std::vector<float> > superimposer(std::vector< std::vector<float> > coord0, std::vector< std::vector<float> > coord1, unsigned int natm){
  typedef std::vector< std::vector<float> > matrix;
  typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;
  //coord0 is the target (probably the input glycine)
  //coord1 is the prediction (alanine backbone minus CB, superpose onto glycine then apply transformation to CB)
  //natm is the number of atoms

  //Given code takes a matrix with the opposite dimensions used, so tranpose 90deg clockwise if necessary
  if (coord0[0].size() == 3){
    coord0 = transpose(coord0);
    coord1 = transpose(coord1);
  }


  matrix mtx;
  mtx.resize(3, std::vector<float>(3, 0.0) );
  std::vector<float> t(3, 0.0);
  std::vector<float> vec(3, 0.0);
  float tolerance = 0.0001;
  float err = 0.0;

  //Error checking
  unsigned int MAXATOM = 200;
  if (natm > MAXATOM){
    std::cerr << "Error! Increase maxatom" << std::endl;
    throw 0;
  }
  if (natm < 4){
    std::cerr << "Error! Too few atoms (must be 4 or greater)" << std::endl;
    throw 0;
  }
  for (unsigned int i = 0; i < natm; ++i){
    if (std::max(coord0[0][i], coord1[0][i]) > 998){
      std::cerr << "Coordinate error!" << coord0[0][i] << " " << coord1[0][i] << std::endl;
      throw 0;
    }
  }

  //---Center on origin---//
  std::vector<float> x0_center(3, 0.0);
  std::vector<float> x1_center(3, 0.0);
  //For each axis, sum all atom coordinates, divide by natm to get avg d to origin
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      x0_center[j] = x0_center[j] + coord0[j][i];
      x1_center[j] = x1_center[j] + coord1[j][i];
    }
  }
  for (int i = 0; i < 3; ++i){
    x0_center[i] = x0_center[i] / float(natm);
    x1_center[i] = x1_center[i] / float(natm);
  }

  //Initialize matrices to hold new centered coordinates
  matrix x0;
  x0.resize(3, std::vector<float>(natm, 0.0));
  matrix x1;
  x1.resize(3, std::vector<float>(natm, 0.0));


  //Subtract avg distance from origin from actual coordinates, fill new matrices
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      x0[j][i] = coord0[j][i] - x0_center[j];
      x1[j][i] = coord1[j][i] - x1_center[j];
    }
  }


  //---End Center on origin---//


  //Not sure what this is but it seems like a bunch of dot products
  matrix aa;
  aa.resize(3, std::vector<float>(3, 0.0) );
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      for (int k = 0; k < 3; ++k){
        aa[j][k] = aa[j][k] + (x1[j][i] * x0[k][i]);
      }
    }
  }


  //Initialize rotation matrix
  matrix rotation;
  rotation.resize(3, std::vector<float>(3, 0));
  rotation[0][0] = 1.0; rotation[1][1] = 1.0; rotation[2][2] = 1.0;


  //---Iterative rotation scheme---//
  int iteration_count = 0;
  bool do51 = true; //"This is a way to deal with those nasty gotos in the FORTRAN code"
  int iflag;
  int ix;
  int iy;
  int iz;
  double sigma;
  double gamma;
  double sg;
  double bb;
  double cc;

  while (true){
    if (do51){
      iflag = 0;
      ix = 0;
    }

    //If the number of iterations exceeds 1000, give up
    ++iteration_count;
    if (iteration_count > 500){ break; }

    //ix,iy,iz are only allowed to be 0, 1, or 2?
    iy = ix + 1;
    if (iy == 3){
      iy = 0;
    }
    iz = 3 - ix - iy;

    //What are sigma and gamma?
    sigma = aa[iz][iy] - aa[iy][iz];
    gamma = aa[iy][iy] + aa[iz][iz];
    sg = sqrt( sigma * sigma + gamma * gamma );


    //If root of (sigma^2 + gamma^2) = 0. this is tripped
    if (sg == 0){
      //"Goto 50 in FORTRAN code"
      //If do51 was true at the beginning of this iteration, this will break loop
      ++ix;
      if (iflag == 0){ break; }
      if (ix < 3){ do51 = false; }
      else{ do51 = true; }
      continue;
    }

    sg = 1.0 / sg;
    if (fabs(sigma) < ( tolerance * fabs(gamma)) ){
      //"Goto 50 in FORTRAN code"
      ++ix;
      if (iflag == 0){ break; }
      if (ix < 3){ do51 = false; }
      else{ do51 = true; }

      continue;
    }

    //Actual rotations (??)
    for (int i = 0; i < 3; ++i){
      bb = (gamma * aa[iy][i]) + (sigma * aa[iz][i]);
      cc = (gamma * aa[iz][i]) - (sigma * aa[iy][i]);
      aa[iy][i] = bb * sg;
      aa[iz][i] = cc * sg;
      bb = gamma * rotation[iy][i] + sigma * rotation[iz][i];
      cc = gamma * rotation[iz][i] - sigma * rotation[iy][i];
      rotation[iy][i] = bb * sg;
      rotation[iz][i] = cc * sg;
    }
    iflag = 1;

    //Goto 50 in FORTRAN
    ++ix;
    if (iflag == 0){ break; }
    if (ix < 3){ do51 = false; }
    else{ do51 = true; }

    continue;


  } //End while loop


  //Change x1 coordinates

  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      t[j] = 0.0;
      for (int k = 0; k < 3; ++k){
        //Seemingly a rotation
        t[j] = t[j] + (rotation[j][k] * x1[k][i]);
        //For some reason everything's coming out rotated 180deg about x(?), so let's rotate it back
      }
    }
    for (int j = 0; j < 3; ++j){
      x1[j][i] = t[j];
    }
  }


  //Calculate RMSD
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; i < 3; ++i){
      err = err + pow( (x0[j][i] - x1[j][i]) , 2 );
    }
  }
  err = sqrt(err / natm);

  //"un-center" points
  for (unsigned int i = 0; i < natm; ++i){
    for (int j = 0; j < 3; ++j){
      coord1[j][i] = x1[j][i] + x0_center[j];
    }
  }

  //Build output matrices
  for (int i = 0; i < 3; ++i){
    t[i] = x0_center[i];
    for (int j = 0; j < 3; ++j){
      mtx[i][j] = rotation[i][j];
      t[i] = t[i] - rotation[i][j] * x1_center[j];
    }
    vec[i] = t[i];
  }


  //std::cout << "RMSD: " << err << std::endl << std::endl;

  return_val ret = std::make_pair(mtx, vec);
  return ret;
}

void superimposer_move(std::vector<float>& x, std::vector< std::vector<float> >& mtx, std::vector<float>& vec){
  std::vector<float> y(3, 0.0);
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      y[i] = y[i] + (x[j] * mtx[i][j]);
    }
  }
  for (int i = 0; i < 3; ++i){
    y[i] = y[i] + vec[i];
  }

  x = y;
  return;
}


std::vector< std::vector<float> > transpose(std::vector< std::vector<float> >& coord0){
  std::vector< std::vector<float> > tmp;
  for (unsigned int i = 0; i < coord0[0].size(); ++i){
    std::vector<float> xyz;
    for (unsigned int j = 0; j < coord0.size(); ++j){
      xyz.push_back(coord0[j][i]);
    }
    tmp.push_back(xyz);
    xyz.clear();
  }
  return tmp;
}


//Matrix output for debugging
void printVec(std::vector< std::vector<float> >& vec){
  for (unsigned int i = 0; i < vec.size(); ++i){
    for (unsigned int j = 0; j < (vec[i]).size(); ++j){
      std::cout << vec[i][j] << " ";
    }
    std::cout << std::endl;
  }
}




//------------------Hashing------------------//

unsigned long int loopHash(std::vector< std::vector<float> > loop){
  float CA_CA = ca_ca_dist(loop);
  float CB_CB = cb_cb_dist(loop);
  int loop_length = loop.size() / 5;
  int key;

  //truncate floats
  char sz[64];
  double lf = CA_CA;
  sprintf(sz, "%.1lf\n", lf);
  double lf2 = atof(sz);
  CA_CA = lf2;

  lf = CB_CB;
  sprintf(sz, "%.1lf\n", lf);
  lf2 = atof(sz);
  CB_CB = lf2;

  //Combine distances in a order-specific way
  //Decrease coefficient to weight data point higher
  float x = 0.5; float y = 0.8; float z = 0.5;

  int sum = int( (CA_CA / x) + (CB_CB / y) + (loop_length / z) );
  //std::hash<float> test;

  //unsigned long int hash_value = test(sum);
  int tablesize = 513;
  key = sum % tablesize ;

  return key;

}



//------------------Output ------------------//
void PDB_out (const std::vector< std::vector<float> >& loop, char* filename){
  //Try to open file for writing
  std::ofstream out_file(filename);
  if ( !out_file.good() ) {
    std::cerr << "Can't open " << filename << " to write." << std::endl;
    exit(1);
  }

  /*
  -----PDB ATOM RECORD-----
  col  |  field
  -------------------------------------------
  0-5    Record name (ATOM  )
  6-10   Serial number (integer)
  12-15  Atom name
  16     Alternate location indicator
  17-19  Residue name
  21     Chain identifier (char)
  22-25  Residue sequence number
  26     Code for insertion of residues
  30-37  Orthogonal x coordinates (angstroms)
  38-45  Orthogonal y coordinates
  46-53  Orthogonal z coordinates
  54-59  Occupancy (float)
  60-65  B factor
  76-77  Element symbol (right-justified)
  78-79  Atom charge
  */

  std::vector<std::string> atoms = {" N  ", " CA ", " C  ", " O  ", " CB "};
  std::vector<std::string> elements = {"N ", "C ", "C ", "O ", "C "};
  int counter1 = 0;
  int counter2 = 1;
  std::string whitespace(1, ' ');

  for (unsigned int i = 0; i < loop.size(); ++i, ++counter1){
    //Non-numeric data
    std::string line(80, ' ');
    line.replace(0, 6, "ATOM  ");
    line.replace(12, 4, atoms[i%5]);
    line.replace(17, 3, "ALA");
    line.replace(21, 1, "A");

    //xyz coordinates
    std::stringstream xyz;
    xyz.precision(3);
    xyz << std::setw(6) << loop[i][0] << "  " << std::setw(6) << loop[i][1] << "  " << std::setw(6) << loop[i][2];
    std::string xyz_str = xyz.str();
    line.replace(32, 20, xyz_str);

    //Atom counting
    int sn = i + 1;
    //Residue counting
    if (counter1 == 5){
      counter1 = 0;
      ++counter2;
    }

    //No convenient way to convert ints to strings for serial number, residue number
    std::stringstream ss1;
    std::stringstream ss2;
    ss1 << sn;
    ss2 << counter2;

    std::string serial_no = ss1.str();
    std::string residue_no = ss2.str();

    whitespace.resize(5 - serial_no.size(), ' ');
    serial_no =  whitespace+ serial_no;

    whitespace.resize(4 - residue_no.size(), ' ');
    residue_no = whitespace + residue_no;
    line.replace(6, 5, serial_no);
    line.replace(22, 4, residue_no);
    line.replace(75, 2, elements[i%5]);

    std::cout << line << std::endl;
    out_file << line << std::endl;
  }//End of ATOM records

  //Terminate Chain
  std::string line(80, ' ');
  line.replace(0, 6, "TER   ");
  out_file << line << std::endl;

  //End of file
  line.replace(0, 6, "END   ");
  std::cout << line << std::endl;
  out_file << line << std::endl;

out_file.close();
}


//------------------Master protein random access file output ------------------//

bool Protein::RAF_out(char* filename){
  /*            FILE STRUCTURE
    Every residue record is 70 bytes:
    4 letter pdb code (4 bytes) - 1 letter aa code (1 byte) - int residue # (4 bytes) - char chain id (1 byte) - 10 bytes total
    Coordinates for N, CA, C, O, CB 12 bytes per atom (4 bytes per float) - 60 total
    This function just appends to the end of the file
  */

  //Try to open file
  std::ofstream iofile(filename, std::fstream::binary | std::fstream::app);
  if ( !iofile ) {
    std::cerr << "Can't open " << filename << " to write." << std::endl;
    return false;
  }


  //Write info
  for (unsigned int i = 0; i < backbone_coordinates.size() ; i += 5){
    iofile.write( identifier.c_str() , sizeof(char) * 4 );                 //4 Letter PDB identifier
    iofile.write( (const char*)& residue_type[i] , sizeof(char) );           //1 letter amino acid code
    iofile.write( (const char*)& atom_residue_numbering[i] , sizeof(int) );  //Residue #
    iofile.write( (const char*)& chain[i], sizeof(char) );                   //Chain ID

    //Write coordinates
    for (unsigned int j = i; j < i + 5; ++j){
      //std::cout << backbone_coordinates[j][0] << " " << backbone_coordinates[j][1] << " " <<backbone_coordinates[j][2] << std::endl;
      iofile.write( (const char*)& backbone_coordinates[j][0], sizeof(float) );
      iofile.write( (const char*)& backbone_coordinates[j][1], sizeof(float) );
      iofile.write( (const char*)& backbone_coordinates[j][2], sizeof(float) );
    }
  }

  iofile.close();
  return true;
}
