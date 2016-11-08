#include "protein.h"

//TODO: add try/catch blocks to substr calls that check token OR test for 'END' as well as 'END   '
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
  std::string tmp_identifier = filename;
  if (tmp_identifier.size() > 3){
    identifier = tmp_identifier.substr(0, 4);
  }
  else{
    identifier = "TEMP";
  }

  //Skip ahead to the ATOM records (if necessary)
  while (token != "ATOM  "){
    std::getline(pdb_str, line);
    token = line.substr(0,6);   //Record identifiers are 6 characters long
  }

  //If the first atom is CA, assume this file is only CAs
  std::string ca_test = line.substr(12, 3);
  if (ca_test == " CA"){
    throw 0;
  }

  while ( !pdb_str.eof() && line.size() > 54 ){
    if (token == "CONECT"){ break; }
    //We don't care about HETATMs or TERs, skip to next line
    else if (token == "HETATM" || token == "TER   "){
      std::getline(pdb_str, line);
      token = line.substr(0,6);
    }

    else if (token == "ATOM  "){
      //Some terminal atoms end in OXT, ignore these
      std::string oxt_test = line.substr(12, 3);
      if (oxt_test == "OXT"){ break; }

      //We want coordinates for backbone atoms N, CA, C, O, CB
      AA = line.substr(17,3);
      residue_number = atoi( line.substr(23,3).c_str() );
      //std::cerr << " Residue: " << residue_number;
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
          token = line.substr(0,6);
        }

        // Skip past the rest of the associated atoms in the residue
        if ( pdb_str.eof() ){ break; }

        // Hopefully this catches the case where the last line is only 3 characters
        int tmp_residue_number = 0;

        try{
          tmp_residue_number = atoi( line.substr(23,3).c_str() );
        }
        catch(std::out_of_range& exception){
          break;
        }

        while (tmp_residue_number == residue_number && (token != "HETATM" && token != "TER   " && token != "END   ") ){
          std::getline(pdb_str, line);
          if ( pdb_str.eof() || line.size() < 6){ break; }
          token = line.substr(0,6);
          tmp_residue_number = atoi( line.substr(23,3).c_str() );
        }
      } //End of Residue
    } //End of chain
  } //End of ATOM record chain

  pdb_str.close();

  /*
    Sometimes there's an extra atom or two from a residue not fully included in the chain.
    Pop these off the back of the vector until only complete residues remain
  */

  while (backbone_coordinates.size() % 5 != 0){
    backbone_coordinates.pop_back();
    atom_residue_numbering.pop_back();
    chain.pop_back();
    residue_type.pop_back();
  }

//std::cout << "Finished" << std::endl;
}//End Protein constructor

//------------------Accessors------------------//
std::vector<std::vector<float> > Protein::getLoop(int start, int end) const{
  if ( start > this->size() || end > this->size() ){
    throw 0;
  }
  std::vector<std::vector<float> > return_loop(backbone_coordinates.begin() + ((start - 1)*5), backbone_coordinates.begin() + (end*5));
  return return_loop;
}


//------------------Measuring------------------//

float ca_ca_dist(std::vector< std::vector<float> > loop){
  if (loop.size() == 0 || loop.size() == 1){
    std::cerr << "ERROR: No loop to measure. " << std::endl;
    return 0;
  }

  //CA1 = loop[1];
  //CA2 = loop[loop.size() - 5];
  float x = pow( (loop[1][0] - loop[loop.size() - 4][0]) , 2);
  float y = pow( (loop[1][1] - loop[loop.size() - 4][1]) , 2);
  float z = pow( (loop[1][2] - loop[loop.size() - 4][2]) , 2);

  return sqrt( x + y + z );
}

float cb_cb_dist(std::vector< std::vector<float> > loop){
  if (loop.size() == 0 || loop.size() == 1){
    std::cerr << "ERROR: No loop to measure. " << std::endl;
    return 0;
  }

  //CB1 = loop[4];
  //CB2 = loop[loop.size() - 1];
  float x =  pow( (loop[4][0] - loop[loop.size() - 1][0]) , 2);
  float y =  pow( (loop[4][1] - loop[loop.size() - 1][1]) , 2);
  float z =  pow( (loop[4][2] - loop[loop.size() - 1][2]) , 2);

  return sqrt( x + y + z );
}

float atom_dist(const std::vector<float>& atom1, const std::vector<float>& atom2){
  float x = pow( atom1[0] - atom2[0] , 2);
  float y = pow( atom1[1] - atom2[1] , 2);
  float z = pow( atom1[2] - atom2[2] , 2);
  return sqrt( x + y + z );
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



//------------------Output ------------------//
//When we actually have the sequence
void PDB_out (const std::vector< std::vector<float> >& loop, const std::vector<char> residues, char* filename){
  ////Initialize map of 3-letter amino acid codes to 1 letter codes
  std::map< char, std::string > amino_codes;
  amino_codes['A'] = "ALA";   amino_codes['C'] = "CYS";
  amino_codes['D'] = "ASP";   amino_codes['E'] = "GLU";
  amino_codes['F'] = "PHE";   amino_codes['G'] = "GLY";
  amino_codes['H'] = "HIS";   amino_codes['I'] = "ILE";
  amino_codes['K'] = "LYS";   amino_codes['L'] = "LEU";
  amino_codes['M'] = "MET";   amino_codes['N'] = "ASN";
  amino_codes['P'] = "PRO";   amino_codes['Q'] = "GLN";
  amino_codes['R'] = "ARG";   amino_codes['S'] = "SER";
  amino_codes['T'] = "THR";   amino_codes['V'] = "VAL";
  amino_codes['W'] = "TRP";   amino_codes['Y'] = "TYR";


  std::ofstream out_file(filename, std::fstream::app);
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
    line.replace(21, 1, " ");

    //xyz coordinates
    std::stringstream xyz;
    xyz.precision(3);
    xyz << std::fixed << std::right << std::setw(7) << loop[i][0] << " " << std::setw(7) << loop[i][1] << " " << std::setw(7) << loop[i][2];
    std::string xyz_str = xyz.str();
    line.replace(31, 20, xyz_str);

    //Atom counting
    int sn = i + 1;
    //Residue counting
    if (counter1 == 5){
      counter1 = 0;
      ++counter2;
    }
    std::string tmp_residue = amino_codes.find(residues[counter2-1])->second;
    line.replace(17, 3, tmp_residue);


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
    line.replace(77, 2, elements[i%5]);

    out_file << line << std::endl;
  }//End of ATOM records

  //Terminate Chain
  std::string line(80, ' ');
  line.replace(0, 6, "TER   ");
  out_file << line << std::endl;

  //End of file
  line.replace(0, 6, "      ");
  //std::cout << line << std::endl;
  out_file << line << std::endl;

out_file.close();
}

/*
  When we don't have the sequence (or don't want to output the original sequence)
  Rosetta automatically builds sidechains when you make a pose from a loop
  Don't need that if you're designing a sequence onto a backbone
*/
void PDB_out (const std::vector< std::vector<float> >& loop, char* filename){

  //Try to open file for writing
  std::ofstream out_file(filename, std::fstream::app);
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
    line.replace(21, 1, " ");

    //xyz coordinates
    std::stringstream xyz;
    xyz.precision(3);
    xyz << std::fixed << std::right << std::setw(7) << loop[i][0] << " " << std::setw(7) << loop[i][1] << " " << std::setw(7) << loop[i][2];
    std::string xyz_str = xyz.str();
    line.replace(31, 20, xyz_str);

    //Atom counting
    int sn = i + 1;
    //Residue counting
    if (counter1 == 5){
      counter1 = 0;
      ++counter2;
    }

    line.replace(17, 3, "ALA");


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
    line.replace(77, 2, elements[i%5]);

    out_file << line << std::endl;
  }//End of ATOM records

  //Terminate Chain
  std::string line(80, ' ');
  line.replace(0, 6, "TER   ");
  out_file << line << std::endl;

  //End of file
  line.replace(0, 6, "      ");
  //std::cout << line << std::endl;
  out_file << line << std::endl;

out_file.close();
}

//------------------        Collision detection         ------------------//
//Detect if any loop atoms (excluding anchors) are within 4A of non-bonded atoms, return true if collision
bool Protein::is_collision (const std::vector<std::vector<float> >& insertion, int start, int end) const{
  bool collision = false;

  //Check all residues before insertion
  //For each backbone atom i before anchor
  for (int i = 0; i < (start * 5) - 10; ++i){
    //For each insertion atom j (excluding anchors)
    for (unsigned int j = 5; j < insertion.size() - 5; ++j){
      if ( atom_dist(insertion[j], backbone_coordinates[i]) < 5.0){
        collision = true;
        //std::cout << "Collision! Atoms are " << atom_dist(insertion[j], backbone_coordinates[i]) << "apart." << std::endl;
        break;
      }
    }
    if (collision == true){ break; }
  }

  //Check all residues after insertion

  for (unsigned int i = (end * 5) + 10; i < backbone_coordinates.size(); ++i){
    //For each insertion atom j (excluding anchors)
    for (unsigned int j = 5; j < insertion.size() - 5; ++j){
      if ( atom_dist(insertion[j], backbone_coordinates[i]) < 3.5){
        collision = true;
        //std::cout << "After-loop collision! Atoms are " << atom_dist(insertion[j], backbone_coordinates[i]) << "apart." << std::endl;
        break;
      }
    }
    if (collision == true){ break; }
  }

  return collision;
}




//------------------PDBselect random access file output ------------------//

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
