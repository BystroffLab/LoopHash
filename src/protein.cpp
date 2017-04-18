#include "protein.h"



/*
    Default constructor
*/
Protein::Protein()
{
  return;
}



/*
    Copy constructor
*/
Protein::Protein(const std::string &new_identifier,
                 const std::vector<std::vector<float> > &new_backbone_coordinates,
                 const std::vector<char> &new_residue_type,
                 const std::vector<int> &new_atom_residue_numbering,
                 const std::vector<char> &new_chain)
{
  identifier = new_identifier;
  backbone_coordinates = new_backbone_coordinates;
  residue_type = new_residue_type;
  atom_residue_numbering = new_atom_residue_numbering;
  chain = new_chain;
  return;
}



/*
    Parses a protein into a backbone representation
*/
Protein::Protein(const char* filename)
{
  //Initialize hardcoded alanine for adding beta carbon to glycine
  //Taken from: https://www.nyu.edu/pages/mathmol/library/life/life1.html
  alanine.resize( 5, std::vector<float>(3, 0.0) );
  alanine[0][0] = 0.039; alanine[0][1] = -0.028; alanine[0][2] = 0.000;  //N
  alanine[1][0] = 1.499; alanine[1][1] = -0.043; alanine[1][2] = 0.000;  //CA
  alanine[2][0] = 2.055; alanine[2][1] =  1.361; alanine[2][2] = 0.000;  //C
  alanine[3][0] = 1.321; alanine[3][1] =  2.356; alanine[3][2] = 0.011;  //O
  alanine[4][0] = 1.956; alanine[4][1] = -0.866; alanine[4][2] = -1.217; //CB

  //Initialize map of 3-letter amino acid codes to 1 letter codes
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

  std::string tmp_identifier = filename;
  if (tmp_identifier.size() > 3){
    identifier = tmp_identifier.substr(0, 4);
  }
  else{
    identifier = "TEMP";
  }

  // Try to read in the pdb file
  try{
    parseBackbone(collectAtoms(filename));
  }
  catch(const std::exception& e){
    if (backbone_coordinates.size() % 5 == 0){
      return;
    }
    else{
      std::cout << "Exception: " << e.what() << std::endl;
      throw std::runtime_error(e.what());
      return;
    }
  }

  return;
}



/*
    Pulls a const loop from a scaffold if anchors exist
*/
std::vector<std::vector<float> > Protein::getLoop(int start, int end) const
{
  if ( start > this->size() || end > this->size() ){
    throw 0;
  }
  std::vector<std::vector<float> > return_loop(backbone_coordinates.begin() + ((start - 1)*5), backbone_coordinates.begin() + (end*5));
  return return_loop;
}



/*
    Measures the distance in angstroms between two alpha carbons of a loop
*/
float alphaCarbonDistance(const std::vector< std::vector<float> >& loop)
{
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



/*
    Measures the distance in angstroms between two beta carbons of a loop
*/
float betaCarbonDistance(const std::vector< std::vector<float> >& loop)
{
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



/*
    Measures the distance between two points in 3D space (here, atoms)
*/
float atomDistance(const std::vector<float>& atom1, const std::vector<float>& atom2)
{
  float x = pow( atom1[0] - atom2[0] , 2);
  float y = pow( atom1[1] - atom2[1] , 2);
  float z = pow( atom1[2] - atom2[2] , 2);
  return sqrt( x + y + z );
}



/*
    Measures the distance between two points in 3D space (here, atoms)
    Skips square rooting to save some time (use squared distance when checking collisions)
*/
float atomDistanceFast(const std::vector<float>& atom1, const std::vector<float>& atom2)
{
  float x = (atom1[0] - atom2[0]) * (atom1[0] - atom2[0]);
  float y = (atom1[1] - atom2[1]) * (atom1[1] - atom2[1]);
  float z = (atom1[2] - atom2[2]) * (atom1[2] - atom2[2]);
  return (x + y + z);
}



/*
    Root mean square deviation
*/
float RMSD(const std::vector< std::vector<float> >& loop1, const std::vector< std::vector<float> >& loop2)
{
  //Loops (or molecules) need to be of the same size
  if (loop1.size() != loop2.size()){
    std::cerr << "Loops not of the same size!" << std::endl;
    std::cerr << "Loop 1: " << loop1.size() << std::endl;
    std::cerr << "Loop 2: " << loop2.size() << std::endl;
    throw 0;
  }

  //Make working copies, Initialize sum and n
  std::vector< std::vector<float> > loop1_copy = loop1;
  std::vector< std::vector<float> > loop2_copy = loop2;
  unsigned int n = loop1.size();
  float sum = 0.0;

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



/*
    Outputs a PDB file with the original (non-alanine) sequence
*/
void pdbOut (const std::vector< std::vector<float> >& loop, const std::vector<char> residues, char* filename)
{
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
  Outputs a loop with its native sequence replaced with alanines
*/
void pdbOut (const std::vector< std::vector<float> >& loop, char* filename)
{

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



/*
    Detect if a given loop is colliding with the scaffold protein
*/
bool Protein::isCollision (const std::vector<std::vector<float> >& insertion, int start, int end, float threshold) const
{

  // Check all residues before insertion
  // Skip residue before the anchor
  for (int i = 0; i < ( (start-2) * 5); ++i){
    for (unsigned int j = 5; j < insertion.size() -5; ++j){
      if (atomDistanceFast(insertion[j], backbone_coordinates[i]) < threshold){
        return true;
      }
    }
  }

  // Check all residues after insertion
  // Skip the residue right after the anchor
  for (unsigned int i = ( (end + 2) * 5); i < backbone_coordinates.size(); ++i){
    for (unsigned int j = 0; j < insertion.size(); ++j){
      if (atomDistanceFast(insertion[j], backbone_coordinates[i]) < threshold){
        return true;
      }
    }
  }

  return false;
}



/*
    Finds and returns the indices where each chain starts.
*/
const std::vector<int> Protein::findChains()
{
  std::vector<int> chain_start_indices;

  // Just return an empty vector if some idiot (read: me) accidently calls this on an empty protein
  if (backbone_coordinates.size() == 0){
    return chain_start_indices;
  }

  char current_chain = chain[0];
  chain_start_indices.push_back(0);

  for (unsigned int i = 1; i < backbone_coordinates.size(); ++i){
    if (chain[i] != current_chain){
      chain_start_indices.push_back(i);
      current_chain = chain[i];
    }
  }

  return chain_start_indices;
}



/*
    Splits a protein into its constituent chains.
*/
std::vector<Protein> Protein::splitChains()
{
  std::vector<Protein> return_chains;
  std::vector<int> chain_start_indices = findChains();

  // Looks messy but it's the easiest way to slice all the vectors
  for (unsigned int i = 0; i < chain_start_indices.size(); ++i){
    if (i == chain_start_indices.size() - 1){
      std::string new_identifier;
      new_identifier = identifier + "_" + chain[ chain_start_indices[i] ];
      std::vector<std::vector<float> > new_backbone_coordinates( backbone_coordinates.begin() + chain_start_indices[i], backbone_coordinates.begin() + backbone_coordinates.size() );
      std::vector<char> new_residue_type( residue_type.begin() + chain_start_indices[i], residue_type.begin() + residue_type.size() );
      std::vector<int> new_atom_residue_numbering( atom_residue_numbering.begin() + chain_start_indices[i], atom_residue_numbering.begin() + atom_residue_numbering.size() );
      std::vector<char> new_chain( chain.begin() + chain_start_indices[i], chain.begin() + chain.size() );
      Protein tmp_protein(new_identifier, new_backbone_coordinates, new_residue_type, new_atom_residue_numbering, new_chain);
      return_chains.push_back(tmp_protein);
    }
    else{
      std::string new_identifier;
      new_identifier = identifier + "_" + chain[ chain_start_indices[i] ];
      std::vector<std::vector<float> > new_backbone_coordinates( backbone_coordinates.begin() + chain_start_indices[i], backbone_coordinates.begin() + chain_start_indices[i+1] );
      std::vector<char> new_residue_type( residue_type.begin() + chain_start_indices[i], residue_type.begin() + chain_start_indices[i+1] );
      std::vector<int> new_atom_residue_numbering( atom_residue_numbering.begin() + chain_start_indices[i], atom_residue_numbering.begin() + chain_start_indices[i+1]);
      std::vector<char> new_chain( chain.begin() + chain_start_indices[i], chain.begin() + chain_start_indices[i+1] );
      Protein tmp_protein(new_identifier, new_backbone_coordinates, new_residue_type, new_atom_residue_numbering, new_chain);
      return_chains.push_back(tmp_protein);
    }
  }

  return return_chains;
}



/*
    Collects backbone atom lines from a pdb file, returns a vector
*/
std::vector<std::string> Protein::collectAtoms(const char* filename)
{
  std::vector<std::string> cleaned_file;
  //Try to open the PDB file for reading
  std::ifstream pdb_str(filename);
  if (!pdb_str.good()) {
    std::string err(filename);
    throw std::runtime_error("Couldn't open PDB file " + err + "\n");
    return cleaned_file;
  }

  //Check to make sure file isn't empty
  if ( pdb_str.peek() == std::ifstream::traits_type::eof() ){
    std::string err(filename);
    throw std::runtime_error("PDB file was empty!" + err + "\n");
    return cleaned_file;
  }

  // Grab backbone atom lines
  std::string line;
  std::getline(pdb_str, line);
  while (!pdb_str.eof() && line.size() > 3){
    if (line.substr(0, 4) == "ATOM"){
      if ((line.substr(13,2) == "N " || line.substr(13,2) == "CA" || line.substr(13,2) == "C " || line.substr(13,2) == "O " || line.substr(13,2) == "CB")
                                                                                                                            && line.substr(16,1) != "B"){
        cleaned_file.push_back(line);
        std::getline(pdb_str, line);
      }
      else{
        std::getline(pdb_str, line);
      }
    }
    else{
      std::getline(pdb_str, line);
    }
  }

  return cleaned_file;
}



/*
    Adds a beta carbon onto a 4-atom residue using reference alanine
*/
void Protein::addBetaCarbon(std::vector<std::vector<float> >& residue)
{
  std::vector<std::vector<float> > truncated_alanine(alanine.begin(), alanine.end() - 1);
  std::vector<float> cb = alanine[4];

  try{
    std::pair< std::vector< std::vector<float> > , std::vector<float> >
    transformation = superimposer(residue, truncated_alanine, 4);
    superimposer_move(cb, transformation.first, transformation.second);
    residue.push_back(cb);
  }
  catch(int){
    throw std::runtime_error("Superimposer failed! \n");
  }

  return;
}



/*
    Parses vector containing pdb backbone atom lines
*/
void Protein::parseBackbone(const std::vector<std::string> &cleaned_file)
{
  std::string AA;
  int residue_number;
  std::string chainID;
  std::vector<float> xyz(3);

  for (unsigned int i = 0; i < cleaned_file.size(); /*Do nothing*/){
    AA = cleaned_file[i].substr(17,3);
    residue_number = atoi( cleaned_file[i].substr(23,3).c_str() );
    chainID = cleaned_file[i].substr(21,1);

    // Glycines get special treatment :)
    if (AA == "GLY"){
      unsigned int j = i + 4;
      std::vector<std::vector<float> > glycine;
      while (i < j){
        xyz[0] =  atof( cleaned_file[i].substr(30,8).c_str() );
        xyz[1] =  atof( cleaned_file[i].substr(38,8).c_str() );
        xyz[2] =  atof( cleaned_file[i].substr(46,8).c_str() );
        glycine.push_back(xyz);
        ++i;
      }

      addBetaCarbon(glycine);

      // Update protein
      for (unsigned int k = 0; k < 5; ++k){
        residue_type.push_back('A');
        chain.push_back(chainID[0]);
        atom_residue_numbering.push_back(residue_number);
        backbone_coordinates.push_back(glycine[k]);
      }
    }

    // Not a glycine, check to make sure everything's there and add a beta carbon if one is missing
    else{
      unsigned int j = i + 5;
      std::vector<std::vector<float> > residue;

      while (i < j){
        // Is this the same residue? If not then it's probably missing a beta carbon
        std::string current_AA = cleaned_file[i].substr(17,3);
        int current_residue_number = atoi( cleaned_file[i].substr(23,3).c_str() );

        if (current_AA == AA && current_residue_number == residue_number){
          xyz[0] =  atof( cleaned_file[i].substr(30,8).c_str() );
          xyz[1] =  atof( cleaned_file[i].substr(38,8).c_str() );
          xyz[2] =  atof( cleaned_file[i].substr(46,8).c_str() );
          residue.push_back(xyz);
          ++i;
        }
        else{
          addBetaCarbon(residue);
          break;
        }
      }

      // Update protein
      for (unsigned int k = 0; k < 5; ++k){
        residue_type.push_back(amino_codes.find(AA)->second);
        chain.push_back(chainID[0]);
        atom_residue_numbering.push_back(residue_number);
        backbone_coordinates.push_back(residue[k]);
      }
    }
  }

  return;
}



/*
    Appends a parsed protein backbone to a random access file (pdblist.dat)
*/
bool Protein::RAFout(char* filename)
{
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
