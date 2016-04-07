#include "loop_generation.h"


/*         --PDBselect random access file--
  Every residue record is 70 bytes:
  4 letter pdb code (4 bytes) - 1 letter aa code (1 byte) - int residue # (4 bytes) - char chain id (1 byte) - 10 bytes total
  Coordinates for N, CA, C, O, CB 12 bytes per atom (4 bytes per float) - 60 total
  This function just reads and prints pdb select file for debugging purposes
*/
void pdbselect_raf_in (char* filename){
  //Try to open
  std::ifstream iofile(filename, std::fstream::binary);
  if (!iofile){
    std::cerr << "Can't open " << filename << " to write." << std::endl;
    return;
  }

  //Get length of file
  unsigned int size;
  iofile.seekg(0, iofile.end);
  size = (unsigned int) iofile.tellg();
  std::cout << "size: " << size << std::endl;
  iofile.seekg(0, iofile.beg); //Return to beginning of file

  //Read and output
  char pdb_code[5];
  char aa;
  int residue_number;
  char chain;

  while (iofile.tellg() < size){
    iofile.read( pdb_code, (sizeof(char) * 4) );
    iofile.read( (char*)& aa, sizeof(char) );
    iofile.read( (char*)& residue_number, sizeof(int) );
    iofile.read( (char*)& chain, sizeof(char) );
    pdb_code[4] = '\0';

    std::cout << "PDB " << pdb_code << "  Residue " << aa << residue_number << " Chain " << chain << std::endl;

    for (int i = 0; i < 5; ++i){
      float x; float y; float z;
      iofile.read( (char*)& x, sizeof(float) );
      iofile.read( (char*)& y, sizeof(float) );
      iofile.read( (char*)& z, sizeof(float) );
      std::cout << x << " " << y << " " << z << std::endl;
    }
  }

  iofile.close();
  return;
}




/*
Counts loops in random access file. Keeps track of loops in 3D array.
Truncate d(CA-CA) and d(CB-CB), multiply by 10 to get index
*/
void countLoops(char* filename, int loop_length, vec_3D& grid){

  //Try to open PDBselect file
  std::ifstream iofile(filename, std::fstream::binary);
  if (!iofile){
    std::cerr << "Can't open " << filename << " to read." << std::endl;
    return;
  }

  //Check to make sure loop is between 2-20 residues long
  if (loop_length == 0 || loop_length == 1){
    std::cerr << "Error: Loop length is too short. " << std::endl;
    return;
  }
  else if (loop_length > 20){
    std::cerr << "Error: Loop length is too long. " << std::endl;
    return;
  }
  else{
    std::cout << "Counting loops of size " << loop_length << std::endl;
  }

  //Get length of file
  unsigned int size;
  iofile.seekg(0, iofile.end);
  size = (unsigned int) iofile.tellg();
  iofile.seekg(0, iofile.beg); //Return to beginning of file


  //Start reading file, initialize some variables first
  unsigned int absolute_position = 0;       //Keep track of where the loop being read starts
  char current_pdb_code[5] = {0};           //Make sure each loop only contains residues from 1 protein
  char temporary_pdb_code[5] = {0};         //^^^
  char aa;                                  //Token to hold residue name
  int residue_number;                       //Token to hold residue number
  char chain;                               //Token to hold chain (eventually will be criteria)
  int counter = 0;                          //How far are we through the current loop?
  std::vector< std::vector<float> > loop;   //Holds backbone coordinates for some calculations

  //While we're not at the end of the file
  while (iofile.tellg() < size){

    counter = 0;
    //Start generating the loop
    while ( counter < loop_length  &&  iofile.tellg() < size ){
      //std::cout << "COUNTER: " << counter << std::endl << std::endl;
      //Check to make sure the residue being read is from the same protein (don't care if this is the start of new loop)
      iofile.read( temporary_pdb_code, (sizeof(char) * 4) );
      temporary_pdb_code[4] = '\0';
      if ( counter > 0 && strcmp(temporary_pdb_code, current_pdb_code) != 0 ){
        break;
      }
      else {
        //arrays aren't assignable, use std copy to accomplish
        //current_pdb_code = temporary_pdb_code;
        std::copy(temporary_pdb_code, temporary_pdb_code + 5, current_pdb_code);
      }

      //Don't really need this data (skip past this with seekg(tellg + 6)?)
      iofile.read( (char*)& aa, sizeof(char) );
      iofile.read( (char*)& residue_number, sizeof(int) );
      iofile.read( (char*)& chain, sizeof(char) );


      //Get coordinates for this residue, push back into loop
      std::vector<float> xyz;
      for (int i = 0; i < 5; ++i){
        float x; float y; float z;
        iofile.read( (char*)& x, sizeof(float) );
        iofile.read( (char*)& y, sizeof(float) );
        iofile.read( (char*)& z, sizeof(float) );

        xyz.push_back(x);
        xyz.push_back(y);
        xyz.push_back(z);
        loop.push_back(xyz);
        xyz.clear();
        //std::cout << x << " " << y << " " << z << std::endl;
      }


    ++counter;
    } //End of current loop read

    //Make sure the whole length requested was read (5 atoms per residue)
    if ( int(loop.size() / 5) == loop_length ){
      //Get d(CA_CA), d(CB_CB), truncate, convert to int
      float CA_CA = ca_ca_dist(loop);
      float CB_CB = cb_cb_dist(loop);

      //Use values to find bucket, store count of each loop type for now

      //Truncate floats
      char sz[64];
      double lf = CA_CA;
      sprintf(sz, "%.1lf\n", lf);
      double lf2 = atof(sz);
      CA_CA = lf2;

      lf = CB_CB;
      sprintf(sz, "%.1lf\n", lf);
      lf2 = atof(sz);
      CB_CB = lf2;

      //Multiply by 10 to get integer indices
      int ca_index = CA_CA * 10;
      int cb_index = CB_CB * 10;


      //Update grid
      if (ca_index < 500 && cb_index < 500){
        ++grid[ca_index][cb_index][loop_length - 1];
      }

    }//Finished doing stuff with loop

    //Get vector ready for next set of coordinates
    loop.clear();

    //Move forward 1 residue from start of loop, update absolute position
    iofile.seekg(absolute_position + 70);
    absolute_position = absolute_position + 70;

  }//End of file

  return;
}




/*
This function takes a 3D array of loop counts, and uses it to construct a 4D vector used
for keeping track of loop writes to a random access file with records of variable length
that contain rows of int(pdbselect fileposition), int(loop length)

After writing to loop file, a random access file version of the 4D vector is written.
*/

void writeLoops(char* pdbselect, char* loopfile, char* gridfile, const vec_3D& grid){
  //Initialize 4D vector to keep track of writes to loopfile
  std::vector<std::vector<std::vector<std::vector<int> > > >
  write_progress(500, std::vector<std::vector<std::vector<int> > >(500, std::vector<std::vector<int > >(20, std::vector<int>(3,0) ) ) );
  // [i][j][k][0] - Total number of loops bin should hold
  // [i][j][k][1] - Number of loops currently in bin
  // [i][j][k][2] - Start of record (in bytes)


  //Keep track of the start of where the loops are going to be stored (in bytes)
  unsigned int file_pos = 0;

  for (unsigned int i = 0; i < grid.size(); ++i){
    for (unsigned int j = 0; j < grid[i].size(); ++j){
      for (unsigned int k = 0; k < grid[i][j].size(); ++k){
        //Grid was generated by first pass to determine how many loops in each bucket.
        //That number is the total record size for this bin
        write_progress[i][j][k][0] = grid[i][j][k];

        //Figure out where these loops are going to be in the loop db file
        if (grid[i][j][k] > 0){
          write_progress[i][j][k][2] = file_pos;
          file_pos = file_pos + (8 * grid[i][j][k]);
        }

      }
    }
  }

  //Build empty output file for loops of size determined by prior loop
  std::ofstream loop_out(loopfile, std::ios::binary | std::ios::out);
  loop_out.seekp(file_pos - 1);
  loop_out.write("", 1);


  //Loop through master protein file for each loop length, write position/length to loop file
  //Increment counter that keeps track of how filled the record is (index 1)
  //Number of final loops record should hold is index 0

  //Try to open PDBselect file
  std::ifstream infile(pdbselect, std::fstream::binary);
  if (!infile){
    std::cerr << "Can't open " << pdbselect << " to read." << std::endl;
    return;
  }


  //Get length of PDBselect
  unsigned int size;
  infile.seekg(0, infile.end);
  size = (unsigned int) infile.tellg();
  infile.seekg(0, infile.beg); //Return to beginning of file


  //Start reading file, initialize some variables first
  unsigned int absolute_position = 0;       //Keep track of where the loop being read starts
  char current_pdb_code[5] = {0};           //Make sure each loop only contains residues from 1 protein
  char temporary_pdb_code[5] = {0};         //^^^
  char aa;                                  //Token to hold residue name
  int residue_number;                       //Token to hold residue number
  char chain;                               //Token to hold chain (eventually will be criteria)
  int counter = 0;                          //How far are we through the current loop?
  std::vector< std::vector<float> > loop;   //Holds backbone coordinates for some calculations


  for (int loop_length = 2; loop_length < 21; ++loop_length){
    std::cout << loop_length << std::endl;
    absolute_position = 0;
    infile.seekg(0, infile.beg); //Return to beginning of file

    //While we're not at the end of the file
    while (infile.tellg() < size){

      counter = 0;
      //Start generating the loop
      while ( counter < loop_length  &&  infile.tellg() < size ){
        //std::cout << "COUNTER: " << counter << std::endl << std::endl;
        //Check to make sure the residue being read is from the same protein (don't care if this is the start of new loop)
        infile.read( temporary_pdb_code, (sizeof(char) * 4) );
        temporary_pdb_code[4] = '\0';
        if ( counter > 0 && strcmp(temporary_pdb_code, current_pdb_code) != 0 ){
          break;
        }
        else {
          //arrays aren't assignable, use std copy to accomplish
          //current_pdb_code = temporary_pdb_code;
          std::copy(temporary_pdb_code, temporary_pdb_code + 5, current_pdb_code);
        }

        //Don't really need this data (skip past this with seekg(tellg + 6)?)
        infile.read( (char*)& aa, sizeof(char) );
        infile.read( (char*)& residue_number, sizeof(int) );
        infile.read( (char*)& chain, sizeof(char) );


        //Get coordinates for this residue, push back into loop
        std::vector<float> xyz;
        for (int i = 0; i < 5; ++i){
          float x; float y; float z;
          infile.read( (char*)& x, sizeof(float) );
          infile.read( (char*)& y, sizeof(float) );
          infile.read( (char*)& z, sizeof(float) );

          xyz.push_back(x);
          xyz.push_back(y);
          xyz.push_back(z);
          loop.push_back(xyz);
          xyz.clear();
          //std::cout << x << " " << y << " " << z << std::endl;
        }


      ++counter;
      } //End of current loop read

      //Make sure the whole length requested was read (5 atoms per residue)
      if ( int(loop.size() / 5) == loop_length ){
        //Get d(CA_CA), d(CB_CB), truncate, convert to int
        float CA_CA = ca_ca_dist(loop);
        float CB_CB = cb_cb_dist(loop);

        //Use values to find bucket

        //Truncate floats
        char sz[64];
        double lf = CA_CA;
        sprintf(sz, "%.1lf\n", lf);
        double lf2 = atof(sz);
        CA_CA = lf2;

        lf = CB_CB;
        sprintf(sz, "%.1lf\n", lf);
        lf2 = atof(sz);
        CB_CB = lf2;

        //Multiply by 10 to get integer indices
        int ca_index = CA_CA * 10;
        int cb_index = CB_CB * 10;


        //Find correct spot for loop, write pdbselect position (absolute_position, bytes) and length (in residues)
        if (ca_index < 500 && cb_index < 500){
          //std::cout << "Writing " << absolute_position << " " << loop_length << " at " << write_progress[ca_index][cb_index][loop_length-1][2] + write_progress[ca_index][cb_index][loop_length-1][1] * 8 << std::endl;


          loop_out.seekp( write_progress[ca_index][cb_index][loop_length-1][2] + (write_progress[ca_index][cb_index][loop_length-1][1] * 8) );
          loop_out.write( (const char*)& absolute_position, sizeof(unsigned int) );
          loop_out.write( (const char*)& loop_length, sizeof(int) );
          ++write_progress[ca_index][cb_index][loop_length-1][1];


        }


      }//Finished doing stuff with current loop

      //Get vector ready for next set of coordinates
      loop.clear();

      //Move forward 1 residue from start of loop, update absolute position
      infile.seekg(absolute_position + 70);
      absolute_position = absolute_position + 70;

    }//End of file

 }

  //Build 3D grid of pairs <record position, length> , write to file
  std::ofstream grid_out(gridfile, std::ios::binary | std::ios::out);
  grid_out.seekp( (500 * 500 * 20 * 8) - 1);
  grid_out.write("", 1);  //Create empty file
  grid_out.seekp(0);      //Return to start of file

  for (int i = 0; i < 500; ++i){
    for (int j = 0; j < 500; ++j){
      for (int k = 0; k < 20; ++k){
        if (write_progress[i][j][k][0] > 0){
          std::cout << "Loops at" << i << " "<< j << " " << k << " " << write_progress[i][j][k][0] << "Pos: " <<  write_progress[i][j][k][2] << std::endl;
        }
        grid_out.write( (const char*)& write_progress[i][j][k][0], sizeof(int) );
        grid_out.write( (const char*)& write_progress[i][j][k][2], sizeof(int) );
      }
    }
  }

  return;
}




//Function for reading loop pointer database (for debugging)
void readLoops(char* loopfile){

  //Try to open file
  std::ifstream infile(loopfile, std::fstream::binary);
  if (!infile){
    std::cerr << "Can't open " << loopfile << " to read." << std::endl;
    return;
  }

  //Get length of file
  unsigned int size;
  infile.seekg(0, infile.end);
  size = (unsigned int) infile.tellg();
  infile.seekg(0, infile.beg); //Return to beginning of file

  unsigned int pdb_position;
  int length;
  int counter = 0;
  while (infile.tellg() < size){
    infile.read( (char*)& pdb_position, sizeof(unsigned int) );
    infile.read( (char*)& length, sizeof(int) );

    std::cout << "Loop of length " << length << " at pdbselect: " << pdb_position << std::endl;
    ++counter;
  }

  std::cout << "Total number of loops: " << counter << std::endl;

}

//Function for reading grid (for debugging)
void readGrid(char* gridfile){
  //Try to open file
  std::ifstream infile(gridfile, std::fstream::binary);
  if (!infile){
    std::cerr << "Can't open " << gridfile << " to read." << std::endl;
    return;
  }

  //Get length of file
  unsigned int size;
  infile.seekg(0, infile.end);
  size = (unsigned int) infile.tellg();
  infile.seekg(0, infile.beg); //Return to beginning of file

  int i = 0;
  int loops;
  int record;

  while (infile.tellg() < size){
    infile.read( (char*)& loops, sizeof(unsigned int) );
    infile.read( (char*)& record, sizeof(int) );

    if (loops != 0){
      std::cout << i << ":  Loops: " << loops << " Record start: " << record << std::endl;
    }

    i = i + 8;
  }

}





//Simple database query function
/*  GRID ACCESS
4D array is of the dimensions 500x500x20x2 with each entry being 4 bytes
To access grid[i][j][k][l] from file use seekg( (80000 * i) + (160 * j) + (8 * k) + (4 * l) )
*/
std::pair< std::vector<std::vector<std::vector<float> > > , std::vector<std::vector<char> > >
db_query(float CA_CA, float CB_CB, int loop_length, char* pdbselect, char* loopfile, char* gridfile){
  std::vector<std::vector<std::vector<float> > > return_val;
  std::vector<std::vector<char> > residues;

  //Truncate floats, multiply by 10 to get int
  char sz[64];
  double lf = CA_CA;
  sprintf(sz, "%.1lf\n", lf);
  double lf2 = atof(sz);
  CA_CA = lf2;

  lf = CB_CB;
  sprintf(sz, "%.1lf\n", lf);
  lf2 = atof(sz);
  CB_CB = lf2;

  int dCA = CA_CA * 10;
  int dCB = CB_CB * 10;


  //Open grid file to find loop pointers in loop file
  std::ifstream grid_in(gridfile, std::fstream::binary);
  if (!grid_in){
    std::cerr << "Can't open " << gridfile << " to read." << std::endl;
    throw 0;
  }

  int number_of_loops;
  int record_start;
  grid_in.seekg( (80000 * dCA) + (160 * dCB) + (8 * (loop_length - 1) ) );
  grid_in.read( (char*)& number_of_loops, sizeof(int) );
  grid_in.read( (char*)& record_start, sizeof(int) );

  grid_in.close();

  //Were any loops found?
  if (number_of_loops == 0){
    std::cout << "No loops found." << std::endl;
    return std::make_pair(return_val, residues);
  }

  //Open loop file to find pdbselect pointers
  std::ifstream loop_in(loopfile, std::fstream::binary);
  if (!loop_in){
    std::cerr << "Can't open " << loopfile << " to read." << std::endl;
    throw 0;
  }

  std::vector<std::vector<int> > loop_pointers;
  loop_in.seekg( record_start );

  //Read in all locations of loops and their lengths
  //This file can actually be half as big because we know the lengths already
  unsigned int pdb_position;
  std::vector<int> tmp;
  for (int i = 0; i < number_of_loops; ++i){
    loop_in.read( (char*)& pdb_position, sizeof(unsigned int) );
    loop_in.read( (char*)& loop_length, sizeof(unsigned int) );

    tmp.push_back(pdb_position);
    tmp.push_back(loop_length);
    loop_pointers.push_back(tmp);
    tmp.clear();
  }

  loop_in.close();

  //Open pdb file to get actual loops
  std::ifstream pdb_in(pdbselect, std::fstream::binary);
  if ( !pdb_in.good() ){
    std::cerr << "Can't open " << pdbselect << " to read." << std::endl;
    throw 0;
  }

  char pdb_code[5];
  char aa;
  int residue_number;
  char chain;
  std::vector<std::vector<float> > loop;
  std::vector<char> loop_residues;
  std::vector<float> xyz;


  //Output
  for (unsigned int i = 0; i < loop_pointers.size(); ++i){
    //std::cout << "Loop at " << loop_pointers[i][0] << " of length " << loop_pointers[i][1] << std::endl;
    pdb_in.seekg(loop_pointers[i][0]);
    //std::cout << "----------------- Loop " << i + 1 << " -----------------" << std::endl;
    for (int j = 0; j < loop_pointers[i][1]; ++ j){

      pdb_in.read( pdb_code, (sizeof(char) * 4) );
      pdb_in.read( (char*)& aa, sizeof(char) );
      pdb_in.read( (char*)& residue_number, sizeof(int) );
      pdb_in.read( (char*)& chain, sizeof(char) );
      pdb_code[4] = '\0';

      //std::cout << "PDB " << pdb_code << "  Residue " << aa << residue_number << " Chain " << chain << std::endl;
      loop_residues.push_back(aa);

      for (int i = 0; i < 5; ++i){
        float x; float y; float z;
        pdb_in.read( (char*)& x, sizeof(float) );
        pdb_in.read( (char*)& y, sizeof(float) );
        pdb_in.read( (char*)& z, sizeof(float) );
        xyz.push_back(x);
        xyz.push_back(y);
        xyz.push_back(z);
        loop.push_back(xyz);
        xyz.clear();

        //std::cout << x << " " << y << " " << z << std::endl;
      }


    } //End loop

    //std::cout << "Real dCA: " << ca_ca_dist(loop) << " dCA Query: " << CA_CA << std::endl;
    //std::cout << "Real dCB: " << cb_cb_dist(loop) << " dCA Query: " << CB_CB << std::endl;

    return_val.push_back(loop);
    residues.push_back(loop_residues);
    loop_residues.clear();
    loop.clear();

    //std::cout << std::endl << std::endl;

  }

  return std::make_pair(return_val, residues);

}

void continuous_query_helper(std::vector<std::vector<std::vector<float> > >& results, float CA_CA, float CB_CB, int loop_length,
  char* pdbselect, char* loopfile, char* gridfile){

    //std::vector<std::vector<std::vector<float> > > return_val;
    std::vector<std::vector<char> > residues;

    //Truncate floats, multiply by 10 to get int
    char sz[64];
    double lf = CA_CA;
    sprintf(sz, "%.1lf\n", lf);
    double lf2 = atof(sz);
    CA_CA = lf2;

    lf = CB_CB;
    sprintf(sz, "%.1lf\n", lf);
    lf2 = atof(sz);
    CB_CB = lf2;

    int dCA = CA_CA * 10;
    int dCB = CB_CB * 10;


    //Open grid file to find loop pointers in loop file
    std::ifstream grid_in(gridfile, std::fstream::binary);
    if (!grid_in){
      std::cerr << "Can't open " << gridfile << " to read." << std::endl;
      throw 0;
    }

    int number_of_loops;
    int record_start;
    grid_in.seekg( (80000 * dCA) + (160 * dCB) + (8 * (loop_length - 1) ) );
    grid_in.read( (char*)& number_of_loops, sizeof(int) );
    grid_in.read( (char*)& record_start, sizeof(int) );

    grid_in.close();

    //Were any loops found?
    if (number_of_loops == 0){
      std::cout << "No loops found." << std::endl;
      return;
    }

    //Open loop file to find pdbselect pointers
    std::ifstream loop_in(loopfile, std::fstream::binary);
    if (!loop_in){
      std::cerr << "Can't open " << loopfile << " to read." << std::endl;
      throw 0;
    }

    std::vector<std::vector<int> > loop_pointers;
    loop_in.seekg( record_start );

    //Read in all locations of loops and their lengths
    //This file can actually be half as big because we know the lengths already
    unsigned int pdb_position;
    std::vector<int> tmp;
    for (int i = 0; i < number_of_loops; ++i){
      loop_in.read( (char*)& pdb_position, sizeof(unsigned int) );
      loop_in.read( (char*)& loop_length, sizeof(unsigned int) );

      tmp.push_back(pdb_position);
      tmp.push_back(loop_length);
      loop_pointers.push_back(tmp);
      tmp.clear();
    }

    loop_in.close();

    //Open pdb file to get actual loops
    std::ifstream pdb_in(pdbselect, std::fstream::binary);
    if ( !pdb_in.good() ){
      std::cerr << "Can't open " << pdbselect << " to read." << std::endl;
      throw 0;
    }

    char pdb_code[5];
    char aa;
    int residue_number;
    char chain;
    std::vector<std::vector<float> > loop;
    std::vector<char> loop_residues;
    std::vector<float> xyz;


    //Output
    for (unsigned int i = 0; i < loop_pointers.size(); ++i){
      //std::cout << "Loop at " << loop_pointers[i][0] << " of length " << loop_pointers[i][1] << std::endl;
      pdb_in.seekg(loop_pointers[i][0]);
      //std::cout << "----------------- Loop " << i + 1 << " -----------------" << std::endl;
      for (int j = 0; j < loop_pointers[i][1]; ++ j){

        pdb_in.read( pdb_code, (sizeof(char) * 4) );
        pdb_in.read( (char*)& aa, sizeof(char) );
        pdb_in.read( (char*)& residue_number, sizeof(int) );
        pdb_in.read( (char*)& chain, sizeof(char) );
        pdb_code[4] = '\0';

        //std::cout << "PDB " << pdb_code << "  Residue " << aa << residue_number << " Chain " << chain << std::endl;
        loop_residues.push_back(aa);

        for (int i = 0; i < 5; ++i){
          float x; float y; float z;
          pdb_in.read( (char*)& x, sizeof(float) );
          pdb_in.read( (char*)& y, sizeof(float) );
          pdb_in.read( (char*)& z, sizeof(float) );
          xyz.push_back(x);
          xyz.push_back(y);
          xyz.push_back(z);
          loop.push_back(xyz);
          xyz.clear();

          //std::cout << x << " " << y << " " << z << std::endl;
        }


      } //End loop

      //std::cout << "Real dCA: " << ca_ca_dist(loop) << " dCA Query: " << CA_CA << std::endl;
      //std::cout << "Real dCB: " << cb_cb_dist(loop) << " dCA Query: " << CB_CB << std::endl;

      results.push_back(loop);
      residues.push_back(loop_residues);
      loop_residues.clear();
      loop.clear();

      //std::cout << std::endl << std::endl;

    }

    return;

  }

  //Get all loops 2-19 resiudes long
  std::vector<std::vector<std::vector<float> > > continuous_query_wrapper(float CA_CA, float CB_CB,
    char* pdbselect, char* loopfile, char* gridfile){

      std::vector<std::vector<std::vector<float> > > results;

      for (int i = 2; i < 20; ++i){
        continuous_query_helper(results, CA_CA, CB_CB, i, pdbselect, loopfile, gridfile);
      }

      return results;


    }
