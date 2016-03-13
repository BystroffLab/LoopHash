#include "loop_generation.h"

void test_raf_in (char* filename){
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

void countLoops(char* filename, int loop_length, vec_3D& grid){
  /*
  Counts loops in random access file. Keeps track of loops in 3D array.
  Truncate d(CA-CA) and d(CB-CB), multiply by 10 to get index
  */

  //Try to open PDBselect file
  std::ifstream iofile(filename, std::fstream::binary);
  if (!iofile){
    std::cerr << "Can't open " << filename << " to write." << std::endl;
    return;
  }

  //Check to make sure loop is between 2-10 residues long
  if (loop_length == 0 || loop_length == 1){
    std::cerr << "Error: Loop length is too short. " << std::endl;
    return;
  }
  else if (loop_length > 10){
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

void writeLoops(char*filename, vec_3D& grid){
  //Loop through every position in grid, look for nonzero items
  //Assign record position to index 2

  //Loop through master protein file for each loop length, write position/length to loop file
  //Increment counter that keeps track of how filled the record is (index 1)
  //Number of final loops record should hold is index 0

  //Build 3d grid of pairs <record position, length> , write to file


  return;
}
