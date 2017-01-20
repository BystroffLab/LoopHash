
#include "protein.h"
#include "loop_generation.h"

/*
  This file is used to construct pdblist.dat
  Appends to the end of the .dat file, run using shell script on folder of pdbs
*/

int main(int argc, char* argv[]){

  // [0]exe [1]filename [2]outfile (master protein list)
  //Construct pdbselect database (Run with with process.sh)
  if (argc == 3){
    try{

    Protein pdb_select(argv[1]);
    assert( pdb_select.getCoordinates().size() % 5 == 0);
    pdb_select.RAFout(argv[2]);
    std::cout << "Wrote out " << argv[1] << std::endl;
    return 0;

    }
    catch (int i){

      std::cout << "Couldn't process " << argv[1] << std::endl;
      return 1;

    }

  }

}
