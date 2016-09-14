
#include "protein.h"
#include "loop_generation.h"

/*
This file is used to construct
-pdbselect database
-the loop database of pointers to the pdb select
 file
-the 3D grid file that points to the loop database. The dimensions of the 3D
grid are 500 x 500 x 20 (x2)

Could be a little more user friendly
*/

int main(int argc, char* argv[]){

  // [0]exe [1]filename [2]outfile (master protein list)
  //Construct pdbselect database (Run with with process.sh)
  if (argc == 3){
    try{

    Protein pdb_select(argv[1]);
    assert( pdb_select.getCoordinates().size() % 5 == 0);
    pdb_select.RAF_out(argv[2]);
    std::cout << "Wrote out " << argv[1] << std::endl;
    return 0;

    }
    catch (int i){

      std::cout << "Couldn't process " << argv[1] << std::endl;
      return 0;

    }

  }

}
