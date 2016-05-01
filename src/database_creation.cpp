
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

  /*

  //Construct pdbselect database (Run with with process.sh)
  if (argc == 4 && argv[2] == "pdbselect"){
    try{

    Protein pdb_select(argv[1]);
    assert( pdb_select.getCoordinates().size() % 5 == 0);
    pdb_select.RAF_out(argv[3]);
    return 0;
    }

    catch (int i){

      std::cout << "Couldn't process " << argv[1] << std::endl;
      return 0;

    }

  }

*/

/*

This block constructs the other database files. Run this after protein random access file is constructed

    vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(20, 0) ) );

    //Count loops of each length
    for (int i = 1; i < 21; ++i){
      countLoops(argv[1], i, grid);
    }

    //Write loops counted in the grid
    writeLoops(argv[1], argv[2], argv[3], grid);
    return 0;

*/

  std::cout << "Finished processing " << argv[1] << std::endl << std::endl;
  return 0;
}
