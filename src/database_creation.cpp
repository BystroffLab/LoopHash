
#include "protein.h"
#include "loop_generation.h"

int main(int argc, char* argv[]){
  /*
  This file is used to construct
  -pdbselect database
  -the loop database of pointers to the pdb select
   file
  -the 3D grid file that points to the loop database. The dimensions of the 3D
  grid are 500 x 500 x 20 (x2)
  */

  /*
  //Construct pdbselect database (used in conjunction with process.sh)
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

  //Construct loop pointers and 3D grid file
*/
    vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(20, 0) ) );

    //Count loops of each length
    for (int i = 1; i < 21; ++i){
      countLoops(argv[1], i, grid);
    }

    //Write loops counted in the grid
    writeLoops(argv[1], argv[2], argv[3], grid);
    return 0;





  /*
  This try/catch block was used to construct the PDBselect random access file.
  Run process.h with the executable compiled with this block uncommented in the
  folder containing all pdbs to be indexed.

  try{

  Protein pdb_select(argv[1]);
  assert( pdb_select.getCoordinates().size() % 5 == 0);
  pdb_select.RAF_out("test.goo");
  //test_raf_in("test.goo");

  }

  catch (int i){

    std::cout << "Couldn't process " << argv[1] << std::endl;
    return 0;

  }

  */





  /*
  //Generate grid used to count loops/write loops

  vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(10, 0) ) );

  //Count loops of each length
  countLoops(argv[1], 4, grid);
  countLoops(argv[1], 5, grid);
  countLoops(argv[1], 6, grid);
  countLoops(argv[1], 7, grid);

  //Write loops counted in the grid
  writeLoops(argv[1], "db.loop", "grid.arr", grid);
  */





  std::cout << "Finished processing." << std::endl << std::endl;
  return 0;
}
