//--- PDB Parser ---//
#include "protein.h"
#include "loop_generation.h"

int main(int argc, char* argv[]){
  //typedef std::vector<std::vector<std::vector<float> > > looplist;
  //typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;

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


  /*
  for (int i = 50; i < 250; ++i){
    db_query(i, i, 4, argv[1], "db.loop", "grid.arr");
    db_query(i, i, 5, argv[1], "db.loop", "grid.arr");
    db_query(i, i, 6, argv[1], "db.loop", "grid.arr");
    db_query(i, i, 7, argv[1], "db.loop", "grid.arr");

  }
  */

  for (int i = 0; i < 500; ++i){
    for (int j = 0; j < 500; ++j){
      for (int k = 4; k < 8; ++k){
        db_query(i, j, k, argv[1], "db.loop", "grid.arr");
      }
    }
  }



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


  std::cout << "Finished processing." << std::endl << std::endl;
  return 0;
}
