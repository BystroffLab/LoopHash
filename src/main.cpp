//--- PDB Parser ---//
#include "protein.h"
#include "loop_generation.h"

int main(int argc, char* argv[]){
  //typedef std::vector<std::vector<std::vector<float> > > looplist;
  //typedef std::pair< std::vector< std::vector<float> > , std::vector<float> > return_val;

  vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(10, 0) ) );
  countLoops(argv[1], 4, grid);
  countLoops(argv[1], 5, grid);
  countLoops(argv[1], 6, grid);
  countLoops(argv[1], 7, grid);

  //Check out the grid
  int loop_count = 0;
  for (unsigned int i = 0; i < grid.size(); ++i){
    for (unsigned int j = 0; j < grid.size(); ++j){
      for (unsigned int k = 0; k < grid[i][j].size(); ++k){
        if (grid[i][j][4] > 0){
          std::cout << grid[i][j][4] << " ";
          ++loop_count;
        }
      }
    }
  }
  std::cout << "Loop count: " << loop_count << std::endl;


  /*
  This try/catch block was used to construct the PDBselect random access file.

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
