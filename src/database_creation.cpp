#include "protein.h"
#include "loop_generation.h"

// Constructs the other database files given processed list of proteins
// [0]exe [1]infile(pdblist.dat) [2]outfile(looplist.dat) [3]outfile(grid.dat)

int main(int argc, char* argv[]){

      vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(20, 0) ) );

      //Count loops of each length
      for (int i = 1; i < 21; ++i){
        countLoops(argv[1], i, grid);
      }

      //Write loops counted in the grid
      writeLoops(argv[1], argv[2], argv[3], grid);


    std::cout << "Finished processing " << argv[1] << std::endl << std::endl;
    return 0;

}
