#include "protein.h"
#include "loop_generation.h"

// Constructs the other database files given processed list of proteins
// [0]exe [1]infile(.pdblist) [2]outfile(.looplist) [3]outfile(.3dgrid)

int main(int argc, char* argv[]){

      vec_3D grid(500, std::vector<std::vector<int> >(500, std::vector<int>(20, 0) ) );

      //Count loops of each length
      for (int i = 1; i < 21; ++i){
        countLoops(argv[1], i, grid);
      }

      //Write loops counted in the grid
      writeLoops(argv[1], argv[2], argv[3], grid);
      return 0;


    std::cout << "Finished processing " << argv[1] << std::endl << std::endl;
    return 0;

}
