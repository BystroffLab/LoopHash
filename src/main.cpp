//--- PDB Parser ---//
#include "protein.h"
#include "loop_generation.h"
#include <windows.h>

int main(int argc, char* argv[]){


  std::cout << "Enter query in form 'query dCA dCB length (in residues)'" << std::endl;
  std::cout << "Enter 'stop' to end query." << std::endl << std::endl;

  while (true){
    float ca;
    float cb;
    int len;
    std::string cmd;

    std::cin >> cmd >> ca >> cb >> len;
    if (cmd == "query"){

      //Start timer
      __int64 ctr1 = 0, ctr2 = 0, freq = 0;
      if ( QueryPerformanceCounter( (LARGE_INTEGER *) &ctr1) != 0){

        //Run query, print results
        db_query(ca, cb, len, argv[1], argv[2], argv[3]);

        //End timer
        QueryPerformanceCounter((LARGE_INTEGER *) &ctr2);
        QueryPerformanceCounter((LARGE_INTEGER *) &freq);

        std::cout << "Query took " << ((ctr2 - ctr1) * 1.0 / freq) << " us" << std::endl;
      }
    }

    else if (cmd == "stop"){
      break;
    }

    else{
      std::cout << "Command not recognized." << std::endl;
      continue;
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





  std::cout << "Finished querying." << std::endl << std::endl;
  return 0;
}
