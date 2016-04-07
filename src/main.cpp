
#include "protein.h"
#include "loop_generation.h"

int main(int argc, char* argv[]){





      // ./query.exe  pdbselect.prot  db.loop  grid.arr  min  max  ca  cb  len  output

      std::pair< std::vector<std::vector<std::vector<float> > > , std::vector<std::vector<char> > >
      pdb_output;
      float ca = atof(argv[6]);
      float cb = atof(argv[7]);
      int len = atoi(argv[8]);

      if (argc < 10 ){
        std::cerr << "Not enough arguments!" << std::endl;
        return 0;
      }
      pdb_output = db_query(ca, cb, len , argv[1], argv[2], argv[3]);





      //Try pdb output
      /*
      if (pdb_output.first.size() > 0){
        PDB_out(pdb_output.first[0], pdb_output.second[0], argv[9]);
      }
      */


      //Output all results
      for (unsigned int i = 0; i < pdb_output.first.size(); ++i){
        PDB_out(pdb_output.first[i], pdb_output.second[i], argv[9]);
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
