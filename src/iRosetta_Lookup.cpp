//Lookup function to be used in conjunction with InteractiveRosetta
#include "protein.h"
#include "loop_generation.h"
#include "lookup.h"

//iRosetta Loop Search. Set parameters, grab results, eliminate junk results, write out pdbs with filenames
//that iRosetta knows to pick up.
int main(int argc, char* argv[]){

  //Command line arguments should be in the format:
  //  ./iRosetta_Lookup.exe   [1]input.pdb   [2]pdblist.dat [3]looplist.dat [4]grid.arr  [5]anchor start  [6]anchor end
  //  [7]Min length           [8]Max length  [9]Min results [10]Max results

  //Parse input protein, grab loop
  Protein protein(argv[1]);
  Protein* protein_pointer = &protein;

  /*
  --DEBUGGING--
  std::vector<std::vector<float>> coords = protein.getCoordinates();
  std::string dout = "debug_out.pdb";
  char* debug_out = const_cast<char*>( dout.c_str() );
  PDB_out(coords, debug_out);
  */

  // Are all the command line parameters there?
  if (argc < 11){
    std::cout << 0;
    std::cerr << "Error: Not enough command line arguments!";
    return 1;
  }


  // Optional RMSD threshold
  float threshold = 1.0;
  if (argc == 12){
    threshold = atoi(argv[11]);
  }


  // Create lookup, change default params
  Lookup lookup(protein_pointer, atoi(argv[5]), atoi(argv[6]));
  lookup.setDB(argv[2], argv[3], argv[4]);
  lookup.setRange(atoi(argv[7]), atoi(argv[8]));
  lookup.setMin(atoi(argv[9]));
  lookup.setMax(atoi(argv[10]));
  lookup.setCutoff(threshold);
  lookup.run();


  //Output all loops
  std::list<Loop> results = lookup.getResults();
  std::list<Loop>::iterator itr = results.begin();
  for (int i = 1; itr != results.end(); ++itr, ++i){
    //output named so that irosetta can pick up results
    std::stringstream ss;
    ss << i;
    std::string fout = "loopout_" + ss.str() + ".pdb";
    char* filename = const_cast<char*>( fout.c_str() );
    PDB_out(itr->coordinates, filename);
  }


  //This should be the only stdout in the lookup, necessary for iRosetta to figure out how many results were output
  std::cout << results.size();
  return 0;
}
