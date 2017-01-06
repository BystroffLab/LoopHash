//Lookup function to be used in conjunction with InteractiveRosetta
#include "protein.h"
#include "loop_generation.h"
#include "lookup.h"

void irosettaOutput(const Lookup& lookup)
{
  // Output all loops
  std::list<Loop> results = lookup.getResults();
  std::list<Loop>::iterator itr = results.begin();

  for (int i = 1; itr != results.end(); ++itr, ++i){
    // Output named so that irosetta can pick up results
    std::stringstream ss;
    ss << i;
    std::string fout = "loopout_" + ss.str() + ".pdb";
    char* filename = strdup(fout.c_str());
    PDB_out(itr->coordinates, filename);
  }

  // Let irosetta know how many results there were
  std::cout << results.size();

  return;
}



//iRosetta Loop Search. Set parameters, grab results, eliminate junk results, write out pdbs with filenames
//that iRosetta knows to pick up.
int main(int argc, char* argv[])
{

  //Command line arguments should be in the format:
  //  ./iRosetta_Lookup.exe   [1]input text file

  // Create lookup, change default params
  Lookup lookup(argv[1]);
  lookup.run();
  irosettaOutput(lookup);

  return 0;
}
