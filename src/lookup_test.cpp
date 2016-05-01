//Script to characterize loop search
#include "protein.h"
#include "loop_generation.h"
#include "lookup.h"

int main(int argc, char* argv[]){

  //Create protein
  Protein protein(argv[1]);
  Protein* protein_pointer = &protein;

  //Create lookup object
  Lookup test_lookup( protein_pointer, atoi(argv[2]), atoi(argv[3]) );
  test_lookup.run();

  //Check sorting and cutoff violator removal
  std::cout << "Good loops" << test_lookup.size() << std::endl;

  std::list<Loop> test_lookup_results = test_lookup.getResults();
  std::list<Loop>::iterator itr = test_lookup_results.begin();
  for( ; itr != test_lookup_results.end(); ++itr){
    std::cout << "RMSD: " << itr->rmsd <<  "Length: " << itr->coordinates.size() / 5 << std::endl;
  }

  return 0;
}
