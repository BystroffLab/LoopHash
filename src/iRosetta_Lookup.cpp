/*
Lookup function to be used in conjunction with InteractiveRosetta
Will Hooper 2016
Bystroff Lab
*/

#include "protein.h"
#include "loop_generation.h"
#include "lookup.h"

/*
   iRosetta Loop Search. Set parameters, grab results, eliminate junk results,
   write out pdbs with filenames that iRosetta knows to pick up.
*/
int main(int argc, char* argv[])
{

  // Create lookup, change default params
  Lookup lookup(argv[1]);
  lookup.run();
  lookup.iRosettaOutput();
  lookup.writeLog();

  return 0;
}
