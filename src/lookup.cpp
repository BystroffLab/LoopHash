#include "lookup.h"

Lookup::Lookup(char* input_file)
{
  //Set default parameters
  length_range[0] = 3;   length_range[1] = 19;
  min_results = 10;      max_results = 25;
  rmsd_cutoff = 1.5;     sequence_filter = "";
  filter = false;        sequence_identity_cutoff = 0.0;
  duplicate_threshold = 1.0;
  symmetry = 1;

  //Default files
  char* db1 = strdup("pdblist.dat");  database_files.push_back(db1);
  char* db2 = strdup("looplist.dat"); database_files.push_back(db2);
  char* db3 = strdup("grid.dat");     database_files.push_back(db3);
  logfile = strdup("indel_log.txt");

  //Parse input, try to find the specified anchors
  parse(input_file);

  try{
    original_loop = scaffold.getLoop(scaffold_start, scaffold_end);
  }
  catch(int){
    logmsg("Couldn't find the requested anchor points. Check inputs and try again.\n");
    writeLog();
    exit(EXIT_FAILURE);
  }

}



/*
    Set a sequence to filter for. Not implemented (yet?)
*/
void Lookup::setSequence(std::string s, float identity)
{
  assert(identity <= 1);
  sequence_filter = s;
  sequence_identity_cutoff = identity;
}



/*
    Set range of loops to look for. If inputs are wrong just use default settings
*/
void Lookup::setRange(int min_length, int max_length)
{
  if (min_length > max_length) {
    logmsg("Specified minimum loop length larger than maximum. Using default lengths.. \n");
    return;
  }
  if (min_length < 3 || max_length > 19) {
    logmsg("Specified loop lengths outside database range. Using default lenghts..\n");
    return;
  }
  length_range[0] = min_length;
  length_range[1] = max_length;
  logmsg("Set search for loops of length " + std::to_string(min_length) + " to " + std::to_string(max_length) + ".\n");
  return;
}



/*
    Comparison used by the sort function
*/
static bool rmsdSort(const Loop& a, const Loop& b)
{
  return (a.rmsd < b.rmsd);
}



/*
    Perform a search with previously set parameters
    TODO: Good lord please separate this into more readable chunks
*/
void Lookup::run()
{
  float CA_CA = ca_ca_dist(original_loop);
  float CB_CB = cb_cb_dist(original_loop);

  // If CA_CA or CB_CB are too big, the database is too small so just quit
  if (CA_CA > 49.9 || CB_CB > 49.9){
    logmsg("ERROR: Anchors are too far apart! Pick closer residues. \n");
    return;
  }


  // Run lookup until we have the minimum number of results asked for
  // Vary the CA-CA and CB-CB distances by .1 angstroms until we have enough results or we hit the end of the database
  float change = 0.0;
  bool firstRun = true;
  while(int(results.size()) < min_results){
    //Pull Loops
    if (firstRun){
      for (int loop_length = length_range[0]; loop_length < length_range[1] + 1; ++loop_length){
        runHelper(CA_CA, CB_CB, loop_length);
      }
    }
    else{
      for (int loop_length = length_range[0]; loop_length < length_range[1] + 1; ++loop_length){
        // Diagonals
        runHelper(CA_CA + change, CB_CB + change, loop_length);
        runHelper(CA_CA - change, CB_CB - change, loop_length);
        runHelper(CA_CA + change, CB_CB - change, loop_length);
        runHelper(CA_CA - change, CB_CB + change, loop_length);

        // Across axes
        runHelper(CA_CA + change, CB_CB,          loop_length);
        runHelper(CA_CA,          CB_CB + change, loop_length);
        runHelper(CA_CA - change, CB_CB         , loop_length);
        runHelper(CA_CA         , CB_CB - change, loop_length);
      }

    }
    firstRun = false;



    // Superpose, check rmsd. If it's below the cutoff, throw out loop. Otherwise, record rmsd in struct
    bool decrement = false;
    std::list<Loop>::iterator loops_itr;
    for (loops_itr = results.begin(); loops_itr != results.end(); ++loops_itr){
      if (decrement == true && loops_itr != results.begin()) {
         --loops_itr;
         decrement = false;
      }
      //Collect anchor atoms for superposition
      std::vector<std::vector<float> > list_superposer_in;
      std::vector<std::vector<float> > original_superposer_in;

      //Get first residue for loop list
      for (unsigned int j = 0; j < 5; ++j){
        list_superposer_in.push_back(loops_itr->coordinates[j]);
      }
      //Get last residue for loop list (can't do both at once because loops not always same size)
      for (unsigned int j = loops_itr->coordinates.size() - 5; j < loops_itr->coordinates.size() ; ++j){
        list_superposer_in.push_back(loops_itr->coordinates[j]);
      }

      //Get first residue for original loop
      for (unsigned int j = 0; j < 5; ++j){
        original_superposer_in.push_back(original_loop[j]);
      }
      //Get last residue for original loop
      for (unsigned int j = original_loop.size() - 5; j < original_loop.size() ; ++j){
        original_superposer_in.push_back(original_loop[j]);
      }

      //Sanity check
      assert(list_superposer_in.size() == 10);
      assert(original_superposer_in.size() == 10);

      //Superpose
      std::pair< std::vector< std::vector<float> > , std::vector<float> > transformation;
      transformation = superimposer(original_superposer_in, list_superposer_in, 10);
      for (unsigned int j = 0; j < list_superposer_in.size(); ++j){
        superimposer_move(list_superposer_in[j], transformation.first, transformation.second);
      }
      for (unsigned int j = 0; j < loops_itr->coordinates.size(); ++j){
        superimposer_move(loops_itr->coordinates[j], transformation.first, transformation.second);
      }

      //Record RMSD, throw out if below cutoff
      loops_itr->rmsd = RMSD(original_superposer_in, list_superposer_in);
      if ( loops_itr->rmsd > rmsd_cutoff && loops_itr != results.begin() ){
        loops_itr = results.erase(loops_itr);
        --loops_itr;
        continue;
      }
      else if ( loops_itr->rmsd > rmsd_cutoff && loops_itr == results.begin() ){
        //Avoid decrementing off the beginning of the list
        loops_itr = results.erase(loops_itr);
        decrement = true;
        continue;
      }


      // Check collisions, throw out if there are any
      if ( scaffold.is_collision(loops_itr->coordinates, scaffold_start, scaffold_end) ){
        if (loops_itr != results.begin()){
          loops_itr = results.erase(loops_itr);
          --loops_itr;
          continue;
        }
        else if(loops_itr == results.begin()){
          loops_itr = results.erase(loops_itr);
          decrement = true;
          continue;
        }
      }


    }

  // Prune duplicates
  cleanDuplicates();

  // Start to modify search
  change += 0.1;
  if (CA_CA + change > 49.9 || CB_CB - change < 0.1 || CB_CB + change > 49.9 || CA_CA - change < 0.1){
    return;
  }

  }

  //Sort by rmsd (lowest first)
  results.sort(rmsdSort);

  return;
}



/*
    Run's helper function. Tries to abstract away some of the lower level code
*/
void Lookup::runHelper(float CA_CA, float CB_CB, int loop_length)
{

  //std::vector<std::vector<std::vector<float> > > return_val;
  //std::vector<std::vector<char> > residues;

  //Truncate floats, multiply by 10 to get int
  char sz[64];
  double lf = CA_CA;
  sprintf(sz, "%.1lf\n", lf);
  double lf2 = atof(sz);
  CA_CA = lf2;

  lf = CB_CB;
  sprintf(sz, "%.1lf\n", lf);
  lf2 = atof(sz);
  CB_CB = lf2;

  int dCA = CA_CA * 10;
  int dCB = CB_CB * 10;


  //Open grid file to find loop pointers in loop file
  std::ifstream grid_in(database_files[2], std::fstream::binary);
  if (!grid_in){
    std::cerr << "Can't open " << database_files[2] << " to read." << std::endl;
    throw 0;
  }

  int number_of_loops;
  int record_start;
  grid_in.seekg( (80000 * dCA) + (160 * dCB) + (8 * (loop_length - 1) ) );
  grid_in.read( (char*)& number_of_loops, sizeof(int) );
  grid_in.read( (char*)& record_start, sizeof(int) );

  grid_in.close();

  //Were any loops found?
  if (number_of_loops == 0){
    //std::cout << "No loops found." << std::endl;
    return;
  }

  //Open loop file to find pdbselect pointers
  std::ifstream loop_in(database_files[1], std::fstream::binary);
  if (!loop_in){
    std::cerr << "Can't open " << database_files[1] << " to read." << std::endl;
    throw 0;
  }

  std::vector<std::vector<int> > loop_pointers;
  loop_in.seekg( record_start );

  //Read in all locations of loops and their lengths
  //This file can actually be half as big because we know the lengths already
  unsigned int pdb_position;
  std::vector<int> tmp;
  for (int i = 0; i < number_of_loops; ++i){
    loop_in.read( (char*)& pdb_position, sizeof(unsigned int) );
    loop_in.read( (char*)& loop_length, sizeof(unsigned int) );

    tmp.push_back(pdb_position);
    tmp.push_back(loop_length);
    loop_pointers.push_back(tmp);
    tmp.clear();
  }

  loop_in.close();

  //Open pdb file to get actual loops
  std::ifstream pdb_in(database_files[0], std::fstream::binary);
  if ( !pdb_in.good() ){
    std::cerr << "Can't open " << database_files[0] << " to read." << std::endl;
    throw 0;
  }

  char pdb_code[5];
  char aa;
  int residue_number;
  char chain;
  std::vector<std::vector<float> > loop;
  std::vector<char> loop_residues;
  std::vector<float> xyz;


  //Output
  for (unsigned int i = 0; i < loop_pointers.size(); ++i){
    //std::cout << "Loop at " << loop_pointers[i][0] << " of length " << loop_pointers[i][1] << std::endl;
    pdb_in.seekg(loop_pointers[i][0]);
    //std::cout << "----------------- Loop " << i + 1 << " -----------------" << std::endl;
    for (int j = 0; j < loop_pointers[i][1]; ++ j){

      pdb_in.read( pdb_code, (sizeof(char) * 4) );
      pdb_in.read( (char*)& aa, sizeof(char) );
      pdb_in.read( (char*)& residue_number, sizeof(int) );
      pdb_in.read( (char*)& chain, sizeof(char) );
      pdb_code[4] = '\0';

      //std::cout << "PDB " << pdb_code << "  Residue " << aa << residue_number << " Chain " << chain << std::endl;
      loop_residues.push_back(aa);

      for (int i = 0; i < 5; ++i){
        float x; float y; float z;
        pdb_in.read( (char*)& x, sizeof(float) );
        pdb_in.read( (char*)& y, sizeof(float) );
        pdb_in.read( (char*)& z, sizeof(float) );
        xyz.push_back(x);
        xyz.push_back(y);
        xyz.push_back(z);
        loop.push_back(xyz);
        xyz.clear();

        //std::cout << x << " " << y << " " << z << std::endl;
      }


    } //End loop

    Loop return_loop;
    return_loop.coordinates = loop;
    return_loop.sequence = loop_residues;
    results.push_back(return_loop);

    loop_residues.clear();
    loop.clear();


  }

  return;

}



/*
    Check if a given loop is similar to any of the reuslts
*/
bool Lookup::isDuplicate(const Loop &candidate)
{
  std::list<Loop>::iterator result_itr;
  for (result_itr = results.begin(); result_itr != results.end(); ++result_itr){

    // Loops already aren't the same if they have different lengths (in residues), and a loop shouldn't be compared to itself
    if (result_itr->coordinates.size() == candidate.coordinates.size() && &candidate != &(*result_itr)){
      float rmsd = RMSD(result_itr->coordinates, candidate.coordinates);
      if (rmsd < duplicate_threshold){
        return true;
      }
    }

  }
  return false;
}



/*
    Clean out duplicates using pairwise comparisons
*/
void Lookup::cleanDuplicates(){
  std::list<Loop>::iterator itr;

  for (itr = results.begin(); itr != results.end(); /*Do nothing*/ ){
    if (isDuplicate(*itr)){
      itr = results.erase(itr);
      //std::cout << "Cleaned duplicate" << std::endl;
    }
    else{
      ++itr;
    }

  }
  return;
}



/*
   Parse top-level input file
   TODO: Compress PROTEIN/DNA to MOLECULE
         PRESERVESEQUENCE
*/
bool Lookup::parse(char* input_file)
{
  std::ifstream in(input_file);
  std::string token = "";

  if (!in){
    std::string errmsg = "Can't open master input file to read. \n";
    logdump.push_back(errmsg);
    throw 0;
  }

  while(!in.eof()){
    in >> token;

    if (token == "SCAFFOLD"){
      in >> token;
      try{
        scaffold = Protein(token.c_str());
      }
      catch(const std::exception& e){
        logmsg("ERROR: Scaffold parsing failed with exception: ");
        logmsg(e.what());
        writeLog();
        exit(EXIT_FAILURE);
      }
    }

    else if (token == "ANCHORS"){
      in >> token;
      scaffold_start = atoi(token.c_str());
      in >> token;
      scaffold_end = atoi(token.c_str());

    }

    else if (token == "RANGE"){
      in >> token;
      int range_min = atoi(token.c_str());
      in >> token;
      int range_max = atoi(token.c_str());
      setRange(range_min, range_max);
    }

    else if (token == "SYMMETRY"){
      in >> token;
      symmetry = atoi(token.c_str());
    }

    else if (token == "DUPLICATECUTOFF"){
      in >> token;
      duplicate_threshold = atof(token.c_str());
    }

    else if (token == "PDBDATA"){
      in >> token;
      database_files[0] = strdup(token.c_str());
    }

    else if (token == "LOOPDATA"){
      in >> token;
      database_files[1] = strdup(token.c_str());
    }

    else if (token == "GRIDDATA"){
      in >> token;
      database_files[2] = strdup(token.c_str());
    }

    else if (token == "LOG"){
      in >> token;
      logfile = strdup(token.c_str());
    }

    else if (token == "COMPLEX"){
      while (token != "END"){
        // TODO: General "Molecule" class to read these coordinates into
        // For now, just do nothing
        in >> token;
      }
    }

  }

  logmsg("Succesfully parsed input file \n");
  return true;

}



/*
   Output lookup log file
*/
void Lookup::writeLog()
{
  std::ofstream out(logfile);
  out << "---BEGIN INDEL LOG---\n \n";
  for (unsigned int i = 0; i < logdump.size(); ++i){
    out << logdump[i] << "\n";
  }
  out << "---END INDEL LOG---\n";
}



/*
   Output lookup results, let IR daemon know how many loops were output
*/
void Lookup::iRosettaOutput()
{
  logmsg("Found " + std::to_string(results.size()) + " loops. \n");
  std::list<Loop>::iterator itr = results.begin();

  for (int i = 1; itr != results.end(); ++itr, ++i){
    // Output named so that irosetta can pick up results
    std::string fout = "loopout_" + std::to_string(i) + ".pdb";
    char* filename = strdup(fout.c_str());
    PDB_out(itr->coordinates, filename);
  }

  // Let irosetta know how many results there were
  logmsg("Succesfully wrote " + std::to_string(results.size()) + " loops. \n");
  std::cout << results.size();
  return;
}
