#include "lookup.h"

Lookup::Lookup(char* input_file)
{
  //Set default parameters
  length_range[0] = 3;   length_range[1] = 19;
  min_results = 10;      max_results = 25;
  rmsd_cutoff = 2.0;     sequence_filter = "";
  filter = false;        sequence_identity_cutoff = 0.0;
  symmetry = 1;          duplicate_threshold = 1.0;
  database_hits = 0;     scaffold_colliding_loops = 0;
  redundant_loops = 0;   complex_colliding_loops = 0;
  bad_fits = 0;          preserve_sequence = false;

  //Default filenames
  char* db1 = strdup("pdblist.dat");  database_files.push_back(db1);
  char* db2 = strdup("looplist.dat"); database_files.push_back(db2);
  char* db3 = strdup("grid.dat");     database_files.push_back(db3);
  logfile = strdup("indel_log.txt");

  //Parse input, try to find the specified anchors
  parse(input_file);

  try{
    original_loop = scaffold[0].getLoop(scaffold_start, scaffold_end);
    original_loop_anchors = collectAnchors(original_loop);
  }
  catch(int){
    logmsg("Couldn't find the requested anchor points. Check inputs and try again.\n");
    writeLog();
    exit(EXIT_FAILURE);
  }
  logmsg("Set anchors to " + std::to_string(scaffold_start) + " and " + std::to_string(scaffold_end) + "\n");

}



/*
    Set the minimum number of results to finish with
*/
void Lookup::setMin(int x){
  if (x <= 0){
     min_results = 1;
  }
  else {
    min_results = x;
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
  else if (min_length < 3 || max_length > 19) {
    logmsg("Specified loop lengths outside database range. Using default lenghts..\n");
    return;
  }
  length_range[0] = min_length;
  length_range[1] = max_length;
  logmsg("Set search for loops of length " + std::to_string(min_length) + " to " + std::to_string(max_length) + ".\n");
  return;
}



/*
    Comparison used by std::sort
*/
static bool rmsdSort(const Loop& a, const Loop& b)
{
  return (a.rmsd < b.rmsd);
}



/*
    Perform a search with set parameters
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

  // Vary the CA-CA and CB-CB distances by .1 angstroms until we have enough results or we hit the end of the database
  float change = 0.0;
  bool firstRun = true;
  while(int(results.size()) < min_results){

    //Pull raw loops from database
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
    database_hits += results_buffer.size();


    // Superimpose, record RMSD, clean out bad results
    std::list<Loop>::iterator loops_itr;
    for (loops_itr = results_buffer.begin(); loops_itr != results_buffer.end(); ++loops_itr){
      superimposeUsingAnchors(*loops_itr);
    }

    cleanBadFits(results_buffer);
    cleanCollisions(results_buffer);
    updateBuffer();
    cleanDuplicates();


    // Start to modify search
    change += 0.1;
    if (CA_CA + change > 49.9 || CB_CB - change < 0.1 || CB_CB + change > 49.9 || CA_CA - change < 0.1){
      logmsg("Exhausted the database. \n");
      break;
    }

  }

  // Output some simple statistics about the search and quit
  results.sort(rmsdSort);
  logmsg("                  Total database hits: " + std::to_string(database_hits)   + "\n");
  logmsg("        Number of self-colliding hits: " + std::to_string(scaffold_colliding_loops) + "\n");
  logmsg("     Number of complex-colliding hits: " + std::to_string(complex_colliding_loops) + "\n");
  logmsg("             Number of redundant hits: " + std::to_string(redundant_loops) + " \n");
  logmsg("Number of bad fits (high anchor RMSD): " + std::to_string(bad_fits)        + "\n");
  logmsg("               Total results returned: " + std::to_string(results.size())  + "\n");
  return;
}



/*
    Run's helper function. Abstracts away some of the ugly random access file operations
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
    std::string grid_err_msg(database_files[2]);
    logmsg("Can't open " + grid_err_msg + " to read. \n");
    writeLog();
    exit(EXIT_FAILURE);
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
    std::string loop_err_msg(database_files[1]);
    logmsg("Can't open " + loop_err_msg + " to read. \n");
    writeLog();
    exit(EXIT_FAILURE);
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
    std::string pdb_err_msg(database_files[0]);
    logmsg("Can't open " + pdb_err_msg + " to read. \n");
    writeLog();
    exit(EXIT_FAILURE);
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
    return_loop.rmsd = 0.0;
    results_buffer.push_back(return_loop);

    loop_residues.clear();
    loop.clear();


  }

  return;

}



/*
    Copy over loops that aren't filtered out into final results vector, clear buffer
*/
void Lookup::updateBuffer()
{
  std::list<Loop>::iterator itr;
  for (itr = results_buffer.begin(); itr != results_buffer.end(); ++itr){
    results.push_back(*itr);
  }
  results_buffer.clear();
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
void Lookup::cleanDuplicates()
{
  std::list<Loop>::iterator itr;
  for (itr = results.begin(); itr != results.end(); /*Do nothing*/ ){
    if (isDuplicate(*itr)){
      itr = results.erase(itr);
      ++redundant_loops;
    }
    else{
      ++itr;
    }
  }
  return;
}



/*
    Clean out loops that collide with the scaffold
*/
void Lookup::cleanCollisions(std::list<Loop>& results)
{
  std::list<Loop>::iterator itr;
  for (itr = results.begin(); itr != results.end(); /*Do nothing*/){
    if (scaffold[0].is_collision(itr->coordinates, scaffold_start, scaffold_end)){
      itr = results.erase(itr);
      ++scaffold_colliding_loops;
    }
    else if (complex.isCollision(itr->coordinates)){
      itr = results.erase(itr);
      ++complex_colliding_loops;
    }
    else{
      ++itr;
    }
  }
}


/*
    Remove (superimposed) loops if their RMSD is above the cutoff
*/
void Lookup::cleanBadFits(std::list<Loop> &results)
{
  std::list<Loop>::iterator itr;
  for (itr = results.begin(); itr != results.end(); /* Do nothing*/ ){
    if (itr->rmsd > rmsd_cutoff) {
      itr = results.erase(itr);
      ++bad_fits;
    }
    else{
      ++itr;
    }
  }
  return;
}



/*
    Collects anchor residue atom coordinates from a loop (minimum 2 residues)
    Used for anchor superimposition
*/
std::vector<std::vector<float> >
Lookup::collectAnchors(const std::vector<std::vector<float> > &loop)
{
  std::vector<std::vector<float> > anchors;
  // N terminus
  for (unsigned int i = 0; i < 5; ++i){
    anchors.push_back(loop[i]);
  }

  //C terminus
  for (unsigned int i = loop.size() - 5; i < loop.size(); ++i){
    anchors.push_back(loop[i]);
  }

  return anchors;
}



/*
    Superimposes database hit's anchors onto the scaffold target anchors and records anchor rmsd
*/
void Lookup::superimposeUsingAnchors(Loop &database_loop){
  std::vector< std::vector<float> > database_loop_anchors = collectAnchors(database_loop.coordinates);

  std::pair< std::vector< std::vector<float> > , std::vector<float> >
  transformation = superimposer(original_loop_anchors, database_loop_anchors, 10);

  for (unsigned int i = 0; i < database_loop.coordinates.size(); ++i){
    superimposer_move(database_loop.coordinates[i], transformation.first, transformation.second);
  }

  database_loop.rmsd = RMSD(original_loop_anchors, collectAnchors(database_loop.coordinates));
  return;
}



/*
   Parse top-level input file
   TODO: Compress PROTEIN/DNA to COMPLEX
         PRESERVE_SEQUENCE
*/
void Lookup::parse(char* input_file)
{
  std::ifstream in(input_file);
  std::string token = "";

  if (!in){
    logmsg("Can't open input file: " + std::string(input_file) + "\n");
    writeLog();
    exit(EXIT_FAILURE);
  }

  while(!in.eof()){
    in >> token;

    if (token == "SCAFFOLD"){
      in >> token;
      try{
        scaffold.push_back(Protein(token.c_str()));
      }
      catch(const std::exception &e){
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

    else if (token == "MIN_RESULTS"){
      in >> token;
      setMin(atoi(token.c_str()));
    }

    else if (token == "PRESERVE_SEQUENCE"){
      in >> token;
      if (token == "false" || token == "False" || token == "FALSE"){
        preserve_sequence = false;
      }
      else {
        preserve_sequence = true;
      }
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

    else if (token == "DUPLICATE_CUTOFF"){
      in >> token;
      duplicate_threshold = atof(token.c_str());
    }

    else if (token == "PDB_DATA"){
      in >> token;
      database_files[0] = strdup(token.c_str());
    }

    else if (token == "LOOP_DATA"){
      in >> token;
      database_files[1] = strdup(token.c_str());
    }

    else if (token == "GRID_DATA"){
      in >> token;
      database_files[2] = strdup(token.c_str());
    }

    else if (token == "LOG"){
      in >> token;
      logfile = strdup(token.c_str());
    }

    else if (token == "COMPLEX"){
      in >> token;
      while (token != "END"){
        try{
          complex.addMolecule(strdup(token.c_str()));
        }
        catch(const std::exception &e){
          logmsg("ERROR: Couldn't parse Complex molecule " + token + "\n");
        }
        in >> token;

      }
    }

  }

  logmsg("Succesfully parsed input file \n");
  return;

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
  out.close();
  return;
}



/*
   Output lookup results, let IR daemon know how many loops were output
*/
void Lookup::iRosettaOutput()
{
  std::list<Loop>::iterator itr = results.begin();

  for (int i = 1; itr != results.end(); ++itr, ++i){
    // Output named so that irosetta can pick up results
    std::string fout = "loopout_" + std::to_string(i) + ".pdb";
    char* filename = strdup(fout.c_str());

    if (preserve_sequence){
      PDB_out(itr->coordinates, itr->sequence, filename);
    }
    else {
      PDB_out(itr->coordinates, filename);
    }

  }

  // Let irosetta know how many results there were
  logmsg("Succesfully wrote " + std::to_string(results.size()) + " loops. \n");
  std::cout << results.size();
  return;
}
