//Script to characterize loop search
#include "protein.h"
#include "loop_generation.h"
#include <random>
#include <chrono>

std::vector<std::vector<float> > chooseRandomLoop(Protein protein, int& start, int& end ){

  //Initialize random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> start_distribution(1,300);
  std::uniform_int_distribution<int> length_distribution(3,19);

  std::vector<std::vector<float> > anchor_loop;
  while (true){

    //Get random numbers start_rn and length_rn, use to choose anchors
    int start_rn = start_distribution(gen);
    int length_rn = length_distribution(gen);
    std::cout << start_rn << " " << length_rn << std::endl;


    //Try to get anchors, first make sure loop is within protein, then check that
    //anchors are within 15A of each other.
    try{
      anchor_loop = protein.getLoop(start_rn, (start_rn + length_rn) );
    }
    catch(int){
      //std::cout << "Loop with start " << start_rn << " and end " << start_rn + length_rn << " was not contained in protein." << std::endl;
      continue;
    }

    if( ca_ca_dist(anchor_loop) > 15.1 && cb_cb_dist(anchor_loop) > 15.1 ){
      anchor_loop.clear();
      continue;
    }
    else{
      start = start_rn;
      end = start_rn + length_rn - 1;
      break;
    }

  }


  return anchor_loop;

}


std::vector<std::pair<int, float> >
compareAnchors(std::vector<std::vector<std::vector<float> > >& loop_list, std::vector<std::vector<float> > original_loop){

  std::vector<std::pair<int, float> > output;

  //for each loop in loop list, superpose its first and last residues onto the original
  //loop's anchors, record <length, rmsd> in output vector

  for (unsigned int i = 0; i < loop_list.size(); ++i){
    std::vector<std::vector<float> > list_superposer_in;
    std::vector<std::vector<float> > original_superposer_in;


    //Get first residue for loop list
    for (unsigned int j = 0; j < 5; ++j){
      list_superposer_in.push_back(loop_list[i][j]);
    }
    //Get last residue for loop list (can't do both at once because loops not always same size)
    for (unsigned int j = loop_list[i].size() - 5; j < loop_list[i].size() ; ++j){
      list_superposer_in.push_back(loop_list[i][j]);
    }

    //Get first residue for original
    for (unsigned int j = 0; j < 5; ++j){
      original_superposer_in.push_back(original_loop[j]);
    }
    //Get last residue for original
    for (unsigned int j = original_loop.size() - 5; j < original_loop.size() ; ++j){
      original_superposer_in.push_back(original_loop[j]);
    }

    //Sanity check
    assert(list_superposer_in.size() == 10);
    assert(original_superposer_in.size() == 10);

    //Calculate superposition, superimpose list loop on original loop, get rmsd
    //std::cout << "RMSD before:: " << RMSD(original_superposer_in, list_superposer_in) << std::endl;

    std::pair< std::vector< std::vector<float> > , std::vector<float> > transformation;
    transformation = superimposer(original_superposer_in, list_superposer_in, 10);
    for (unsigned int j = 0; j < list_superposer_in.size(); ++j){
      superimposer_move(list_superposer_in[j], transformation.first, transformation.second);
    }
    for (unsigned int j = 0; j < loop_list[i].size(); ++j){
      superimposer_move(loop_list[i][j], transformation.first, transformation.second);
    }


    float rmsd = RMSD(original_superposer_in, list_superposer_in);
    int len = loop_list[i].size()/5;
    //If we get a crazy high RMSD, write out template and insertion to the same pdb file for inspection
    if (rmsd > 5){
      std::cout << "RMSD of " << rmsd << " when trying to insert " << len << std::endl;
      //Use int of rmsd to generate filename
      int i_rmsd = rmsd * 100;
      std::stringstream ss;
      ss << i_rmsd;
      std::string s_rmsd = ss.str();
      s_rmsd = s_rmsd + ".pdb";
      char* filename = const_cast<char*>( s_rmsd.c_str() );
      PDB_out(loop_list[i], filename);
      PDB_out(original_loop, filename);
    }
    if (rmsd < 1.3){
      //std::cout << "RMSD of " << rmsd << " when trying to insert " << len << std::endl;
      //Use int of rmsd to generate filename
      int i_rmsd = rmsd * 100;
      std::stringstream ss;
      ss << i_rmsd;
      std::string s_rmsd = ss.str();
      s_rmsd = s_rmsd + ".pdb";
      char* filename = const_cast<char*>( s_rmsd.c_str() );
      //PDB_out(loop_list[i], filename);
      //PDB_out(original_loop, filename);
    }
    //std::cout << "RMSD after:: " << rmsd << std::endl;


    output.push_back( std::make_pair(len, rmsd) ) ;
  }


  return output;
}


//Outputs .csv in the form Loop Length,Anchor RMSD
void excelOutput(char* outfile, std::vector<std::pair<int, float> > statistics){

  //Try to open file for output (append onto end of file)
  std::ofstream out_file(outfile, std::fstream::app);
  if ( !out_file.good() ) {
    std::cerr << "Can't open " << out_file << " to write." << std::endl;
    exit(1);
  }

  //Write out data from vector of <int,float> pairs
  for (unsigned int i = 0; i < statistics.size(); ++i){
    out_file << statistics[i].first << ',' << statistics[i].second << "\n";
  }

  return;
}

//Outputs .csv in the form Loop Length,Anchor RMSD,Sequence Separation
void excelOutput(char* outfile, std::vector<std::pair<int, float> > statistics, int sequence_separation){

  //Try to open file for output (append onto end of file)
  std::ofstream out_file(outfile, std::fstream::app);
  if ( !out_file.good() ) {
    std::cerr << "Can't open " << out_file << " to write." << std::endl;
    exit(1);
  }

  //Write out data from vector of <int,float> pairs
  for (unsigned int i = 0; i < statistics.size(); ++i){
    out_file << statistics[i].first << ',' << statistics[i].second << ',' << sequence_separation << "\n";
  }

  return;
}

//----------------------------------------------------------
//----------------------------------------------------------


int main(int argc, char* argv[]){

  // ./characterize.exe xxxx.pdb output.txt int(# anchor sets to try)
  assert(argc == 4);

  //Create protein object to operate on
  Protein test_set(argv[1]);

  //Start loop iterating through anchor sets
  int anchor_sets = atoi(argv[3]);

  //Track timing
  std::vector<std::pair<int, float> > timing;

  for(int i = 0; i < anchor_sets; ++i){
    //Choose random loop takes 2 int references, updates when it's chosen a suitable loop
    //end -start + 1 == sequence separation
    int start = 0; int end = 0;
    std::vector<std::vector<float> > myLoop = chooseRandomLoop(test_set, start, end);
    int sequence_separation = end - start + 1;
    //--------------------------------------------------------
    //Query database for all loops 2 < L < 20 given dCA and dCB
    //--------------------------------------------------------
    auto time_start = std::chrono::steady_clock::now();
    std::vector<std::vector<std::vector<float> > > query_return;
    query_return = continuous_query_wrapper(ca_ca_dist(myLoop), cb_cb_dist(myLoop), "pdbselect.prot", "dbv2.loop", "gridv2.arr");

    //Superpose and get statistics, then check for collisions
    //If a loop doesn't collide, push its statistics to an output vector
    std::vector<std::pair<int, float> > statistics = compareAnchors(query_return, myLoop);

    std::vector<std::pair<int, float> > output;
    int counter = 0;
    for(unsigned int i = 0; i < query_return.size(); ++i){
      if (test_set.is_collision(query_return[i], start, end) ){
        ++counter;
      }
      else{
        output.push_back(statistics[i]);
      }


    }

    //--------------------------------------------------------
    //End of query
    //--------------------------------------------------------

    std::cout << "Loops found: " << query_return.size() << std::endl;
    std::cout << "Loops with collisions: " << counter << std::endl;
    std::cout << "Start: " << start << "End: " << end << std::endl;
    auto time_end = std::chrono::steady_clock::now();
    auto difference = time_end - time_start;
    //timing.push_back(std::make_pair( query_return.size(), std::chrono::duration<double, std::milli>(difference).count() ));
    //std::cout << std::chrono::duration<double, std::milli>(difference).count() << " ms" << std::endl;
    //excelOutput(argv[2] , statistics, sequence_separation);
    timing.clear();

  }


  //std::cout << "Finished test." << std::endl << std::endl;
  return 0;
}
