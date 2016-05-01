Loop database generation and query tools
By Will Hooper 2016

Compiling:
   Database Creation
   g++ -std=c++11 database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe -Wall

   Lookup Characterization
   g++ -std=c++11 characterize.cpp protein.cpp superimposer.cpp loop_generation.cpp -o characterize.exe -Wall

   iRosetta Lookup
   g++ -std=c++11 iRosetta_Lookup.cpp protein.cpp superimposer.cpp loop_generation.cpp -o iRosetta_Lookup.exe -Wall
