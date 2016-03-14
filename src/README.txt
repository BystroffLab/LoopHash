Loop database generation and query tools
By Will Hooper 2016

Compile using:
   g++ -std=c++11 main.cpp protein.cpp superimposer.cpp loop_generation.cpp -o test.exe -Wall

Database generation:
-Master PDB random access file generation:
  -Parse PDB
  -Write PDB

-Indexing of PDB file:
  -Count loops
  -Write loops

Database querying:
  -Query 
