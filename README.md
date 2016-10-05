#Loop database generation and query tools
#By Will Hooper 2016
Dependencies: None if not using with InteractiveROSETTA, requires InteractiveROSETTA and PyRosetta packages if so


------------------------------------------------


#Compilation:

   -Database Creation-
   make database
        or
   g++ -std=c++11 database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe -Wall
   g++ -std=c++11 pdblist_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o pdblist_creation.exe -Wall


   -Lookup Characterization-
   make characterize
        or
   g++ -std=c++11 characterize.cpp protein.cpp superimposer.cpp loop_generation.cpp -o characterize.exe -Wall

   -iRosetta Lookup-
   make lookup
        or
   g++ -std=c++11 iRosetta_Lookup.cpp protein.cpp superimposer.cpp lookup.cpp loop_generation.cpp -o iRosetta_Lookup.exe -Wall


------------------------------------------------


#Making database files from a folder of .pdb files:
  -compile database creation files using make database
  -move [database_creation.exe, pdblist_creation.exe, make_db.sh]
  -run make_db.sh
    -outputs three files

#Querying database:
  -compile lookup using make lookup
  -place [pdblist.dat, looplist.dat, grid.dat] in same folder as iRosetta_Lookup.exe
    ./iRosetta_Lookup.exe (Protein scaffold pdb) pdblist.dat looplist.dat grid.dat (anchor start) (anchor end)
  -outputs a list of matches, do whatever you want with them, anything! I don't care.


http://www.bioinfo.rpi.edu/bystrc/
