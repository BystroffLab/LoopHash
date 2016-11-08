#Loop database generation and query tools
#By Will Hooper 2016
Dependencies: None if not using with InteractiveROSETTA, requires InteractiveROSETTA for full functionality


#------------------------------------------------


#Compilation:

   -Database Creation-
   make database

   -Lookup Characterization-
   make characterize

   -iRosetta Lookup-
   make lookup



#------------------------------------------------


#Making database files from a folder of .pdb files:
  -compile database creation files using make database
  -move [database_creation.exe, pdblist_creation.exe, make_db.sh]
  -run make_db.sh
    -outputs three files

#Querying database:
  -compile lookup using make lookup
  -place [pdblist.dat, looplist.dat, grid.dat] in same folder as iRosetta_Lookup.exe
    ./iRosetta_Lookup.exe (Protein scaffold pdb) pdblist.dat looplist.dat grid.dat (anchor start) (anchor end)
  -outputs a list of matches, do whatever you want with them


#http://www.bioinfo.rpi.edu/bystrc/
