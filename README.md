#Loop database generation and query tools
#Will Hooper 2017
Dependencies: None, requires InteractiveROSETTA for full functionality

#Compilation:

   -Database Creation-
   make database

   -Lookup Characterization-
   make characterize

   -Lookup-
   make lookup

#Making your own database from a folder of .pdb files:
  -compile database creation files
  -move [database_creation.exe, pdblist_creation.exe, make_db.sh] to a folder of pdb files you'd like to compile a database from
  -run make_db.sh
    -outputs three database files

#Querying database:
  -compile lookup
  -run lookup with an input file (see sample in scripts)
  -outputs a list of matches, do whatever you want with them (or use this program how it was intended, with InteractiveROSETTA)

#http://www.bioinfo.rpi.edu/bystrc/
