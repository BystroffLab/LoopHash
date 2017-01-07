#Loop database generation and query tools
Will Hooper 2017 \n
Dependencies: None, requires InteractiveROSETTA for full functionality

#Compilation:

   -Database Creation- \n
   make database \n

   -Lookup Characterization- \n
   make characterize \n

   -Lookup- \n
   make lookup \n

#Making your own database from a folder of .pdb files:
  -compile database creation files \n
  -move [database_creation.exe, pdblist_creation.exe, make_db.sh] to a folder of pdb files you'd like to compile a database from \n
  -run make_db.sh \n
    -outputs three database files \n

#Querying database:
  -compile lookup \n
  -run lookup with an input file (see sample in scripts) \n
  -outputs a list of matches, do whatever you want with them (or use this program how it was intended, with InteractiveROSETTA) \n

#http://www.bioinfo.rpi.edu/bystrc/
