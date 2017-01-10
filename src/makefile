

CC= g++
CFLAGS= -Wall -std=c++11 -O3

all:
	$(CC) $(CFLAGS) database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe
	$(CC) $(CFLAGS) pdblist_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o pdblist_creation.exe
	$(CC) $(CFLAGS) characterize.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o characterize.exe
	$(CC) $(CFLAGS) iRosetta_Lookup.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o iRosetta_Lookup.exe

database:
	$(CC) $(CFLAGS) pdblist_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o pdblist_creation.exe
	$(CC) $(CFLAGS) database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe


lookup:
	$(CC) $(CFLAGS) iRosetta_Lookup.cpp protein.cpp superimposer.cpp lookup.cpp complex.cpp -o iRosetta_Lookup.exe

# Cleaning up executables
lookup_clean:
	rm iRosetta_Lookup.exe

db_clean:
	rm database_creation.exe
	rm pdblist_creation.exe

clean:
	rm database_creation.exe
	rm pdblist_creation.exe
	rm iRosetta_Lookup.exe
