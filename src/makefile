

CC=g++
CFLAGS= -Wall -std=c++11

all:
	$(CC) $(CFLAGS) database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe
	$(CC) $(CFLAGS) characterize.cpp protein.cpp superimposer.cpp loop_generation.cpp -o characterize.exe
	$(CC) $(CFLAGS) iRosetta_Lookup.cpp protein.cpp superimposer.cpp loop_generation.cpp -o iRosetta_Lookup.exe

database:
	$(CC) $(CFLAGS) database_creation.cpp protein.cpp lookup.cpp superimposer.cpp loop_generation.cpp -o database_creation.exe

characterize:
	$(CC) $(CFLAGS) characterize.cpp protein.cpp superimposer.cpp loop_generation.cpp -o characterize.exe

lookup:
	$(CC) $(CFLAGS) iRosetta_Lookup.cpp protein.cpp superimposer.cpp lookup.cpp loop_generation.cpp -o iRosetta_Lookup.exe
