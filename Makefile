
# Generic settings
SHELL = /bin/sh

# Targets
MYBINS = bamgrouper

# Building binaries
CC = gcc
CXX = g++
CFLAGS = -Wall -O3
CXXFLAGS = -Wall -O3
#LDENV = -Llib/
#LDFLAGS = -lbamtools
LIBS = -lz -lm

OBJECTS=BamAlignment.o \
		BamMultiReader.o \
		BamReader.o \
		BamIndex.o \
		BamWriter.o \
		BGZF.o \
		Fasta.o \
		BandedSmithWaterman.o \
		SmithWatermanGotoh.o \
		split.o

all: ogap

clean:
	rm *.o ogap

.PHONY: all clean

BamAlignment.o: BamAlignment.cpp BamAlignment.h BamAux.h
	$(CC) $(CFLAGS) -c BamAlignment.cpp 

BamMultiReader.o: BamMultiReader.cpp BamMultiReader.h BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamMultiReader.cpp 	

BamReader.o: BamReader.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamReader.cpp 	

BamWriter.o: BamWriter.cpp BamReader.h BamAux.h
	$(CC) $(CFLAGS) -c BamWriter.cpp 	

BamIndex.o: BamIndex.cpp BamIndex.h BamAux.h
	$(CC) $(CFLAGS) -c BamIndex.cpp 	

BGZF.o: BGZF.cpp BGZF.h
	$(CC) $(CFLAGS) -c BGZF.cpp

Fasta.o: Fasta.h Fasta.cpp
	$(CC) $(CFLAGS) -c Fasta.cpp

BandedSmithWaterman.o: BandedSmithWaterman.cpp BandedSmithWaterman.h
	$(CC) $(CFLAGS) -c BandedSmithWaterman.cpp

SmithWatermanGotoh.o: SmithWatermanGotoh.cpp SmithWatermanGotoh.h
	$(CC) $(CFLAGS) -c SmithWatermanGotoh.cpp

split.o: split.h split.cpp
	$(CC) $(CFLAGS) -c split.cpp

ogap: ogap.o $(OBJECTS)
	$(LDENV) $(CXX) $(CXXFLAGS) -o $@ ogap.o $(OBJECTS) $(LDFLAGS) $(LIBS)
