
# Generic settings
SHELL = /bin/sh

# Targets
MYBINS = bamgrouper

# Building binaries
CC = gcc
CXX = g++
CFLAGS = -Wall -O3 -D_FILE_OFFSET_BITS=64
#CXXFLAGS = -Wall -O3
#LDENV = -Llib/
#LDFLAGS = -lbamtools
LIBS = -lz -lm -lbamtools

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib

CXXFLAGS=-Wall -O3 -D_FILE_OFFSET_BITS=64

SMITHWATERMAN = smithwaterman/SmithWatermanGotoh.o
REPEATS = smithwaterman/Repeats.o
INDELALLELE = smithwaterman/IndelAllele.o
DISORDER = smithwaterman/disorder.c
LEFTALIGN = smithwaterman/LeftAlign.o
FASTAHACK = fastahack/Fasta.o

OBJECTS=split.o

all: ogap

clean:
	rm -f ogap *.o
	cd fastahack && $(MAKE) clean
	cd bamtools/build && $(MAKE) clean
	cd smithwaterman && $(MAKE) clean

.PHONY: all clean



# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && make
	cp bamtools/lib/libbamtools.a ./


$(SMITHWATERMAN):
	cd smithwaterman && $(MAKE)

$(DISORDER): $(SMITHWATERMAN)

$(REPEATS): $(SMITHWATERMAN)

$(LEFTALIGN): $(SMITHWATERMAN)

$(INDELALLELE): $(SMITHWATERMAN)

$(FASTAHACK):
	cd fastahack && $(MAKE)

split.o: split.h split.cpp
	$(CXX) $(CFLAGS) -c split.cpp

ogap: ogap.cpp $(OBJECTS) $(SMITHWATERMAN) $(FASTAHACK) $(DISORDER) $(REPEATS) $(LEFTALIGN) $(INDELALLELE) libbamtools.a
	$(CXX) $(CXXFLAGS) -o $@ ogap.cpp $(OBJECTS) $(LDFLAGS) $(SMITHWATERMAN) $(FASTAHACK) $(DISORDER) $(REPEATS) $(LEFTALIGN) $(INDELALLELE) $(LIBS) 
