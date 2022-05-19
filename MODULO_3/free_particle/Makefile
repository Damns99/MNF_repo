.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

free_particle: free_particle.cpp
	g++ -std=c++2a -I .. -I $(INCLUDEPATH) -L .. -Wl,-rpath=.. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o free_particle free_particle.cpp ../lattice.o ../root_lattice.o -lm -O2 $(ROOT_CONFIG) -lran2 -llattice


test: free_particle
	#./free_particle

clean:
	rm -f free_particle
	rm -f lattice*.txt
	rm -f rand_*.txt
	rm -f *.o
	rm -f *.so