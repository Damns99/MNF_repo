.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: compile

double_hole: double_hole.cpp ../lattice.o ../root_lattice.o
	g++ -std=c++2a -I .. -I $(INCLUDEPATH) -L .. -Wl,-rpath=.. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o double_hole double_hole.cpp ../lattice.o ../root_lattice.o -lm -O2 $(ROOT_CONFIG) -lran2 -llattice

plot_measures: plot_measures.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_measures plot_measures.cpp -lm -O2 $(ROOT_CONFIG)

calc_all: calc_all.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_all calc_all.cpp -lm -O2 -lran2

two_point: two_point.cpp
	g++ -std=c++2a -I .. -I $(INCLUDEPATH) -L .. -Wl,-rpath=.. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o two_point two_point.cpp ../lattice.o ../root_lattice.o -lm -O2 $(ROOT_CONFIG) -lran2 -llattice

plot_two_point: plot_two_point.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_two_point plot_two_point.cpp -lm -O2 $(ROOT_CONFIG)

compile: double_hole plot_measures calc_all two_point plot_two_point

test: compile
	rm -f -r measures/test
	./double_hole 2 -nm 1000 -nc 1000 -ns 100 -len 1000 -npart 1 -snapshots 1 -folder test -beta 3 -lambda 10
	./plot_measures test
	./calc_all test double_hole_meas.txt
	./two_point test double_hole_conf 1000
	./plot_two_point test -beta 30

clean:
	rm -f double_hole
	rm -f lattice*.txt
	rm -f rand_*.txt
	rm -f *.o
	rm -f *.so
	rm -f plot_measures
	rm -f two_point