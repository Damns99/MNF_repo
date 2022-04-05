.PHONY: all clean

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

timing_test: timing_test.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o timing_test timing_test.cpp -lm -lran2

metro_ising: metro_ising.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o metro_ising metro_ising.cpp -lm -O2 $(ROOT_CONFIG) -lran2

test: timing_test metro_ising
	# ./timing_test
	./metro_ising

clean:
	rm -f timing_test
	rm -f metro_ising
	rm -f -r measures*
	rm -f lattice*.txt
	rm -f rand_*.txt