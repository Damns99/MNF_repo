.PHONY: all clean

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

timing_test: timing_test.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o timing_test timing_test.cpp -lm -lran2

lattice.o: lattice.h lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic lattice.cpp

root_lattice.o: lattice.h root_lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic root_lattice.cpp $(ROOT_CONFIG)

cuda_lattice.o: lattice.h cuda_lattice.cu
	/usr/local/cuda-11.4/bin/nvcc -std=c++17 -I $(INCLUDEPATH) -c -Xcompiler -fpic cuda_lattice.cu

liblattice.so: lattice.o root_lattice.o cuda_lattice.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L $(CUDALIBPATH1) -Wl,-rpath=$(CUDALIBPATH1) -L $(CUDALIBPATH2) -Wl,-rpath=$(CUDALIBPATH2) -o liblattice.so lattice.o root_lattice.o cuda_lattice.o -lm $(ROOT_CONFIG) -lran2 -lcuda -lcudart

metro_ising: metro_ising.cpp liblattice.so
	g++ -std=c++2a -I $(INCLUDEPATH) -L . -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o metro_ising metro_ising.cpp -lm -O2 $(ROOT_CONFIG) -llattice -lran2

test: timing_test metro_ising
	# ./timing_test
	./metro_ising

clean:
	rm -f timing_test
	rm -f metro_ising
	rm -f -r measures*
	rm -f lattice*.txt
	rm -f rand_*.txt