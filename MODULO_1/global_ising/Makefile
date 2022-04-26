.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
NVCC_ROOT_CONFIG := -Xcompiler -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Xcompiler \"-Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib\" -Xcompiler -pthread -lm -ldl -Xcompiler -rdynamic
#`root-config --cflags --libs`

BETA_STRIDE_FOLDER := _beta_stride

all: test

lattice.o: lattice.h lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic lattice.cpp

root_lattice.o: lattice.h root_lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic root_lattice.cpp $(ROOT_CONFIG)

cuda_lattice.o: lattice.h cuda_lattice.cu cuda_lattice.cuh
	/usr/local/cuda-11.4/bin/nvcc -std=c++17 -I $(INCLUDEPATH) -c -Xcompiler -fpic cuda_lattice.cu

liblattice.so: lattice.o root_lattice.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o liblattice.so lattice.o root_lattice.o -lm $(ROOT_CONFIG) -lran2

libcuda_lattice.so: lattice.o root_lattice.o cuda_lattice.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L $(CUDALIBPATH1) -Wl,-rpath=$(CUDALIBPATH1) -L $(CUDALIBPATH2) -Wl,-rpath=$(CUDALIBPATH2) -o libcuda_lattice.so lattice.o root_lattice.o cuda_lattice.o -lm $(ROOT_CONFIG) -lran2 -lcuda -lcudart

metro_ising: metro_ising.cpp liblattice.so
	g++ -std=c++2a -I $(INCLUDEPATH) -L . -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o metro_ising metro_ising.cpp lattice.o root_lattice.o -lm -O2 $(ROOT_CONFIG) -lran2 -llattice

cuda_metro_ising: cuda_metro_ising.cu libcuda_lattice.so
	/usr/local/cuda-11.4/bin/nvcc -arch=sm_75 -gencode=arch=compute_75,code=sm_75 -std=c++17 -I $(INCLUDEPATH) -L . -Xcompiler \"-Wl,-rpath=.\" -L $(LIBPATH) -Xcompiler \"-Wl,-rpath=$(LIBPATH)\" -o cuda_metro_ising cuda_metro_ising.cu -lm -O2 $(NVCC_ROOT_CONFIG) -lran2 -lcuda_lattice

calc_all: calc_all.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_all calc_all.cpp -lm -O2 -lran2

plot_measures: plot_measures.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_measures plot_measures.cpp -lm -O2 $(ROOT_CONFIG)

plot_results: plot_results.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_results plot_results.cpp -lm -O2 $(ROOT_CONFIG)

test: metro_ising cuda_metro_ising calc_all plot_measures plot_results
	# ./metro_ising
	# ./cuda_metro_ising

clean:
	rm -f metro_ising
	rm -f cuda_metro_ising
	rm -f -r measures2*
	rm -f lattice*.txt
	rm -f rand_*.txt
	rm -f *.o
	rm -f *.so
	rm -f calc_all
	rm -f plot_measures
	rm -f plot_results

beta_stride: metro_ising calc_all plot_measures plot_results
	rm -f -r measures$(BETA_STRIDE_FOLDER)/* 
	for ii in 0.300 ; do \
		./metro_ising -folder $(BETA_STRIDE_FOLDER) -len 100 -nm 1024 -nc 100 -ns 128 -beta $$ii 1 meas_beta_$$ii.txt; \
		./calc_all $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt $$ii -append 0; \
		./plot_measures $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt -outname beta_$$ii; \
	done
	for ii in 0.325 0.350 0.375 0.400 0.425 0.450 0.475 0.500 ; do \
		./metro_ising -folder $(BETA_STRIDE_FOLDER) -len 100 -nm 1024 -nc 100 -ns 128 -beta $$ii 1 meas_beta_$$ii.txt; \
		./calc_all $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt $$ii -append 1; \
		./plot_measures $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt -outname beta_$$ii; \
	done
	./plot_results $(BETA_STRIDE_FOLDER) results.txt -outname beta -xname beta