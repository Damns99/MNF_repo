.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
NVCC_ROOT_CONFIG := -Xcompiler -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Xcompiler \"-Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib\" -Xcompiler -pthread -lm -ldl -Xcompiler -rdynamic
#`root-config --cflags --libs`

BETA_STRIDE_FOLDER := beta_stride
LATTICE_LENGTH := 100
N_MEASURES := 8192

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

fit_results: fit_results.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o fit_results fit_results.cpp -lm -O2 $(ROOT_CONFIG)

test: metro_ising cuda_metro_ising calc_all plot_measures plot_results fit_results
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
	rm -f -r measures/$(BETA_STRIDE_FOLDER)/* 
	# for ii in 0.300 ; do \
	for ii in 0.400 ; do \
		./metro_ising -folder $(BETA_STRIDE_FOLDER) -len $(LATTICE_LENGTH) -nm $(N_MEASURES) -nc 100 -ns 128 -beta $$ii -geom 0 1 meas_beta_$$ii.txt; \
		./calc_all $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt $$ii -append 0; \
		./plot_measures $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt -outname beta_$$ii; \
	done
	# for ii in 0.325 0.350 0.375 0.400 0.425 0.450 0.475 0.500 ; do \
	# for ii in 0.405 0.410 0.415 0.420 0.425 0.430 0.435 0.440 0.445 0.450 0.455 0.460 0.465 0.470 0.475 0.480 0.485 0.490 0.495 0.500 ; do \
	for ii in 0.401 0.402 0.403 0.404 0.405 0.406 0.407 0.408 0.409 0.410 0.411 0.412 0.413 0.414 0.415 0.416 0.417 0.418 0.419 0.420 0.421 0.422 0.423 0.424 0.425 0.426 0.427 0.428 0.429 0.430 0.431 0.432 0.433 0.434 0.435 0.436 0.437 0.438 0.439 0.440 0.441 0.442 0.443 0.444 0.445 0.446 0.447 0.448 0.449 0.450 0.451 0.452 0.453 0.454 0.455 0.456 0.457 0.458 0.459 0.460 0.461 0.462 0.463 0.464 0.465 0.466 0.467 0.468 0.469 0.470 0.471 0.472 0.473 0.474 0.475 0.476 0.477 0.478 0.479 0.480 0.481 0.482 0.483 0.484 0.485 0.486 0.487 0.488 0.489 0.490 0.491 0.492 0.493 0.494 0.495 0.496 0.497 0.498 0.499 0.500 ; do \
		./metro_ising -folder $(BETA_STRIDE_FOLDER) -len $(LATTICE_LENGTH) -nm $(N_MEASURES) -nc 100 -ns 128 -beta $$ii -geom 0 1 meas_beta_$$ii.txt; \
		./calc_all $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt $$ii -append 1; \
		./plot_measures $(BETA_STRIDE_FOLDER) meas_beta_$$ii.txt -outname beta_$$ii; \
	done
	./plot_results $(BETA_STRIDE_FOLDER) results.txt -outname beta -xname beta