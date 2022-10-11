.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

lattice.o: lattice.h lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic lattice.cpp

root_lattice.o: lattice.h root_lattice.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -fpic root_lattice.cpp $(ROOT_CONFIG)

liblattice.so: lattice.o root_lattice.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o liblattice.so lattice.o root_lattice.o -lm $(ROOT_CONFIG) -lran2

metro_ising: metro_ising.cpp liblattice.so
	g++ -std=c++2a -I $(INCLUDEPATH) -L . -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o metro_ising metro_ising.cpp lattice.o root_lattice.o -lm -O2 $(ROOT_CONFIG) -lran2 -llattice

plot_measures: plot_measures.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_measures plot_measures.cpp -lm -O2 $(ROOT_CONFIG)

calc_corr_meas: calc_corr_meas.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_corr_meas calc_corr_meas.cpp -lm -O2 $(ROOT_CONFIG)

calc_corr_conf: calc_corr_conf.cpp liblattice.so
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L . -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_corr_conf calc_corr_conf.cpp -lm -O2 $(ROOT_CONFIG) -lran2 -llattice

calc_means: calc_means.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_means calc_means.cpp -lm -O2 $(ROOT_CONFIG)

calc_all: calc_all.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o calc_all calc_all.cpp -lm -O2 -lran2

test: metro_ising plot_measures calc_corr_meas calc_corr_conf calc_means calc_all
	# ./metro_ising

plot_all_measures: plot_measures
	cd measures; \
	for d in *; do \
		if [ -f "$$d/metro_ising_meas.txt" ]; then \
			#if [ ! -f "$$d/m_hist.pdf" ]; then \
				cd ..; \
				./plot_measures $$d; \
				cd measures; \
			#fi; \
		fi; \
	done; \
	cd ..;

calc_corr_meas_all: calc_corr_meas
	cd measures; \
	for d in *; do \
		if [ -f "$$d/metro_ising_meas.txt" ]; then \
			#if [ ! -f "$$d/corr.txt" ]; then \
				cd ..; \
				./calc_corr_meas $$d; \
				cd measures; \
			#fi; \
		fi; \
	done; \
	cd ..;

calc_corr_conf_all: calc_corr_conf
	cd measures; \
	for d in *; do \
		if [ -f "$$d/metro_ising_conf.txt" ]; then \
			#if [ ! -f "$$d/spatial_corr.pdf" ]; then \
				cd ..; \
				./calc_corr_conf $$d; \
				cd measures; \
			#fi; \
		fi; \
	done; \
	cd ..;

calc_means_all: calc_means
	cd measures; \
	for d in *; do \
		if [ -f "$$d/metro_ising_meas.txt" ]; then \
			#if [ ! -f "$$d/m_blocking.pdf" ]; then \
				cd ..; \
				./calc_means $$d; \
				cd measures; \
			#fi; \
		fi; \
	done; \
	cd ..;

calc_all_all: calc_all
	cd measures; \
	for d in *; do \
		if [ -f "$$d/metro_ising_meas.txt" ]; then \
			#if [ ! -f "$$d/results.txt" ]; then \
				cd ..; \
				./calc_all $$d; \
				cd measures; \
			#fi; \
		fi; \
	done; \
	cd ..;

clean:
	rm -f metro_ising
	rm -f plot_measures
	rm -f *.o
	rm -f *.so