.PHONY: all clean
SHELL := /bin/bash # Use bash syntax

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG := -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
NVCC_ROOT_CONFIG := -Xcompiler -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Xcompiler \"-Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib\" -Xcompiler -pthread -lm -ldl -Xcompiler -rdynamic
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

cuda_lattice.o: lattice.h cuda_lattice.cu cuda_lattice.cuh
	/usr/local/cuda-11.4/bin/nvcc -std=c++17 -I $(INCLUDEPATH) -c -Xcompiler -fpic cuda_lattice.cu

libcuda_lattice.so: lattice.o root_lattice.o cuda_lattice.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L $(CUDALIBPATH1) -Wl,-rpath=$(CUDALIBPATH1) -L $(CUDALIBPATH2) -Wl,-rpath=$(CUDALIBPATH2) -o libcuda_lattice.so lattice.o root_lattice.o cuda_lattice.o -lm $(ROOT_CONFIG) -lran2 -lcuda -lcudart

cuda_metro_ising: cuda_metro_ising.cu liblattice.so libcuda_lattice.so
	/usr/local/cuda-11.4/bin/nvcc -arch=sm_75 -gencode=arch=compute_75,code=sm_75 -std=c++17 -I $(INCLUDEPATH) -L . -Xcompiler \"-Wl,-rpath=.\" -L $(LIBPATH) -Xcompiler \"-Wl,-rpath=$(LIBPATH)\" -o cuda_metro_ising cuda_metro_ising.cu -lm -O2 $(NVCC_ROOT_CONFIG) -lran2 -lcuda_lattice

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

retrieve_and_fit: retrieve_and_fit.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o retrieve_and_fit retrieve_and_fit.cpp -lm -O2 $(ROOT_CONFIG)

fit_critical: fit_critical.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o fit_critical fit_critical.cpp -lm -O2 $(ROOT_CONFIG)

plot_scalings: plot_scalings.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_scalings plot_scalings.cpp -lm -O2 $(ROOT_CONFIG)

retrieve_and_fit_2: retrieve_and_fit_2.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o retrieve_and_fit_2 retrieve_and_fit_2.cpp -lm -O2 $(ROOT_CONFIG)

plot_fsscaling: plot_fsscaling.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -I .. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o plot_fsscaling plot_fsscaling.cpp -lm -O2 $(ROOT_CONFIG)

test: metro_ising plot_measures calc_corr_meas calc_corr_conf calc_means calc_all retrieve_and_fit fit_critical plot_scalings cuda_metro_ising retrieve_and_fit_2 plot_fsscaling
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

analyze: plot_measures calc_corr_meas calc_corr_conf calc_means calc_all
	./plot_measures $(MEAS_DIR) -iscuda $(IS_CUDA)
	./calc_corr_meas $(MEAS_DIR) -iscuda $(IS_CUDA)
	./calc_corr_conf $(MEAS_DIR) -iscuda $(IS_CUDA)
	./calc_means $(MEAS_DIR) -iscuda $(IS_CUDA)
	./calc_all $(MEAS_DIR) -iscuda $(IS_CUDA)

measure_standard: test
	for ii in 3 35 36 37 38 39 4 405 41 415 42 425 43 435 44 445 45 46 47 48 53 ; do \
		./metro_ising -beta 0.$$ii -nm 65536 -nc 16 -ns 256 -len $(LENGTH) -snapshots 8192 -folder L$(LENGTH)b0$$ii; \
	done

analyze_standard: test
	for ii in 3 35 36 37 38 39 4 405 41 415 42 425 43 435 44 445 45 46 47 48 53 ; do \
		make analyze MEAS_DIR=L$(LENGTH)b0$$ii IS_CUDA=0; \
	done

measure_temporary: test
	for ii in 4175 4225 4275 4325 4375 4425 4475 ; do \
		./metro_ising -beta 0.$$ii -nm 16384 -nc 8 -ns 256 -len $(LENGTH) -snapshots 4096 -folder L$(LENGTH)b0$$ii; \
		make analyze MEAS_DIR=L$(LENGTH)b0$$ii IS_CUDA=0; \
	done

cuda_measure_standard: test
	for ii in 3 35 36 37 38 39 4 405 41 415 4175 42 4225 425 4275 43 4325 435 4375 44 4425 445 4475 45 46 47 48 53 ; do \
		./cuda_metro_ising -beta 0.$$ii -nm 16384 -nc 16 -ns 256 -len $(LENGTH) -snapshots 8192 -folder L$(LENGTH)b0$$ii; \
	done

cuda_analyze_standard: test
	for ii in 3 35 36 37 38 39 4 405 41 415 4175 42 4225 425 4275 43 4325 435 4375 44 4425 445 4475 45 46 47 48 53 ; do \
		make analyze MEAS_DIR=L$(LENGTH)b0$$ii IS_CUDA=1; \
	done

table_of_measures: test
	bbbb=(); \
	for ff in `ls -d measures/L*b*` ; do \
		read -a aaaa < $${ff}/metro_ising_meas.txt; \
		bbbb+=($${aaaa[15]}); \
	done; \
	u_bbbb=(`for ii in "$${bbbb[@]}"; do echo "$$ii"; done | sort -u`); \
	printf "" > measures/table_of_measures.txt ; \
	for bb in $${u_bbbb[@]} ; do \
		s_bb=$${bb//./}; \
		printf "\multicolumn{2}{c|}{$${bb}} " >> measures/table_of_measures.txt ; \
		for ll in 8 16 32 64 128; do \
			read -a cccc < measures/L$${ll}/beta_results_K_fit.txt; \
			if [ -f "measures/L$${ll}b$${s_bb}/metro_ising_meas.txt" ]; then \
				read -a aaaa < measures/L$${ll}b$${s_bb}/metro_ising_meas.txt; \
				if [ `echo "$${aaaa[15]} >= $${cccc[3]} && $${aaaa[15]} <= $${cccc[5]}" | bc` -eq 1 ]; then \
					printf "& \checkmark " >> measures/table_of_measures.txt ; \
				else \
					printf "& x " >> measures/table_of_measures.txt ; \
				fi; \
			else \
				printf "& " >> measures/table_of_measures.txt ; \
			fi; \
		done; \
		printf "\\" >> measures/table_of_measures.txt ; printf "\\" >> measures/table_of_measures.txt ; \
		printf "\n" >> measures/table_of_measures.txt ; \
	done

clean:
	rm -f metro_ising
	rm -f plot_measures
	rm -f calc_corr_meas
	rm -f calc_corr_conf
	rm -f calc_means
	rm -f calc_all
	rm -f *.o
	rm -f *.so