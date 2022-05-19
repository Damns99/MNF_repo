.PHONY: all clean

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs
CUDALIBPATH1 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/
CUDALIBPATH2 := /usr/local/cuda-11.4/targets/x86_64-linux/lib/stubs/

ROOT_CONFIG:= -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
NVCC_ROOT_CONFIG:= -Xcompiler -pthread -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Xcompiler \"-Wl,-rpath=/home/damns99/Programs/ROOT/builddir/lib\" -Xcompiler -pthread -lm -ldl -Xcompiler -rdynamic
#`root-config --cflags --libs`

all: test

timing_test: timing_test.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o timing_test timing_test.cpp -lm -lran2

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

test_cuda: test_cuda.cu
	/usr/local/cuda-11.4/bin/nvcc -g -lineinfo -Xcompiler -rdynamic -arch=sm_75 -gencode=arch=compute_75,code=sm_75 -std=c++17 -o test_cuda test_cuda.cu

test_cuda2: test_cuda2.cu
	/usr/local/cuda-11.4/bin/nvcc -g -lineinfo -Xcompiler -rdynamic -arch=sm_75 -gencode=arch=compute_75,code=sm_75 -std=c++17 -o test_cuda2 test_cuda2.cu

test: timing_test metro_ising cuda_metro_ising test_cuda test_cuda2
	# ./timing_test
	# ./metro_ising
	# ./cuda_metro_ising
	# ./test_cuda
	# ./test_cuda2

clean:
	rm -f timing_test
	rm -f test_cuda
	rm -f test_cuda2
	rm -f metro_ising
	rm -f cuda_metro_ising
	rm -f -r measures*
	rm -f lattice*.txt
	rm -f rand_*.txt
	rm -f *.o
	rm -f *.so