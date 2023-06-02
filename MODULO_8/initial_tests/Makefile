.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

INCLUDEPATH := ../shared_libs
LIBPATH := ../shared_libs

all: compile

integrate_LIF.o: integrate_LIF.h integrate_LIF.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -Wall -pedantic -fpic integrate_LIF.cpp -O4 -ffast-math

libintegrate_LIF.so: integrate_LIF.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o libintegrate_LIF.so integrate_LIF.o -lm -O4 -ffast-math

integrate_LIF1: integrate_LIF1.cpp
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF1 integrate_LIF1.cpp -lm -O4 $(ROOT_CONFIG)

integrate_LIF2: integrate_LIF2.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF2 integrate_LIF2.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF
	
integrate_LIF3: integrate_LIF3.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF3 integrate_LIF3.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math
	
integrate_LIF4: integrate_LIF4.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF4 integrate_LIF4.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math
	
integrate_LIF5: integrate_LIF5.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF5 integrate_LIF5.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math

integrate_sin: integrate_sin.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_sin integrate_sin.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math

compile: integrate_LIF1 integrate_LIF2 integrate_LIF3 integrate_LIF4 integrate_LIF5 integrate_sin

test: integrate_LIF1 integrate_LIF2 integrate_LIF3 integrate_LIF4 integrate_LIF5 integrate_sin
	./integrate_LIF1
	./integrate_LIF2
	./integrate_LIF3
	./integrate_LIF4
	./integrate_LIF5
	./integrate_sin

clean: 