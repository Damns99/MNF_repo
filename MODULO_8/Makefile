.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

INCLUDEPATH := ../shared_libs
LIBPATH := ../shared_libs

all: compile

integrate_LIF: integrate_LIF.cpp
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o integrate_LIF integrate_LIF.cpp -lm -O2 $(ROOT_CONFIG)

integrate_LIF2: integrate_LIF2.cpp
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o integrate_LIF2 integrate_LIF2.cpp -lm -O2 $(ROOT_CONFIG)

compile: integrate_LIF integrate_LIF2

test: integrate_LIF integrate_LIF2
	./integrate_LIF
	./integrate_LIF2

clean: 