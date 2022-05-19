.PHONY: all clean

INCLUDEPATH := ../shared_libs
LIBPATH := ../shared_libs

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

bootstrap_test: bootstrap_test.cpp bootstrap.h libbox_muller.so
	g++ -std=c++2a -I $(INCLUDEPATH) -L . -Wl,-rpath=. -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o bootstrap_test bootstrap_test.cpp -lm -O2 $(ROOT_CONFIG) -lbox_muller -lran2

box_muller.o: box_muller.h box_muller.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -Wall -pedantic -fpic box_muller.cpp

libbox_muller.so: box_muller.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o libbox_muller.so box_muller.o -lm -lran2

test: bootstrap_test
	./bootstrap_test

clean:
	rm -f bootstrap_test
	rm -f bootstrap*.pdf