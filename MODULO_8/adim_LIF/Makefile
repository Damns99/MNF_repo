.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

INCLUDEPATH := ../../shared_libs
LIBPATH := ../../shared_libs

all: compile

box_muller.o: box_muller.h box_muller.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -Wall -pedantic -fpic box_muller.cpp

libbox_muller.so: box_muller.o
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -o libbox_muller.so box_muller.o -lm -lran2

integrate_LIF.o: integrate_LIF.h integrate_LIF.cpp
	g++ -std=c++2a -I $(INCLUDEPATH) -c -Wall -pedantic -fpic integrate_LIF.cpp -O4 -ffast-math

libintegrate_LIF.so: integrate_LIF.o libbox_muller.so
	g++ -shared -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -o libintegrate_LIF.so integrate_LIF.o -lm -O4 -ffast-math -lran2 -lbox_muller
	
integrate_LIF3: integrate_LIF3.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF3 integrate_LIF3.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated
	
integrate_LIF4: integrate_LIF4.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF4 integrate_LIF4.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated
	
integrate_LIF5: integrate_LIF5.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF5 integrate_LIF5.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated

integrate_LIF5_2: integrate_LIF5_2.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o integrate_LIF5_2 integrate_LIF5_2.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated

take_measure: take_measure.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o take_measure take_measure.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated

then_plot: then_plot.cpp libintegrate_LIF.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -L . -Wl,-rpath=. -I $(INCLUDEPATH) -o then_plot then_plot.cpp -lm -O4 $(ROOT_CONFIG) -lintegrate_LIF -ffast-math -lran2 -lbox_muller -Wno-deprecated

compile: integrate_LIF3 integrate_LIF4 integrate_LIF5 integrate_LIF5_2 take_measure then_plot

test: integrate_LIF3 integrate_LIF4 integrate_LIF5 integrate_LIF5_2 take_measure then_plot
	./integrate_LIF3
	./integrate_LIF4
	./integrate_LIF5
	./integrate_LIF5_2

clean: 