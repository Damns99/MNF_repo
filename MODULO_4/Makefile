.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: compile

bcv_test: bcv_test.cpp bound_cond_vecs.h
	g++ -std=c++2a -I ../shared_libs -o bcv_test bcv_test.cpp -lm -O2

der_test: der_test.cpp bound_cond_vecs.h derivators.h
	g++ -std=c++2a -I ../shared_libs -o der_test der_test.cpp -lm -O2

compile: bcv_test der_test

test: bcv_test der_test
	./bcv_test
	./der_test

clean:
	rm -f bcv_test
	rm -f der_test