.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

rk_test: rk_test.cpp runge_kutta.h
	g++ -std=c++2a -I ../shared_libs -o rk_test rk_test.cpp -lm -O2 $(ROOT_CONFIG)

ms_test: ms_test.cpp multistep.h runge_kutta.h
	g++ -std=c++2a -I ../shared_libs -o ms_test ms_test.cpp -lm -O2 $(ROOT_CONFIG)

comparison: comparison.cpp  multistep.h runge_kutta.h
	g++ -std=c++2a -I ../shared_libs -o comparison comparison.cpp -lm -O2 $(ROOT_CONFIG)

compile: rk_test ms_test comparison

test: rk_test ms_test comparison
	./rk_test
	./ms_test
	./comparison

clean:
	rm -f rk_test
	rm -f rk*txt
	rm -f rkroot.pdf
	rm -f ms_test
	rm -f ms*.txt
	rm -f msroot.pdf
	rm -f comparison
	rm -f comparison.pdf