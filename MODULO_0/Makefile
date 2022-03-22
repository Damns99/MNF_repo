.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

all: test

rkms_test: rkms_test.cpp runge_kutta.h
	g++ -std=c++2a -I ../shared_libs -o rkms_test rkms_test.cpp -lm -O2 $(ROOT_CONFIG)

test: rkms_test
	./rkms_test

clean:
	rm rkms_test
	rm rk.txt
	rm rk_richardson.txt
	rm rkroot*.txt
	rm rkroot.pdf