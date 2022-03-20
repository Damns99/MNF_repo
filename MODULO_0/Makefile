.PHONY: all

all: test

rkms_test: rkms_test.cpp runge_kutta.h
	g++ -std=c++2a -I ../shared_libs -o rkms_test rkms_test.cpp -lm -O2 # `root-config --cflags --libs`

test: rkms_test
	./rkms_test