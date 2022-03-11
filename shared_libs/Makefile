.PHONY: all

all: libran2.so libmyfft.so libtridiag.so

ran2.o: ran2.cpp ran2.h
	g++ -std=c++2a -c -Wall -pedantic -fpic ran2.cpp

libran2.so: ran2.o
	g++ -shared -o libran2.so ran2.o

myfft.o: myfft.cpp myfft.h
	g++ -std=c++2a -c -Wall -pedantic -fpic myfft.cpp

libmyfft.so: myfft.o
	g++ -shared -o libmyfft.so myfft.o -lfftw3 -lm

tridiag.o: tridiag.cpp tridiag.h
	g++ -std=c++2a -c -Wall -pedantic -fpic tridiag.cpp

libtridiag.so: tridiag.o
	g++ -shared -o libtridiag.so tridiag.o
