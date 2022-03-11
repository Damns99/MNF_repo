.PHONY: all

all: libran2.so

ran2.o: ran2.cpp ran2.h
	g++ -std=c++2a -c -Wall -pedantic -fpic ran2.cpp

libran2.so: ran2.o
	g++ -shared -o libran2.so ran2.o
