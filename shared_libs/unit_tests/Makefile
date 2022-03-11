.PHONY: all

INCLUDEPATH := ..
LIBPATH := ..

all: test

text_io_unit_test: text_io_unit_test.cpp $(INCLUDEPATH)/text_io.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o text_io_unit_test text_io_unit_test.cpp

test: text_io_unit_test
	./text_io_unit_test