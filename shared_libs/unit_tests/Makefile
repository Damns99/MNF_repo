.PHONY: all

INCLUDEPATH := ..
LIBPATH := ..

all: test

text_io_unit_test: text_io_unit_test.cpp $(INCLUDEPATH)/text_io.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o text_io_unit_test text_io_unit_test.cpp

cmdline_parser_unit_test: cmdline_parser_unit_test.cpp $(INCLUDEPATH)/cmdline_parser.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o cmdline_parser_unit_test cmdline_parser_unit_test.cpp

test: text_io_unit_test cmdline_parser_unit_test
	./text_io_unit_test
	./cmdline_parser_unit_test 31 -b 23.456 23.456 23.456 -c c