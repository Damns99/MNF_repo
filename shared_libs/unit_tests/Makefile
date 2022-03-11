.PHONY: all

INCLUDEPATH := ..
LIBPATH := ..

all: test

$(LIBPATH)/libran2.so: $(INCLUDEPATH)/ran2.cpp $(INCLUDEPATH)/ran2.h
	cd $(LIBPATH) && make

$(LIBPATH)/libmyfft.so: $(INCLUDEPATH)/myfft.cpp $(INCLUDEPATH)/myfft.h
	cd $(LIBPATH) && make

text_io_unit_test: text_io_unit_test.cpp $(INCLUDEPATH)/text_io.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o text_io_unit_test text_io_unit_test.cpp

cmdline_parser_unit_test: cmdline_parser_unit_test.cpp $(INCLUDEPATH)/cmdline_parser.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o cmdline_parser_unit_test cmdline_parser_unit_test.cpp

vec_hist_unit_test: vec_hist_unit_test.cpp $(INCLUDEPATH)/vec_hist.h
	g++ -std=c++2a -I $(INCLUDEPATH) -o vec_hist_unit_test vec_hist_unit_test.cpp

ran2_unit_test: ran2_unit_test.cpp $(LIBPATH)/libran2.so $(INCLUDEPATH)/vec_hist.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o ran2_unit_test ran2_unit_test.cpp -lran2

myfft_unit_test: myfft_unit_test.cpp $(LIBPATH)/libmyfft.so
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o myfft_unit_test myfft_unit_test.cpp -lmyfft -lm

test: text_io_unit_test cmdline_parser_unit_test vec_hist_unit_test ran2_unit_test myfft_unit_test
	./text_io_unit_test
	./cmdline_parser_unit_test 31 -b 23.456 23.456 23.456 -c c
	./vec_hist_unit_test
	./ran2_unit_test
	./myfft_unit_test