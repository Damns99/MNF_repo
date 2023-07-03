.PHONY: all clean

ROOT_CONFIG:=-pthread -std=c++2a -m64 -I/home/damns99/Programs/ROOT/builddir/include -L/home/damns99/Programs/ROOT/builddir/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -Wl,-rpath,/home/damns99/Programs/ROOT/builddir/lib -pthread -lm -ldl -rdynamic
#`root-config --cflags --libs`

INCLUDEPATH := ../shared_libs
LIBPATH := ../shared_libs

all: compile

bcv_test: bcv_test.cpp bound_cond_vecs.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o bcv_test bcv_test.cpp -lm -O2 -Wno-deprecated

der_test: der_test.cpp bound_cond_vecs.h derivators.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o der_test der_test.cpp -lm -O2 -lmyfft -Wno-deprecated

derivator_scaling: derivator_scaling.cpp bound_cond_vecs.h derivators.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o derivator_scaling derivator_scaling.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -Wno-deprecated

int_test: int_test.cpp bound_cond_vecs.h derivators.h integrators.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o int_test int_test.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -Wno-deprecated

int_grav_waves: int_grav_waves.cpp bound_cond_vecs.h derivators.h integrators.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o int_grav_waves int_grav_waves.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -Wno-deprecated

cn_linear: cn_linear.cpp bound_cond_vecs.h derivators.h integrators.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o cn_linear cn_linear.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -ltridiag -Wno-deprecated

int_HH: int_HH.cpp bound_cond_vecs.h derivators.h integrators.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o int_HH int_HH.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -ltridiag -Wno-deprecated

int_test_2: int_test_2.cpp bound_cond_vecs.h derivators.h integrators.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o int_test_2 int_test_2.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -ltridiag -Wno-deprecated

HH_max_plot_and_fit: HH_max_plot_and_fit.cpp bound_cond_vecs.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o HH_max_plot_and_fit HH_max_plot_and_fit.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -Wno-deprecated

HH_max_plot_fitres: HH_max_plot_fitres.cpp bound_cond_vecs.h wave_plots.h
	g++ -std=c++2a -L $(LIBPATH) -Wl,-rpath=$(LIBPATH) -I $(INCLUDEPATH) -o HH_max_plot_fitres HH_max_plot_fitres.cpp -lm -O2 $(ROOT_CONFIG) -lmyfft -Wno-deprecated

compile: bcv_test der_test derivator_scaling int_test int_grav_waves cn_linear int_HH int_test_2 HH_max_plot_and_fit HH_max_plot_fitres

test: bcv_test der_test derivator_scaling int_test int_grav_waves cn_linear int_HH int_test_2 HH_max_plot_and_fit HH_max_plot_fitres
	./bcv_test
	./der_test
	./derivator_scaling
	./int_test

clean:
	rm -f bcv_test
	rm -f der_test
	rm -f derivator_scaling
	rm -f int_test
	rm -f int_grav_waves
	rm -f cn_linear

measure_HH:
	for ss in 40 35 30 25 20 15 10 5 ; do \
		for aa in -5500 -6000 -6500 -7000 -7500 -8000 -8500 -9000 ; do \
			./int_HH -sigma $$ss -ampl $$aa ; \
			./HH_max_plot_and_fit ; \
		done ; \
	done