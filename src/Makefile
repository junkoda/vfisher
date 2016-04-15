#
# Fisher matrix code for peculiar velocity surveys
#

#CXX        = g++ # set C++ compiler if necessary
CXXFLAGS  := -O3 -Wall
LIBS      := -lm -lboost_program_options -lgsl -lgslcblas

# In some environment BOOST library name is -lboost_program_options-mt

# Location of libraries BOOST/program_options and GSL
# /include and /lib are added below

BOOST_DIR ?= #e.g. /opt/local 
GSL_DIR   ?= #e.g. $(HOME)/Research/opt/gcc/gsl
DIR_PATH   = $(BOOST_DIR) $(GSL_DIR)

CXXFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

EXEC    = vfisher
all: $(EXEC)

OBJS = fisher_matrix.o volume.o print_fisher_matrix.o power_spectrum_camb.o read_prior.o


vfisher: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

# g++ -MM *.cpp
fisher_matrix.o: fisher_matrix.cpp fisher_matrix.h sigma_redshift_space.h \
  power_spectrum_camb.h volume.h print_fisher_matrix.h read_prior.h
power_spectrum_camb.o: power_spectrum_camb.cpp fisher_matrix.h \
  power_spectrum_camb.h
print_fisher_matrix.o: print_fisher_matrix.cpp fisher_matrix.h \
  print_fisher_matrix.h
read_prior.o: read_prior.cpp read_prior.h
volume.o: volume.cpp volume.h fisher_matrix.h


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)

