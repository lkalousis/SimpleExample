
#
# makefile
#
# author: Leonidas N. Kalousis
#

# files

srcs = example.cc

objs = $(srcs:.cc=.o)

exe = $(srcs:.cc=)-exec

# libs

cxx = g++

cxxflags = -g -W -O -Wall -Wno-deprecated -fPIC -std=c++11

cxxflags += -I$(shell root-config --prefix)/include/ -I$(shell gsl-config --prefix)/include/gsl/

rootlibs = $(shell root-config --glibs) -lMinuit -lMinuit2 -lfftw3 -lm -L/$(shell gsl-config --prefix)/ -lgsl -lgslcblas 

# targets

all	: start $(exe) end

start	: 
	@echo ''		
	@echo ' * make ...'
	@echo ''
	@rm -f ./#* ./*~ ./*.*~	

$(exe)	: $(objs) 
	@$(cxx) $(rootlibs) $(cxxflags) -o $@ $^
	@echo ' * building ( exe )  :' $@
	@echo ''

%.o	: %.cc	
	@$(cxx) $(cxxflags) -c -o $@ $<
	@echo ' * building ( objs ) :' $@
	@echo ''

end	:
	@echo ' * make done !'
	@echo ''

clean	:	
	@echo ''	
	@echo ' * make clean ...'
	@echo ''		
	@rm -f ./#* ./*~ ./*.*~ *.o $(exe)

fresh	: clean all
