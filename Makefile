# COMPILER

#IBM cluster compiler and edu215 compiler
#CXX= g++
CXX= mpicxx

#Cygwin compilers
#CXX= /usr/bin/g++ (this was also used on edu215)
#CXX= c:/MinGW/bin/g++
#CXX= c:/MinGW/bin/c++
#CXX= c:/MinGW/bin/gcc
#CXX= c:/MinGW/bin/mingw32-g++
#CXX= c:/MinGW/bin/mingw32-c++
#CXX= c:/MinGW/bin/mingw32-gcc

# FILES

OBJS =	main.o \
	BB_FORM_EXT.o \
	global_variables.o \
	global_functions.o \
	graph.o \
	digraph.o \
	solution.o 

# CPLEX VERSION (LIBS and INCLUDE files)

#cplex 12.6 HERE
CPLEXLIBDIR =  /opt/ibm/ILOG/CPLEX_Studio126/cplex/lib/x86-64_linux/static_pic
CPXLP_INCLUDE= /opt/ibm/ILOG/CPLEX_Studio126/cplex/include
CONCERTLIBDIR =  /opt/ibm/ILOG/CPLEX_Studio126/concert/lib/x86-64_linux/static_pic
ILOLP_INCLUDE= /opt/ibm/ILOG/CPLEX_Studio126/concert/include

#CPLEXLIBDIR =  /media/646bcabe-12bb-400d-87cb-475f895d4667/fabio/ILOG/CPLEX_Studio_AcademicResearch126/cplex/lib/x86-64_linux/static_pic
#LP_INCLUDE= /media/646bcabe-12bb-400d-87cb-475f895d4667/fabio/ILOG/CPLEX_Studio_AcademicResearch126/cplex/include/ilcplex

# Nothing should be changed

#LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread
LP_LIBS = -L$(CONCERTLIBDIR) -L$(CPLEXLIBDIR) -fopenmp -lilocplex -lcplex -lconcert -lm -lpthread 
#LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm -lpthread -liberty -lbfd


#LP_LIBS    = c:/ilog/cplex100/lib/x86_.net2005_8.0/stat_mda/cplex100.lib
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_RHEL3.0_3.2/static_pic/libcplex.a
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_sles9.0_3.3/static_pic/libcplex.a
#LP_LIBS    = /usr/ilog/cplex100/lib/x86-64_rhel4.0_3.4/static_pic/libcplex.a

#LP_INCLUDE= c:/ilog/cplex100/include/ilcplex
#LP_INCLUDE= /usr/ilog/cplex100/include/ilcplex

#DBG= -O3
DBG= -g

#DEFS = $(OS_VERSION) $(COMPILER) $(LP_SOLVERS)

#INCDIR = -I. -I$(LP_INCLUDE) -I$(CONCORDE_INCLUDE)
INCDIR = -I. -I$(ILOLP_INCLUDE) -I$(CPXLP_INCLUDE)

#COMPILER FLAGS

#CXXFLAGS =  $(DBG) $(DEFS) $(INCDIR)
#IBM cluster compiler flags
#CXXFLAGS =  $(DBG) $(INCDIR) -mcpu=powerpc64 -maix64
#edu215 compiler flags
CXXFLAGS = -std=c++0x -O $(DBG) $(INCDIR) -fPIC -DNDEBUG -fexceptions -DIL_STD

main.sh: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main.sh  $(OBJS) $(LP_LIBS)

#$(CXX) -c $(CXXFLAGS) $< -o $@
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

global_variables.o: global_variables.cpp global_variables.h
	$(CXX) $(CXXFLAGS) -c global_variables.cpp

global_functions.o: global_functions.cpp global_functions.h
	$(CXX) $(CXXFLAGS) -c global_functions.cpp
	
BB_FORM_EXT.o: BB_FORM_EXT.cpp BB_FORM_EXT.h
	$(CXX) $(CXXFLAGS) -c BB_FORM_EXT.cpp

graph.o: graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -c graph.cpp

digraph.o: digraph.cpp digraph.h
	$(CXX) $(CXXFLAGS) -c digraph.cpp

solution.o: solution.cpp solution.h
	$(CXX) $(CXXFLAGS) -c solution.cpp

#LDLIBS = $(CONCORDE_LIBS) $(MY_LIBS) $(LP_LIBS)
#LDLIBS = $(CONCORDE_LIBS) $(LP_LIBS)
LDLIBS = $(LP_LIBS)

all: BB

#BB: $(OBJS)
#	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS)

#$(OBJS): Makefile

clean:
	rm -f $(OBJS)
	rm -f main.sh
#	rm BB

