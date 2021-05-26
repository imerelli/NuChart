CC					= gcc
CXX					= g++
AS					= nasm
LINK_OPT			=
#OPTIMIZE_FLAGS		= -O3 -finline-functions #-DTRACE_FASTFLOW
CXXFLAGS			= -Wall -std=c++0x -DWEIGHIT #-DPRINTLOG -DVERBOSE -DTEST3
CFLAGS				=
LDFLAGS				=
INCS				= -I /home/dagostino/gsl/include #-I /include #-I #glm/ #-I../../../fastflow-dev/
LIBS				= -L /home/dagostino/gsl/lib -lcurl -lm -lgsl -lgslcblas -lpthread
#ALLOC_STATS          = -DALLOCATOR_STATS

MACHINE := $(shell uname -n)

# compile for debug or release
ifeq ($(DEBUG), 1)
OPTIMIZE_FLAGS		= -O0 -g #-DTRACE_FASTFLOW
else
OPTIMIZE_FLAGS		= -O3 -finline-functions #-DTRACE_FASTFLOW
endif

ifeq ($(EXPINTER), 1)
CXXFLAGS			+= -DEXPINTER
endif

ifeq ($(NUMA), 1)
CXXFLAGS			+= -DUSE_NUMA
LIBS				+= -lnuma
endif

ifeq ($(NORM), 1)
CXXFLAGS			+= -DWEIGHIT
endif

ifeq ($(TEST_EDGES), 1)
CXXFLAGS			+= -DTEST_EDGES
endif

ifeq ($(TEST2), 1)
CXXFLAGS			+= -DTEST2
endif

ifeq ($(TEST3), 1)
CXXFLAGS			+= -DTEST3
endif

ifeq ($(TEST4), 1)
CXXFLAGS			+= -DTEST4
endif

ifeq ($(OPENMP), 1)
OPTIMIZE_FLAGS		+= -fopenmp
CXXFLAGS			+= -DUSE_OPENMP
endif

ifeq ($(TBB), 1)
LIBS				+= -ltbb
CXXFLAGS			+= -DUSE_TBB
endif

#ifneq (,$(findstring paracool.di.unito.it, $(MACHINE)))
#FF_PATH=/home/tordini/FF/mc-fastflow
#else ifneq (,$(findstring fabioXPS, $(MACHINE)))
#FF_PATH=/space/repos/Fastflow/mc-fastflow
#else ifneq (,$(findstring cms47xl.to.infn.it, $(MACHINE)))
#FF_PATH=/home/tordini/fastflow
#else ifneq (,$(findstring paranoia, $(MACHINE)))
#FF_PATH=/home/tordini/mc-fastflow
#else ifneq (,$(findstring titanic, $(MACHINE)))
#FF_PATH=/usr/local/fastflow-2.0.4
#endif
FF_PATH=../fastflow

.SUFFIXES: .c .cpp .o
.PHONY: all clean

INCLUDES			= -I. $(INCS) -I $(FF_PATH)
COMMON				= common.hpp Vertex.hpp Edg.hpp Timings.hpp FileIO.hpp
_DEPS				= Fragment.hpp Gene.hpp SamData.hpp Finder.hpp Edge.hpp CMap.hpp Parsers.hpp PrintOut.hpp glm/Fit.hpp glm/GLM.hpp \
					SamDataT.hpp FragmentT.hpp GeneT.hpp EdgeT.hpp
DEPS				= $(COMMON) $(_DEPS)
OBJ				= Finder.o MemoryMapped.o testth.o
TARGET				= nuchart

all: $(TARGET)

%.o: %.cpp	$(DEPS)
	$(CXX) $(INCLUDES) -c $(OPTIMIZE_FLAGS) $< -o $@ $(CXXFLAGS)

nuchart:	$(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(OPTIMIZE_FLAGS) $(LIBS)

clean:
	rm -fr *.o *~ $(TARGET) *.d
	cd glm && $(MAKE) clean

#cleanall: clean
#	-rm -fr $(TARGET) *.d
