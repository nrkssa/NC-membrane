EXENAME :=free_ener
	
CXX :=mpic++
	
FC :=mpif90
	
CXXFLAGS :=  -Wall -O2 -Wuninitialized -std=c++11

FCFLAGS :=  -O2 -Wall  -Wno-tabs -fimplicit-none  -ffree-line-length-0  --free-form -fno-second-underscore  -fno-underscoring -fno-range-check -Wuninitialized

LDCXX :=

CXXF2C :=-lgfortran 
	
CC_LOCALINCLUDEPATH :=-I./CC_INCLUDE
FC_LOCALINCLUDEPATH :=-I./FC_INCLUDE  
vpath %.cpp ./CC_SOURCE
vpath %.f ./FC_SOURCE

CPP_SOURCE_DIR=./CC_SOURCE
FC_SOURCE_DIR=./FC_SOURCE

_CPP_SOURCE := debugger.cpp parser.cpp  datainput.cpp interaction_parameters.cpp nanocarrier.cpp receptor.cpp membrane.cpp setup.cpp potential.cpp selfavoidance.cpp init.cpp linklist.cpp distribution.cpp movement.cpp histogram.cpp datawriter.cpp restart.cpp fortran_common_functions.cpp  Maincode-freeenergy.cpp
CPP_SOURCE = $(patsubst %,$(CPP_SOURCE_DIR)/%,$(_CPP_SOURCE))

_FC_SOURCE := module_string_utility.f module_randomnumber.f module_datastruct.f module_initialize_curvature_field.f module_curvcalc.f module_writedata.f module_makesurface.f module_mcsmoves.f
FC_SOURCE = $(patsubst %,$(FC_SOURCE_DIR)/%,$(_FC_SOURCE))

COBJS := $(CPP_SOURCE:.cpp=.o)
FOBJS :=$(FC_SOURCE:.f=.o)
FMOD :=$(FC_SOURCE:.f=.mod)
FMODS :=$(FC_SOURCE:.f=.mod)
OBJS :=$(FOBJS) $(COBJS)

all : $(EXENAME)

$(EXENAME):$(OBJS)
	@echo " "
	@echo " "
	@echo "Building executable with GNU libraries"
	$(CXX) -o $(EXENAME)  $(CXXFLAGS) $(CC_LOCALINCLUDEPATH) $(OBJS) $(CXXF2C) $(LDCXX)
	@echo " "

$(COBJS): %.o:%.cpp
	$(CXX) $(CXXFLAGS)  $(CC_LOCALINCLUDEPATH) -c $< -o $*.o
	@echo "-------------------------------------------->"
	

$(FOBJS):%.o:%.f
	$(FC) $(FCFLAGS) $(FC_LOCALINCLUDEPATH) -c $< -o $*.o
	@mv *.mod $(FC_SOURCE_DIR)
	@echo "-------------------------------------------->"

.PHONY: clean
clean:
	@rm -fv *.vtu *.dat *.pdb *.tcl *.mod $(FMOD) $(OBJS) $(EXENAME) $(FMODS) ./SYS_STATE/MEMB/*   ./SYS_STATE/NC-DATA/*.vtu ./SYS_STATE/ANTIGEN/* ./SYS_STATE/SYSTEM/*.vtu

.PHONY: cleanall
cleanall:mvdump
	@rm -fv *.vtu *.dat *.pdb *.tcl *.mod $(FMOD) $(OBJS) $(EXENAME) $(FMODS) ./SYS_STATE/MEMB/*   ./SYS_STATE/NC-DATA/*.vtu ./SYS_STATE/ANTIGEN/* ./SYS_STATE/SYSTEM/*.vtu


.PHONY: cleandata
cleandata:
	@echo $(COBJS)
	@echo $(FOBJS)
	@rm -fv *.vtu *.dat *.pdb *.tcl ./SYS_STATE/MEMB/* ./SYS_STATE/NC-DATA/*.vtu ./SYS_STATE/ANTIGEN/* ./SYS_STATE/SYSTEM/*.vtu ./SYS_STATE/DUMP/*.dat ./SYS_STATE/DUMP/*.dump

.PHONY: mvdump
mvdump:
	@mv -v ./SYS_STATE/RESTART ./SYS_STATE/RESTART-`date +%d%m%y-%k%M`
	@mv -v ./SYS_STATE/DUMP   ./SYS_STATE/RESTART
	@mkdir -v ./SYS_STATE/DUMP

