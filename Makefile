NETCDF_LIB = /users/afanasm/lib/monch #/apps/monch/netcdf/4.3.1/intel/14.0.1/nonparallel
EXODUS_LIB = /users/afanasm/lib/monch #/apps/monch/exodus/6.02/intel/14.0.1/nonparallel
HDF5_LIB   = /users/afanasm/lib/monch

EXODUS_INC = /users/afanasm/include/monch
NETCDF_INC = /users/afanasm/include/monch

CXX      = icpc
CC       = icc
CXXFLAGS = -O3 -std=c++11 -fopenmp -I$(EXODUS_INC) -I$(NETCDF_INC) 
CFLAGS   = -O3

LDFLAGS = -fopenmp -L$(EXODUS_LIB) -lexodus -L$(NETCDF_LIB) -lnetcdf_c++4 -lnetcdf -L$(HDF5_LIB) -lhdf5_hl -lhdf5 -lz

OBJS_GENERIC     = \
 	./src/Exodus_file.o \
	./src/Mesh.o \
	./src/Driver.o \
	./src/Model_file.o \
	./src/Utilities.o \
	./src/Kdtree.o \
	./src/Interpolater.o \
	./src/Discontinuity.o \
	./src/Background_models.o \
	./src/Attenuation.o									 								 
								 
OBJS_INTERPOLATE = ./src/Interpolate.o

OBJS_EXTRACT_SES = ./src/Extract_ses3d.o

OBJS_EXTRACT_SPC = ./src/Extract_specfem.o

OBJS_CRUST       = ./src/Add_crust.o

OBJS_REFINE      = ./src/Refine.o

OBJS_TOPO        = ./src/Add_topography.o

###################
all: construct extract_spc extract_s3d crust refine topo

construct: $(OBJS_INTERPOLATE) $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o bin/interpolate $(OBJS_INTERPOLATE)  $(OBJS_GENERIC) $(LDFLAGS)
	
extract_s3d:   $(OBJS_EXTRACT_SES) $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o bin/extract_s3d $(OBJS_EXTRACT_SES)  $(OBJS_GENERIC) $(LDFLAGS)
	
extract_spc:   $(OBJS_EXTRACT_SPC) $(OBJS_GENERIC)   
	$(CXX) $(LDFLAGS) -o bin/extract_spec $(OBJS_EXTRACT_SPC) $(OBJS_GENERIC) $(LDFLAGS)
	
crust:     $(OBJS_CRUST)   $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o bin/add_crust   $(OBJS_CRUST)        $(OBJS_GENERIC) $(LDFLAGS)
	
refine:    $(OBJS_REFINE)  $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o bin/refine      $(OBJS_REFINE)       $(OBJS_GENERIC) $(LDFLAGS)
		
topo:      $(OBJS_TOPO)  $(OBJS_GENERIC)
	$(CXX) $(LDFLAGS) -o bin/topo        $(OBJS_TOPO)         $(OBJS_GENERIC) $(LDFLAGS)
	
###################	
clean:
	$(RM) ./src/*.o ./bin/*
