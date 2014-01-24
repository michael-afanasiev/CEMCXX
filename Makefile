CXXFLAGS = -g -std=c++11
CFLAGS   = -g
LDFLAGS  = -L/usr/local/lib -lexoIIv2c -lnetcdf -lhdf5

OBJS= ./src/Exodus_file.o \
      ./src/Mesh.o \
			./src/Read_driver.o \
			./src/Model_file.o \
			./src/Utilities.o \
			./src/Kdtree.o \
			./src/Interpolater.o \
			./src/Construct.o \
			
###################
all: Construct

Construct: $(OBJS)
	$(CXX) $(LDFLAGS) -o Construct $(OBJS)
	
clean:
	$(RM) ./src/*.o ./Construct
