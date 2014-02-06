CXXFLAGS = -O3 -std=c++11
CFLAGS   = -O3
LDFLAGS  = -L/Users/michaelafanasiev/Development/lib -lexoIIv2c -L/usr/local/lib -lnetcdf

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
