CXX = /usr/local/bin/g++-4.8
CC  = /usr/local/bin/gcc-4.8

CXXFLAGS = -O3 -std=c++11
CFLAGS   = -O3
LDFLAGS  = -L/usr/local/lib -lexoIIv2c -lnetcdf

OBJS= \
			./src/Construct.o \
      ./src/Exodus_file.o \
      ./src/Mesh.o \
			./src/Read_driver.o \
			./src/Model_file.o \
			./src/Utilities.o \
			./src/Kdtree.o
			
###################
all: Construct

Construct: $(OBJS)
	$(CXX) $(LDFLAGS) -o Construct $(OBJS)
	
clean:
	$(RM) ./src/*.o ./Construct
