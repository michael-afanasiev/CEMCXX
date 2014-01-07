CPPFLAGS=-g -std=c++11
LDFLAGS=-L/usr/local/lib -lexoIIv2c -lnetcdf

OBJS= \
      ./src/Exodus_file.o \
      ./src/Mesh.o \
			./src/read_driver.o \
      ./src/construct.o \
			
###################
%.o: ./src/%.cpp
	$(CXX) $(CPPFLAGS) -c $<

all: Construct

Construct: $(OBJS)
	$(CXX) $(LDFLAGS) -o Construct $(OBJS)
	
clean:
	$(RM) ./src/*.o
