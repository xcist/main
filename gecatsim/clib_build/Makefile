all: libcatsim.so

include ../MakeVariables_1
include ../MakeVariables_2
include ../MakeVariables_3

CFLAGS = -fPIC -O3
CXXFLAGS = -Wall -fPIC -O3

LDFLAGS = -fPIC -lpthread

64bitLibrary = libcatsim.so

libcatsim.so: $(OBJECTS)
	@echo ............................................
	@echo .
	@echo . Building system-specific library ...
	@echo .
	@echo .
	$(CXX) $(CXXFLAGS) -static-libgcc -o $(64bitLibrary) -shared $(OBJECTS) $(LDFLAGS)
	mv libcatsim.so ../../lib/
	@echo .
	@echo .
	@echo . Done building system-specific library.
	@echo .
	@echo ............................................

clean:
	@echo ............................................
	@echo .
	@echo . Cleaning ...
	@echo .
	@echo .
	rm -f *.o
	@echo .
	@echo .
	@echo . Done cleaning.
	@echo .
	@echo ............................................

# DO NOT DELETE
