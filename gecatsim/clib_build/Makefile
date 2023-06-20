all: libcatsim.so

include ../MakeVariables_1
include ../MakeVariables_2
include ../MakeVariables_3

CFLAGS = -fPIC -O3
CXXFLAGS = -Wall -fPIC -O3

LDFLAGS = -fPIC -lpthread

p_nlog_inline.o: p_nlog_inline.c
	gcc -O0 -g -Wall -msse -fPIC -c p_nlog_inline.c

64bitLibrary = libcatsim.so

libcatsim.so: $(OBJECTS)
	@echo ............................................
	@echo .
	@echo . Building system-specific library ...
	@echo .
	@echo .
	$(CXX) $(CXXFLAGS) -static-libgcc -static-libstdc++ -o $(64bitLibrary) -shared $(OBJECTS) $(LDFLAGS)
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
