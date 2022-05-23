# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

include ../MakeVariables_1
include ../MakeVariables_2
include ../MakeVariables_3

CFLAGS = -fPIC -O3
CXXFLAGS = -Wall -fPIC -O3

LDFLAGS = -fPIC


64bitLibrary = libcatsim.so

$(64bitLibrary): $(OBJECTS)
	@echo ............................................
	@echo .
	@echo . Building system-specific library ...
	@echo .
	@echo .
	$(CXX) $(CXXFLAGS) -static-libgcc -static-libstdc++ -o $(64bitLibrary) -shared $(OBJECTS) $(LDFLAGS)
	mv libcatsim.so ../../catsim/lib
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
