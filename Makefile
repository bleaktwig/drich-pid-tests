# C++ compiler + version.
CXX         := g++
CXX_STD     := -std=c++20

# C++ flags.
CFLAGS_DBG  := -g -Wall -Wcast-align -Wcast-qual \
			   -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 \
			   -Winit-self -Wlogical-op -Wmissing-declarations \
			   -Wmissing-include-dirs -Wnoexcept -Wold-style-cast \
			   -Woverloaded-virtual -Wredundant-decls -Wshadow \
			   -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel \
			   -Wstrict-overflow=4 -Wswitch-default -Wundef -Werror -Wno-unused
CFLAGS_PROD := -O3
CXX         := $(CXX) $(CFLAGS_PROD)

# ROOT.
ROOTCFLAGS  := -pthread $(CXX_STD) -m64 -isystem$(shell root-config --incdir)
RXX         := $(CXX) $(ROOTCFLAGS)
RLIBS       := $(shell root-config --libs) -lEG

# DD4HEP.
# DD4HEPCFLAGS := -I/path/to/dd4hep/include
DD4HEPCFLAGS := -I/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/dd4hep-1.27-luzqi3jwjmeycjkbed3j6pbuq4en45j7/include
DXX          := $(RXX) $(DD4HEPCFLAGS)
DD4HEPLIBS   := -L/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/dd4hep-1.27-luzqi3jwjmeycjkbed3j6pbuq4en45j7/lib -lDDCore

# Targets.
all: root2csv

root2csv: root2csv.c
	$(DXX) -o bin/root2csv root2csv.c $(RLIBS) $(DD4HEPLIBS)

clean:
	@rm root2csv
