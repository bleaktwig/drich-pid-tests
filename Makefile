# C++ compiler + version.
CXX     := g++
CXX_STD := -std=c++20

CFLAGS_DBG  := -g -Wall -Wcast-align -Wcast-qual \
			   -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 \
			   -Winit-self -Wlogical-op -Wmissing-declarations \
			   -Wmissing-include-dirs -Wnoexcept -Wold-style-cast \
			   -Woverloaded-virtual -Wredundant-decls -Wshadow \
			   -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel \
			   -Wstrict-overflow=4 -Wswitch-default -Wundef -Werror -Wno-unused
CFLAGS_PROD := -O3
CXX         := $(CXX) $(CFLAGS_DBG)

# ROOT.
ROOTCFLAGS  := -pthread $(CXX_STD) -m64 -isystem$(shell root-config --incdir)
RLIBS       := $(shell root-config --libs) -lEG
RXX         := $(CXX) $(ROOTCFLAGS)

# Targets.
all: root2csv

root2csv: root2csv.c
	$(RXX) -o bin/root2csv root2csv.c $(RLIBS)

clean:
	@echo "Removing all build files and binaries."
	@rm root2csv
