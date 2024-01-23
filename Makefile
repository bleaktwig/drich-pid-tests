# C++ compiler + version.
CXX         := g++
CXX_STD     := -std=c++20

# C++ flags.
# NOTE. Flags not included in this list:
#   * -pedantic and -Wextra: I don't agree with the C++ specification that va-
#         riable-length arrays are a bad thing, and I don't feel like adding the
#         `__extension__` statement every time I want to use one. In any case, a
#         variable-length array is a much better solution memory-wise than an
#         std::vector.
#   * -Wold-style-cast: Similar logic to -pedantic, with the added reason that
#         this flag breaks dd4hep.
#   * -Wsign-conversion, -Wcast-qual: These flags break dd4hep.
# In addition, I had to add the flag -Wno-sign-compare because dd4hep compares
#     int32_t with size_t type.
CFLAGS_DBG  := -g -Wall -Wcast-align -Wctor-dtor-privacy \
			   -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op \
			   -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept \
			   -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo \
			   -Wstrict-null-sentinel -Wstrict-overflow=4 -Wswitch-default \
			   -Wundef -Werror -Wno-unused -Wno-sign-compare
CFLAGS_PROD := -O3
CXX         := $(CXX) $(CFLAGS_DBG)

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
	@rm bin/root2csv
