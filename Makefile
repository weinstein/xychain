# We'll expect that ./bin, ./build, and ./doc exist later
GEN_FOLDERS = ./bin ./build ./doc


#INCLUDE_PATH_FLAGS := -I"$(shell pwd)/include"

#LINK_PATH_FLAGS := -L/usr/local/lib/

# === Are we on Linux (Jack) or Mac-OSX (Kay)? ===

UNAME = $(shell uname -s)

# === "I'll just have one of everything on the menu" ===

ifdef AWESOME
VISUALS = "yes"
PNG_VIS = "oh, yes"
override undefine DEBUGGING
endif

# === External Library Flags ===

ifdef VISUALS
SDL_LINK_FLAGS = $(shell sdl-config --libs)
SDL_COMPILER_FLAGS = $(shell sdl-config --cflags)
  ifdef PNG_VIS
  PNG_LINKER_FLAGS = $(shell libpng-config --ldflags)
  PNG_COMPILER_FLAGS = $(shell libpng-config --cflags)
  endif
endif



GSL_LINK_FLAGS = $(shell gsl-config --libs)
GSL_COMPILER_FLAGS = $(shell gsl-config --cflags)

# Possibly ensure that the GFLAGS_DIR is the location where gflags was installed
GFLAGS_DIR = /usr/local/lib
GFLAGS = -L$(GFLAGS_DIR) -lgflags -Wl,-rpath -Wl,$(GFLAGS_DIR)

MPI_COMPILER_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)

ifeq ($(UNAME),Linux)
GL_FLAGS = -lGL -lGLU
endif
ifeq ($(UNAME),Darwin)
GL_FLAGS =  -framework OpenGL
endif

# === Compilation Flags ===

COMPILER = g++
COMPILER_FLAGS = -c -std=c++0x -pthread $(PNG_COMPILER_FLAGS) $(SDL_COMPILER_FLAGS) $(GL_FLAGS) $(GSL_COMPILER_FLAGS) $(MPI_COMPILER_FLAGS)#-ldl

ifeq ($(UNAME),Linux)
COMPILER_FLAGS += -DLINUX
endif
ifeq ($(UNAME),Darwin)
COMPILER_FLAGS += -DOSX
endif

NODEBUG = -O3
DEBUG = -O0 -g
ifdef DEBUGGING
COMPILER_FLAGS += $(DEBUG)
endif
ifndef DEBUGGING
COMPILER_FLAGS += $(NODEBUG)
endif

ifdef VISUALS
COMPILER_FLAGS += -DVISUALS
  ifdef PNG_VIS
  COMPILER_FLAGS += -DPNG_VIS
  endif
endif

LINKER = g++
LINKER_FLAGS = $(GFLAGS) -pthread $(PNG_LINKER_FLAGS) $(SDL_LINK_FLAGS) $(GL_FLAGS) $(GSL_LINK_FLAGS) $(MPI_LINK_FLAGS) $(LINK_PATH_FLAGS)

# === Common Files ===

COMMON = build/threadpool.o $(GEN_FOLDERS)
COMMON_OBJECTS = build/threadpool.o

ifdef VISUALS
COMMON += build/viewer.o
COMMON_OBJECTS += build/viewer.o
endif

# ===================================================================

.PHONY : clean tidy backup test test-valgrind

all : xychain backup

build/%.o : source/%.cpp $(GEN_FOLDERS)
	$(COMPILER) $(COMPILER_FLAGS) -c $< -o $@

# ====== Testing =======

TESTS = threadpool-test spins-test random-test xychain-test

test-valgrind : $(TESTS)
	$(foreach var,$(TESTS),valgrind --leak-check=full ./bin/$(var);)

test : $(TESTS)
	$(foreach var,$(TESTS),./bin/$(var);)

# ====== Backup  =======

ALL_SOURCES = $(shell echo ./source/*) $(shell echo ./scripts/*) ./Makefile
TIMESTAMP = $(shell date +%s)
backup : ./backups ./backups/latest
backups/latest : $(ALL_SOURCES)
	zip -r ./backups/$(TIMESTAMP).zip $(ALL_SOURCES)
	ln -s -f $(TIMESTAMP).zip ./backups/latest

# ====== Folders =======

./bin :
	mkdir -p bin

./build : 
	mkdir -p build

./doc :
	mkdir -p doc

./backups :
	mkdir -p backups

# ====== Executables =======

threadpool-test : build/threadpool.o build/threadpool-test.o $(GEN_FOLDERS)
	$(LINKER) build/threadpool.o build/threadpool-test.o $(LINKER_FLAGS) -o bin/$@

spins-test: build/spins-test.o $(GEN_FOLDERS)
	$(LINKER) build/spins-test.o $(LINKER_FLAGS) -o bin/$@

random-test: build/random-test.o $(GEN_FOLDERS)
	$(LINKER) build/random-test.o $(LINKER_FLAGS) -o bin/$@

xychain-test: build/xychain.o build/xychain-test.o $(COMMON)
	$(LINKER) build/xychain.o build/xychain-test.o $(COMMON_OBJECTS) $(LINKER_FLAGS) -o bin/$@

xychain-mpi-test-main: build/xychain-mpi.o build/xychain-mpi-test-main.o $(COMMON)
	$(LINKER) build/xychain-mpi.o build/xychain-mpi-test-main.o $(COMMON_OBJECTS) $(LINKER_FLAGS) -o bin/$@

xychain : build/xychain.o build/xychain-cli.o build/xychain-cli-main.o $(COMMON)
	$(LINKER) build/xychain.o build/xychain-cli.o build/xychain-cli-main.o $(COMMON_OBJECTS) $(LINKER_FLAGS) -o bin/$@

xychain-mpi : build/xychain-mpi.o build/xychain-cli-mpi.o build/xychain-cli-mpi-main.o $(COMMON)
	$(LINKER) build/xychain-mpi.o build/xychain-cli-mpi.o build/xychain-cli-mpi-main.o $(COMMON_OBJECTS) $(LINKER_FLAGS) -o bin/$@

avg-runs: scripts/AverageRuns.cpp
	g++ --std=c++0x scripts/AverageRuns.cpp -o bin/avg-runs

bin2ascii: scripts/bin2ascii.cpp
	g++ --std=c++0x scripts/bin2ascii.cpp -o bin/bin2ascii

# ======= Clean up ============

clean :
	-rm -f -r ./*.o ./build/*.o ./bin/* ./doc/*

tidy :
	-rm -f *.png *.root *.dat
