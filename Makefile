PROG_VERSION := "0.0.1"

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: densityfold

CC          := g++

SRCDIR      := .

LIBS        := -lz -lm 
CFLAGS      := -w -fno-pic
CXXFLAGS    := -w -fno-pic -DPROG_VERSION=\"$(PROG_VERSION)\" $(INCS) -std=c++11
LDFLAGS     := #-static

CLASPOBJ     = $(CLASPDIR)/*.o

BWAOBJ       = $(BWADIR)/*.o

densityfold: 
	$(CC) densityfold.cpp $(CXXFLAGS) -o densityfold
	$(CC) densityfold_lin.cpp $(CXXFLAGS) -o densityfold-lin

clasplib: 
	@$(MAKE) -C $(CLASPDIR)

bwalib:
	@$(MAKE) -C $(BWADIR)

clean:
	@rm -f $(LORDFASTOBJ)
	@rm -f densityfold densityfold-lin

DEBUG_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O0 -ggdb)
	$(eval LIBS = $(LIBS) -O0 -ggdb)

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O3)

PROFILE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)

