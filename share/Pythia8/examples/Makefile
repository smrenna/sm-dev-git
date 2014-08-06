# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2014 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Include the configuration and handle GZIP support.
-include Makefile.inc
ifeq ($(GZIP_USE),true)
  CXX_COMMON+=-L$(BOOST_LIB) -L$(GZIP_LIB) -lboost_iostreams -lz
endif

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make mainXX"

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be installed, please run "make install" \
                in the top PYTHIA directory)

# Examples without external dependencies.
main% : main%.cc $(PREFIX_LIB)/libpythia8.a $(PREFIX_LIB)/libpythia8lhapdf5.a
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -o $@ $^

# HEPMC2.
main41 main42: $$@.cc $(PREFIX_LIB)/libpythia8.a \
	       $(PREFIX_LIB)/libpythia8lhapdf5.a
ifeq ($(HEPMC2_USE),true)
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -I$(HEPMC2_INCLUDE) -o $@ $^ \
	-L$(HEPMC2_LIB) -Wl,-rpath $(HEPMC2_LIB) -lHepMC
else
	@echo "Error: $@ requires HEPMC2"
endif

# LHAPDF5.
main51 main52 main53 main54: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(LHAPDF5_USE),true)
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -o $@ $^ \
	-L$(LHAPDF5_LIB) -Wl,-rpath $(LHAPDF5_LIB) -lLHAPDF -lgfortran
else
	@echo "Error: $@ requires LHAPDF5"
endif

# HEPMC2 and LHAPDF5.
main61 main62 main85 main86 main87 main88: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(HEPMC2_USE)$(LHAPDF5_USE),truetrue)
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -I$(HEPMC2_INCLUDE) -o $@ $^ \
	-L$(HEPMC2_LIB) -Wl,-rpath $(HEPMC2_LIB) -lHepMC \
	-L$(LHAPDF5_LIB) -Wl,-rpath $(LHAPDF5_LIB) -lLHAPDF -lgfortran
else
	@echo "Error: $@ requires HEPMC2 and LHAPDF5"
endif

# FASTJET3.
main71 main72: $$@.cc $(PREFIX_LIB)/libpythia8.a \
               $(PREFIX_LIB)/libpythia8lhapdf5.a
ifeq ($(FASTJET3_USE),true)
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -I$(FASTJET3_INCLUDE) -o $@ $^\
	-L$(FASTJET3_LIB) -Wl,-rpath $(FASTJET3_LIB) -lfastjet
else
	@echo "Error: $@ requires FASTJET3"
endif

# FASTJET3, HEPMC2, and LHAPDF5.
main81 main82 main83 main84: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE)$(HEPMC2_USE)$(LHAPDF5_USE),truetruetrue)
	$(CXX) $(CXX_COMMON) -I$(PREFIX_INCLUDE) -I$(FASTJET3_INCLUDE) \
	-I$(HEPMC2_INCLUDE) -o $@ $^ -\
	-L$(HEPMC2_LIB) -Wl,-rpath $(HEPMC2_LIB) -lHepMC \
	-L$(LHAPDF5_LIB) -Wl,-rpath $(LHAPDF5_LIB) -lLHAPDF -lgfortran \
	-L$(FASTJET3_LIB) -Wl,-rpath $(FASTJET3_LIB) -lfastjet
else
	@echo "Error: $@ requires FASTJET3, HEPMC2, and LHAPDF5"
endif

# ROOT (turn off all warnings for readability).
main91: $$@.cc $(PREFIX_LIB)/libpythia8.a $(PREFIX_LIB)/libpythia8lhapdf5.a
ifeq ($(ROOT_USE),true)
	$(CXX) $(CXX_COMMON) -w -I$(PREFIX_INCLUDE) -I$(ROOT_INCLUDE) -o $@ $^ \
	-Wl,-rpath $(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif
main92: $$@.cc $$@.h $$@LinkDef.h $(PREFIX_LIB)/libpythia8.a \
        $(PREFIX_LIB)/libpythia8lhapdf5.a
ifeq ($(ROOT_USE),true)
	export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB); \
	$(ROOT_BIN)rootcint -f $@Dct.cc -c -I$(PREFIX_INCLUDE) $@.h $@LinkDef.h
	$(CXX) $(CXX_COMMON) -w -I$(PREFIX_INCLUDE) -I$(ROOT_INCLUDE) \
	-o $@ $@Dct.cc $^ \
	-Wl,-rpath $(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*