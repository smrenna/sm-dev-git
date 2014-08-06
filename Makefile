# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2014 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
#
# This is is the Makefile used to build PYTHIA on POSIX systems. Example usage 
# is:
#     make -j2
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script
# and the distribution structure.
################################################################################

# Include the configuration and set the local directory structure.
ifeq (,$(findstring clean, $(MAKECMDGOALS)))
  -include Makefile.inc
endif
LOCAL_BIN=bin
LOCAL_EXAMPLES=examples
LOCAL_HTMLDOC=htmldoc
LOCAL_INCLUDE=include
LOCAL_LIB=lib
LOCAL_PDFDOC=pdfdoc
LOCAL_PHPDOC=phpdoc
LOCAL_SRC=src
LOCAL_TMP=tmp
LOCAL_XMLDOC=xmldoc
LOCAL_MKDIRS:=$(shell mkdir -p $(LOCAL_TMP) $(LOCAL_LIB))

# PYTHIA.
OBJECTS=$(patsubst $(LOCAL_SRC)/%.cc,$(LOCAL_TMP)/%.o, \
	$(wildcard $(LOCAL_SRC)/*.cc))
TARGETS=$(LOCAL_LIB)/libpythia8.a
ifeq ($(ENABLE_SHARED),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX)
endif

# LHAPDF5.
TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf5.a
ifeq ($(ENABLE_SHARED),true)
  TARGETS+=$(LOCAL_LIB)/libpythia8lhapdf5$(LIB_SUFFIX)
endif

################################################################################
# RULES: Definition of the rules used to build PYTHIA.
################################################################################

# Rules without physical targets.
.PHONY: all install clean fullclean

# All targets.
all: Makefile.inc $(TARGETS)

# The Makefile configuration.
Makefile.inc:
	./configure

# Auto-generated (with -MD flag) dependencies.
-include $(LOCAL_TMP)/*.d

# PYTHIA.
$(LOCAL_TMP)/%.o: $(LOCAL_SRC)/%.cc
ifeq ($(GZIP_USE),true)
	$(CXX) -MD $(CXX_COMMON) -I$(LOCAL_INCLUDE) -I$(BOOST_INCLUDE) \
	-DGZIPSUPPORT -c -o $@ $^
else
	$(CXX) -MD $(CXX_COMMON) -I$(LOCAL_INCLUDE) -c -o $@ $^
endif
$(LOCAL_LIB)/libpythia8.a: $(OBJECTS)
	ar cru $@ $^
$(LOCAL_LIB)/libpythia8$(LIB_SUFFIX): $(OBJECTS)
	$(CXX) $(CXX_COMMON) $(CXX_SHARED) -o $@ $^ $(CXX_SONAME),$(notdir $@)

# LHAPDF5.
$(LOCAL_TMP)/LHAPDF5.o: $(LOCAL_INCLUDE)/Pythia8Plugins/LHAPDF5.h
	$(CXX) -MD $(CXX_COMMON) -I$(LOCAL_INCLUDE) -x c++ -c -o $@ $^
$(LOCAL_LIB)/libpythia8lhapdf5.a: $(LOCAL_TMP)/LHAPDF5.o
	ar cru $@ $^
$(LOCAL_LIB)/libpythia8lhapdf5$(LIB_SUFFIX): $(LOCAL_TMP)/LHAPDF5.o
	$(CXX) $(CXX_COMMON) $(CXX_SHARED) -o $@ $^ $(CXX_SONAME),$(notdir $@)

# Install (rsync is used for finer control).
install: all
	@mkdir -p $(PREFIX_BIN) $(PREFIX_INCLUDE) $(PREFIX_LIB) $(PREFIX_SHARE)
	@mkdir -p $(PREFIX_SHARE)/doc
	@rsync -a pythia8-config.in $(LOCAL_BIN)
	@rsync -a $(LOCAL_BIN)/* $(PREFIX_BIN) --exclude .svn
	@rsync -a $(LOCAL_INCLUDE)/* $(PREFIX_INCLUDE) --exclude .svn
	@rsync -a $(LOCAL_LIB)/* $(PREFIX_LIB) --exclude .svn
	@rsync -a  Makefile.inc $(LOCAL_EXAMPLES)
	@rsync -a $(LOCAL_EXAMPLES) $(PREFIX_SHARE) --exclude .svn
	@rsync -a README $(PREFIX_SHARE)/doc
	@rsync -a $(LOCAL_HTMLDOC) $(PREFIX_SHARE)/doc --exclude .svn
	@rsync -a $(LOCAL_PDFDOC) $(PREFIX_SHARE)/doc --exclude .svn
	@rsync -a $(LOCAL_PHPDOC) $(PREFIX_SHARE)/doc --exclude .svn
	@rsync -a $(LOCAL_XMLDOC) $(PREFIX_SHARE)/doc --exclude .svn

# Clean.
clean:
	@rm -rf $(LOCAL_TMP) $(LOCAL_LIB)
	@rm -f $(LOCAL_EXAMPLES)/*Dct.*
	@rm -f $(LOCAL_EXAMPLES)/main[0-9][0-9]

# Clean all temporary and generated files.
fullclean: clean
	@find . -type f -name Makefile.inc -print0 | xargs -0 rm -f
	@find . -type f -name "*~" -print0 | xargs -0 rm -f
	@find . -type f -name "#*" -print0 | xargs -0 rm -f