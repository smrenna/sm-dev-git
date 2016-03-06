#!/usr/bin/env bash
# Copyright (C) 2016 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, March 2016.

# This is a script to create a Python interface for Pythia using SWIG.
# It should be run from the main directory, but is located in private/.

# Check SWIG and SED exist.
if ! type "sed" &> /dev/null; then 
    echo "Error: SED not found."; exit; fi
if ! type "swig" &> /dev/null; then 
    echo "Error: SWIG not found."; exit; fi

# Set the SWIG configuration file and the Python interface file.
CFG_FILE="python.cfg"
CXX_FILE=`basename $CFG_FILE .cfg`"_wrap.cxx"
PYTHON_FILE="pythia8.py"
HEADER_FILE="include/Pythia8Plugins/PythonWrapper.h"

# Copy all headers for use with SWIG.
rm -rf python $HEADER_FILE
mkdir -p python/include
DIRS="include/Pythia8 include/Pythia8Plugins"
for DIR in $DIRS; do cp -r $DIR python/$DIR; done

# Handle in-place SED for OSX.
ARCH=$(uname | grep -i -o -e Linux -e Darwin)
ARCH=$(echo $ARCH | awk '{print toupper($0)}')
if [ "$ARCH" == "DARWIN" ]; then SED="sed -r -i .sed"; else SED="sed -r -i"; fi

# Force "std" prefix for all STL classes used.
FILES=`ls python/include/*/*.h`
CLASSES="string vector map deque set ostream ostringstream complex"
for CLASS in $CLASSES; do
    REPLACE="std::$CLASS"
    if [ "$CLASS" == "complex" ]; then REPLACE+="<double> "; fi
    # Middle of line.
    $SED "s|([^[:alnum:]]+)$CLASS([^[:alnum:]]+)|\1$REPLACE\2|g" $FILES
    # Start of line.
    $SED "s|^$CLASS([^[:alnum:]]+)|$REPLACE\1|g" $FILES
    # End if line.
    $SED "s|([^[:alnum:]]+)$CLASS$|\1$REPLACE|g" $FILES
    # Complete line.
    $SED "s|^$CLASS$|$REPLACE|g" $FILES; done
cp include/Pythia8/PythiaComplex.h python/include/Pythia8/PythiaComplex.h

# Write the SWIG configuration header.
cat > $CFG_FILE << BLOCKTEXT
// Set the module name.
%module pythia8

// Include the STL type maps.
%{
#include <cstddef>
%}
%include <typemaps.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_map.i>
%include <std_deque.i>
%include <std_set.i>
%include <std_ios.i>
%include <std_iostream.i>
%include <std_sstream.i>
%include <std_streambuf.i>
%include <complex.i>

// Operators that can be ignored that are handled later.
%ignore Pythia8::Vec4::operator=;
%ignore Pythia8::Vec4::operator[];
%ignore Pythia8::RotBstMatrix::operator=;
%ignore Pythia8::Hist::operator=;
%ignore Pythia8::HEPRUP::operator=;
%ignore Pythia8::HEPEUP::operator=;
%ignore Pythia8::LHmatrixBlock::operator=;
%ignore Pythia8::LHtensor3Block::operator=;
%ignore Pythia8::DecayChannel::operator=;
%ignore Pythia8::ParticleDataEntry::operator=;
%ignore Pythia8::ParticleData::operator=;
%ignore Pythia8::Particle::operator=;
%ignore Pythia8::Junction::operator=;
%ignore Pythia8::Event::operator=;
%ignore Pythia8::Event::operator[];
%ignore Pythia8::FlavContainer::operator=;
%ignore Pythia8::ColConfig::operator[];
%ignore Pythia8::SingleClusterJet::operator=;
%ignore Pythia8::SingleSlowJet::operator=;
%ignore Pythia8::BeamParticle::operator[];

// Methods that should be renamed.
%rename(list) print;

// Members that must be ignored because that cannot be handled by SWIG.
%ignore Pythia8::Writer::headerStream;
%ignore Pythia8::Writer::initStream;
%ignore Pythia8::Writer::eventStream;
%ignore Pythia8::CoupSUSY::LsddX;
%ignore Pythia8::CoupSUSY::LsuuX;
%ignore Pythia8::CoupSUSY::LsduX;
%ignore Pythia8::CoupSUSY::LsudX;
%ignore Pythia8::CoupSUSY::LsvvX;
%ignore Pythia8::CoupSUSY::LsllX;
%ignore Pythia8::CoupSUSY::LsvlX;
%ignore Pythia8::CoupSUSY::LslvX;
%ignore Pythia8::CoupSUSY::RsddX;
%ignore Pythia8::CoupSUSY::RsuuX;
%ignore Pythia8::CoupSUSY::RsduX;
%ignore Pythia8::CoupSUSY::RsudX;
%ignore Pythia8::CoupSUSY::RsvvX;
%ignore Pythia8::CoupSUSY::RsllX;
%ignore Pythia8::CoupSUSY::RsvlX;
%ignore Pythia8::CoupSUSY::RslvX;
%ignore Pythia8::CoupSUSY::rvLLE;
%ignore Pythia8::CoupSUSY::rvLQD;
%ignore Pythia8::CoupSUSY::rvUDD;

// Error flags to allow iterable objects.
%{
#include <assert.h>
static int VEC4_ERROR(0), EVENT_ERROR(0), BEAMPARTICLE_ERROR(0),
  COLCONFIG_ERROR(0);
%}

// Headers and namespace to copy verbatim to interface.
%{
#include "Pythia8/Pythia.h"
using namespace Pythia8;
BLOCKTEXT

# Determine the files.
NAMES=""
BAN="FJCORE FASTJET3 LHAPDF5 LHAPDF6 EXECINFO EVTGEN HEPMC2 LHAPOWHEG"
BAN+=" LHAFORTRAN"
for FILE in $FILES; do
    NAME=`basename $FILE .h | awk '{print toupper($0)}'`
    if [[ $BAN =~ (^| )$NAME($| ) ]]; then continue; fi
    DEPS=`g++ -MM -MG $FILE | grep -o " Pythia8\(Plugins\)\?/.*\.h"`
    eval ${NAME}_DEPS=""; eval ${NAME}_FILE=$FILE
    for DEP in $DEPS; do
	DEP=`basename $DEP .h | awk '{print toupper($0)}'`
	eval ${NAME}_DEPS+=\" ${DEP}\"; done
    DIR="Pythia8"; if [[ $FILE =~ "Pythia8Plugins" ]]; then DIR+="Plugins"; fi
    FILE=`basename $FILE`
    echo "#include \"$DIR/$FILE\"" >> $CFG_FILE
    NAMES+=" $NAME"; done

# Determine the dependencies.
cat >> $CFG_FILE << BLOCKTEXT
%}

// Headers to generate interface for (order is important).
BLOCKTEXT
BAN="PYTHIASTDLIB"
ITER=0; ITERS=`echo $NAMES | wc -w`
while [ $ITER -lt $ITERS ]; do
    for NAME in $NAMES; do
	eval OLD_DEPS=\${${NAME}_DEPS}
	if [ "$OLD_DEPS" == "NONE" ]; then continue; fi
	NEW_DEPS=""
	for DEP in $OLD_DEPS; do
	    eval SUB_DEPS=\${${DEP}_DEPS}
	    if [ "$SUB_DEPS" == "NONE" ]; then continue; fi
	    SIZE=`echo $SUB_DEPS | wc -w`
	    if [ "$SIZE" == "0" ]; then continue
	    else NEW_DEPS+=" $DEP"; fi; done
	if [ `echo $NEW_DEPS | wc -w` -eq "0" ]; then 
	    eval FILE=\${${NAME}_FILE}
	    eval ${NAME}_DEPS=\"NONE\"
	    if [[ $BAN =~ (^| )$NAME($| ) ]]; then
		echo "//%include \"$FILE\"" >> $CFG_FILE
	    else echo "%include \"$FILE\"" >> $CFG_FILE; fi
	else eval ${NAME}_DEPS=\"$NEW_DEPS\"; fi; done
    ITER=$[$ITER+1]; done

# Write the SWIG configuration footer.
cat >> $CFG_FILE << BLOCKTEXT

// Exceptions to handle iterable classes.
%exception Pythia8::Vec4::__getitem__ {
  assert(!VEC4_ERROR); \$action
  if (VEC4_ERROR) {VEC4_ERROR = 0; 
    SWIG_exception(SWIG_IndexError, "Index out of bounds");}
}
%exception Pythia8::Event::__getitem__ {
  assert(!EVENT_ERROR); \$action
  if (EVENT_ERROR) {EVENT_ERROR = 0; 
    SWIG_exception(SWIG_IndexError, "Index out of bounds");}
}
%exception Pythia8::BeamParticle::__getitem__ {
  assert(!BEAMPARTICLE_ERROR); \$action
  if (BEAMPARTICLE_ERROR) {BEAMPARTICLE_ERROR = 0; 
    SWIG_exception(SWIG_IndexError, "Index out of bounds");}
}
%exception Pythia8::ColConfig::__getitem__ {
  assert(!COLCONFIG_ERROR); \$action
  if (COLCONFIG_ERROR) {COLCONFIG_ERROR = 0;
    SWIG_exception(SWIG_IndexError, "Index out of bounds");}
}

// Macro definitions for common class extensions.
%define __STR_OSS__()
  std::string __str__() {
    std::ostringstream oss(std::ostringstream::out); oss << *(\$self);
    return oss.str();}
%enddef
%define __STR_LIST__()
  std::string __str__() {
    std::streambuf* old = cout.rdbuf();
    std::ostringstream oss(std::ostringstream::out);
    cout.rdbuf(oss.rdbuf()); \$self->list(); cout.rdbuf(old);
    return oss.str();}
%enddef
%define __STR_LIST_OSS__()
  std::string __str__() {
    std::ostringstream oss(std::ostringstream::out); \$self->list(oss);
    return oss.str();}
%enddef
%define __STR_PRINT__()
  std::string __str__() {
    std::streambuf* old = cout.rdbuf();
    std::ostringstream oss(std::ostringstream::out);
    cout.rdbuf(oss.rdbuf()); \$self->print(); cout.rdbuf(old);
    return oss.str();}
%enddef
%define __STR_PRINT_OSS__()
  std::string __str__() {
    std::ostringstream oss(std::ostringstream::out); \$self->print(oss);
    return oss.str();}
%enddef

// Class extensions.
%extend Pythia8::BeamParticle {
  __STR_LIST_OSS__()
  Pythia8::ResolvedParton *__getitem__(int i) {
  if (i >= \$self->size()) {
    BEAMPARTICLE_ERROR = 1; return 0;} return &(*(\$self))[i];}
}
%extend Pythia8::CellJet {__STR_LIST_OSS__()}
%extend Pythia8::Clustering {__STR_LIST__()}
%extend Pythia8::ClusterJet {__STR_LIST_OSS__()}
%extend Pythia8::ColConfig {
  __STR_LIST_OSS__()
  Pythia8::ColSinglet *__getitem__(int i) {
    if (i >= \$self->size()) {COLCONFIG_ERROR = 1; return 0;}
    return &(*(\$self))[i];}
}
%extend Pythia8::ColourDipole {__STR_PRINT__()}
%extend Pythia8::ColourJunction {__STR_PRINT__()}
%extend Pythia8::ColourParticle {__STR_LIST__()}
%extend Pythia8::Event {
  __STR_LIST_OSS__()
  Pythia8::Particle *__getitem__(int i) {
    if (i >= \$self->size()) {EVENT_ERROR = 1; return 0;}
    return &(*(\$self))[i];}
}
%extend Pythia8::GammaMatrix {
  __STR_OSS__()
  Pythia8::GammaMatrix __rmul__(std::complex<double> s) {
    Pythia8::GammaMatrix g = *\$self; return s*g;}
  Pythia8::GammaMatrix __rsub__(std::complex<double> s) {
    Pythia8::GammaMatrix g = *\$self; return s-g;}
  Pythia8::GammaMatrix __radd__(std::complex<double> s) {
    Pythia8::GammaMatrix g = *\$self; return s+g;}
}
%extend Pythia8::HardProcess {__STR_LIST__()}
%extend Pythia8::Hist {
  __STR_OSS__()
  Pythia8::Hist __radd__(double f) {
    Pythia8::Hist h = *\$self; return f+h;}
  Pythia8::Hist __rsub__(double f) {
    Pythia8::Hist h = *\$self; return f-h;}
  Pythia8::Hist __rmul__(double f) {
    Pythia8::Hist h = *\$self; return f*h;}
  Pythia8::Hist __rdiv__(double f) {
    Pythia8::Hist h = *\$self; return f/h;}
}
%extend Pythia8::Info {__STR_LIST_OSS__()}
%extend Pythia8::LHAgenerator {__STR_PRINT_OSS__()}
%extend Pythia8::LHAinitrwgt {__STR_PRINT_OSS__()}
%extend Pythia8::LHArwgt {__STR_PRINT_OSS__()}
%extend Pythia8::LHAscales {__STR_PRINT_OSS__()}
%extend Pythia8::LHAweightgroup {__STR_PRINT_OSS__()}
%extend Pythia8::LHAweight {__STR_PRINT_OSS__()}
%extend Pythia8::LHAweights {__STR_PRINT_OSS__()}
%extend Pythia8::LHAwgt {__STR_PRINT_OSS__()}
%extend Pythia8::LHblock {__STR_PRINT_OSS__()}
%extend Pythia8::LHmatrixBlock {__STR_PRINT_OSS__()}
%extend Pythia8::LHtensor3Block {__STR_PRINT_OSS__()}
%extend Pythia8::ParticleData {__STR_LIST_OSS__()}
%extend Pythia8::PartonSystems {__STR_LIST_OSS__()}
%extend Pythia8::RotBstMatrix {__STR_OSS__()}
%extend Pythia8::SlowJet {__STR_LIST_OSS__()}
%extend Pythia8::SpaceShower {__STR_LIST_OSS__()}
%extend Pythia8::Sphericity {__STR_LIST_OSS__()}
%extend Pythia8::Thrust {__STR_LIST_OSS__()}
%extend Pythia8::TimeShower {__STR_LIST_OSS__()}
%extend Pythia8::Vec4 {
  __STR_OSS__()
  double __getitem__(int i) {
    if (i >= 4) {VEC4_ERROR = 1; return 0;} return (*(\$self))[i];}
  Pythia8::Vec4 __rmul__(double f) {
    Pythia8::Vec4 v = *\$self; return f*v;}
}
%extend Pythia8::Wave4 {
  __STR_OSS__()
  Pythia8::Wave4 __rmul__(std::complex<double> s) {
    Pythia8::Wave4 w = *\$self; return s*w;}
  Pythia8::Wave4 __rmul__(double s) {
    Pythia8::Wave4 w = *\$self; return s*w;}
  Pythia8::Wave4 __mul__(Pythia8::GammaMatrix g) {
    Pythia8::Wave4 w = *\$self; return w*g;}
}
%extend Pythia8::XMLTag {__STR_PRINT_OSS__()}
BLOCKTEXT

# Run SWIG and create the C++ wrapper header file.
swig -c++ -python python.cfg
cat > $HEADER_FILE << BLOCKTEXT
// Copyright (C) 2016 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten, March 2016.

// This file contains a Python interface to Pythia 8 generated with SWIG.

BLOCKTEXT
cat $CXX_FILE >> $HEADER_FILE
echo "// PYTHON SOURCE" >> $HEADER_FILE
cat >> $HEADER_FILE << BLOCKTEXT
//# Copyright (C) 2016 Torbjorn Sjostrand.
//# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
//# Please respect the MCnet Guidelines, see GUIDELINES for details.

//# This file contains a Python interface to Pythia 8 generated with SWIG.

BLOCKTEXT

# Insert documentation into the Python file and write to the header file.
python << BLOCKTEXT
files = """$FILES"""
classes = {}; methods = {}; level = [('Pythia8', 1)]; brace = 0; dstrs = '"""\n'
for f in files.split():
    f = open(f[7:])
    for l in f:
        l = l.lstrip()
        brace += l.count('{') - l.count('}')
        if brace < 1: continue
        elif '// ' in l: dstrs += '//' + l.split('// ')[-1]
        elif 'class ' in l:
            name = l.split()[1]; level += [(name, brace)];
            classes[name] = dstrs + '//"""\n'; dstrs = '"""\n'
        elif '(' in l:
            name = l.split('(')[0].split()
            if len(name) == 0: continue
            else: name = name[-1]
            if not level[-1][0] in methods: methods[level[-1][0]] = {}
            methods[level[-1][0]][name] = dstrs + '//"""\n'; dstrs = '"""\n'
        if brace < level[-1][1]: level.pop(); dstrs = '"""\n'
level = 'Pythia8'; i = open('$PYTHON_FILE'); o = open('$HEADER_FILE', 'a')
for l in i:
    try: name = l.split()[1].split('(')[0]
    except: name = None
    if name == '__init__': name = level
    if not name: o.write('//' + l)
    elif l.startswith('class'):
        level = name; o.write('//' + l)
        if name in classes: o.write('//    ' + classes[name])
    elif l.startswith('def'):
        level = 'Pythia8'; o.write('//' + l)
        if level in methods and name in methods[level]:
            o.write('//  ' + methods[level][name])
    elif l.startswith('    def') and ':' in l:
        if level in methods and name in methods[level]:
            l = l.split(':')
            if l[1].lstrip() == '':
                o.write('//' + l[0] + ':\n')
                o.write('//        ' + methods[level][name])
            else:
                o.write('//' + l[0] + ':\n')
                o.write('//        ' + methods[level][name])
                o.write('//        ' + l[1].lstrip())
        else: o.write('//' + l)
    else: o.write('//' + l)
BLOCKTEXT
rm -rf $CFG_FILE $CXX_FILE $PYTHON_FILE python
