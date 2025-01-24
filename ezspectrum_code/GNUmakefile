PROGRAM = ezFCF

# If your armadillo is installed in a non-default location uncomment and adjust
# approprietly the next line.
# ARMADILLO_ROOT = /path/to/armadillo
# And uncomment the next two lines as well.
# INCLUDES_ARMADILLO = -I$(ARMADILLO_ROOT)/include
# LDFLAGS = -L$(ARMADILLO_ROOT)/lib64

LDLIBS = -larmadillo
CXX_WNO_DEPRECATED = -Wno-deprecated
INCLUDES = $(INCLUDES_ARMADILLO) 

CXX = g++
CC = gcc
CXXFLAGS = -O3 --std=c++20
# CXXFLAGS = --std=c++20 -g
CPPFLAGS = $(INCLUDES)

mol_prop_src = $(wildcard molecular_prop/*.C)
mol_prop_obj = $(patsubst %.C,%.o,$(mol_prop_src))

basic_methods_src = $(wildcard basic_methods/*.C)
basic_methods_obj = $(basic_methods_src:.C=.o)

harmonic_pes_src =$(wildcard harmonic_pes/*.C)
harmonic_pes_obj =$(harmonic_pes_src:.C=.o)

sources = main.C $(mol_prop_src) $(basic_methods_src) $(harmonic_pes_src)
objects = $(sources:.C=.o)

INCLUDES += $(addprefix -I ,molecular_prop basic_methods harmonic_pes)

SHELL = /bin/sh
.DELETE_ON_ERROR:

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): $(objects)
	$(CXX) $^ $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -o $(PROGRAM)

# Add headers as prerequisites
%.d: %.C
	set -e; rm -f $@; \
	 $(CXX) -MM $(CPPFLAGS) $(INCLUDES) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$
include $(objects:.o=.d)

# Common prefixes for installation directories.
# Allows the user a modification of the default install paths, e.g.,
#  make prefix=/home/user/.local install
prefix=/usr/local
exec_prefix=$(prefix)
bindir=$(exec_prefix)/bin
datarootdir=$(prefix)/share
datadir=$(datarootdir)

# Variables for conveniece
# DESTDIR is for the user's use, e.g,
#  make DESTDIR=/tmp/staged/ install
ezFCF_bindir = $(DESTDIR)$(bindir)
ezFCF_datadir = $(DESTDIR)$(datadir)/$(PROGRAM)

INSTALL = install -c
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA = $(INSTALL) -m 644

.PHONY: install
install: installdirs
	$(INSTALL_PROGRAM) $(PROGRAM) $(ezFCF_bindir)
	$(INSTALL_DATA) atomicMasses.xml $(ezFCF_datadir)

.PHONY: installdirs
installdirs:
	-mkdir -p $(ezFCF_datadir)
	-mkdir -p $(ezFCF_bindir)

# If some prefix or DESTDIR variables were used during installation
# then these variables must be specified also during uninstallation.
.PHONY: uninstall
uninstall:
	rm -f $(ezFCF_bindir)/$(PROGRAM)
	rm -f $(ezFCF_datadir)/atomicMasses.xml
	rmdir $(ezFCF_datadir)

.PHONY: mostlyclean
mostlyclean:
	rm -f $(PROGRAM) *.o $(objects)	

.PHONY: clean
clean: mostlyclean
	rm -f $(objects:.o=.d)

.PHONY: distclean
distclean: clean
