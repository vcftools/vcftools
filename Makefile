#
# 1) Compiling
#   Type make in this directory
#
# 2) Installation
#   Edit the BINDIR and MODDIR below as necessary or pass the PREFIX variable
#   to the make command. When not set, the programs will be placed in "bin" 
#   and "lib" subdirectories in this directory.
#       PREFIX="/install/to/path/prefix" make install
#
#   Add the MODDIR to your PERL5LIB environment variable:
#       export PERL5LIB=${PREFIX}/lib:${PERL5LIB}
#
#   Add the MANDIR to your MANPATH environment variable:
#       export MANPATH=${PREFIX}/bin:$MANPATH
#

export SRCDIR = $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
ifndef PREFIX
    export PREFIX = ${SRCDIR}
endif
export BINDIR = ${PREFIX}/bin
export MANDIR = ${PREFIX}/bin/man1
export MODDIR = ${PREFIX}/lib/perl5/site_perl

DIRS = cpp perl
install:
	    @mkdir -p $(BINDIR); mkdir -p $(MODDIR); mkdir -p $(MANDIR); \
	    cp ${SRCDIR}/cpp/vcftools.1 $(MANDIR); \
        for dir in $(DIRS); do cd $$dir && $(MAKE) $(MAKEFLAGS) && cd ..; done

clean:
		@for dir in $(DIRS); do cd $$dir && $(MAKE) clean && cd ..; done
