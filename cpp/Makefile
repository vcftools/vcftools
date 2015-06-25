# Make file for vcftools
# Author: Adam Auton
# ($Revision: 230 $)

# Compiler
CC = gcc
CPP = g++
# Output executable
EXECUTABLE = vcftools
# Flag used to turn on compilation of PCA routines
ifndef VCFTOOLS_PCA
	VCFTOOLS_PCA = 0
endif
# Compiler flags
CFLAGS = -O2 -m64
#CFLAGS = -Wall -O2 -pg -m64
CPPFLAGS = -O2 -D_FILE_OFFSET_BITS=64 
#CPPFLAGS = -O2 -Wall -pg -D_FILE_OFFSET_BITS=64
# Included libraries (zlib)
LIB = -lz 
#LIB = -lz -I/opt/local/include/ -L/opt/local/lib/

OBJS = vcftools.o bcf_file.o vcf_file.o variant_file.o \
	header.o bcf_entry.o vcf_entry.o entry.o entry_getters.o entry_setters.o \
	vcf_entry_setters.o bcf_entry_setters.o entry_filters.o variant_file_filters.o variant_file_output.o parameters.o \
	variant_file_format_convert.o variant_file_diff.o output_log.o bgzf.o gamma.o

ifeq ($(VCFTOOLS_PCA), 1)
	# Define flag for PCA routine compilation
	CPPFLAGS += -DVCFTOOLS_PCA
	# Add LAPACK library
	LIB += -llapack	
	# Add PCA source code
	OBJS+= dgeev.o
endif

vcftools: $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) -o vcftools $(LIB)
ifdef BINDIR
	cp $(CURDIR)/$@ $(BINDIR)/$@
endif

bgzf: bgzf.c
	$(CC) -c $(CFLAGS) $(FLAGS) bgzf.c $(LIB) -o bgzf.o

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPP) -MM $(CPPFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f vcftools *.o *.d
	@rm -f $(BINDIR)/vcftools
