# Make file for vcftools
# Author: Adam Auton
# ($Revision: 230 $)

# Compiler
CC = gcc
CXX = g++

# Flag used to turn on compilation of PCA routines
ifndef VCFTOOLS_PCA
	VCFTOOLS_PCA = 0
endif

# Compiler flags
CFLAGS += -O2 -D_FILE_OFFSET_BITS=64
CXXFLAGS += -O2 -D_FILE_OFFSET_BITS=64

# Included libraries (zlib)
LIBS = -lz

OBJS = vcftools.o bcf_file.o vcf_file.o variant_file.o header.o bcf_entry.o     \
	vcf_entry.o entry.o entry_getters.o entry_setters.o vcf_entry_setters.o \
	bcf_entry_setters.o entry_filters.o variant_file_filters.o              \
	variant_file_output.o parameters.o variant_file_format_convert.o        \
	variant_file_diff.o output_log.o bgzf.o gamma.o

ifeq ($(VCFTOOLS_PCA), 1)
	# Define flag for PCA routine compilation
	CXXFLAGS += -DVCFTOOLS_PCA
	# Add LAPACK library
	LIBS += -llapack
	# Add PCA source code
	OBJS += dgeev.o
endif

vcftools: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJS) -o vcftools $(LDFLAGS) $(LIBS)
ifdef BINDIR
	cp $(CURDIR)/$@ $(BINDIR)/$@
endif

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $*.cpp -o $*.o
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f vcftools *.o *.d
	@rm -f $(BINDIR)/vcftools
