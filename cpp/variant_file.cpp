/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	return N_kept_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   int i = 0;
   int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_contigs(const string &contigs_file, vector<string> &contig_vector)
{
	if (contigs_file == "")
		LOG.error("Contig declarations in header are necessary for BCF conversion. Use --contigs <filename> to add contigs to the header.");

	ifstream contigs(contigs_file.c_str());
	if (!contigs.is_open())
		LOG.error("Could not open contigs file: " + contigs_file);

	string line;
	int contig_lines = 0;
	contig_vector.resize(0);

	while (getline(contigs, line))
	{
		if (line.find("##contig=")==string::npos)
			LOG.error("Contigs file must contain only contig header lines.");

		contig_vector.push_back(line);
		contig_lines++;
	}

	contigs.close();
	LOG.printLOG("Including "+header::int2str(contig_lines)+" header lines from the contig file.\n");
}

void variant_file::read_temp_site(ifstream &tmp_file, string &CHROM, int &POS, vector< pair<int,int> > &GTs)
{
	stringstream chr;
	char tmp_char;
	while(true)
	{
		tmp_file.read(&tmp_char,sizeof(char));
		if (tmp_char == '\n')
			break;
		chr << tmp_char;
	}
	CHROM = chr.str();

	tmp_file.read((char*)&POS,sizeof(POS));

	char in_byte, tmp_gt;
	for(unsigned int ui=0; ui<GTs.size(); ui++)
	{
		tmp_file.read(&in_byte,sizeof(in_byte));
		tmp_gt = in_byte & 0x03;
		if (tmp_gt == 0x02)
			GTs[ui].second = -1;
		else
			GTs[ui].second = (int)tmp_gt;
		in_byte = in_byte >> 4;
		tmp_gt = in_byte & 0x03;
		if (tmp_gt == 0x02)
			GTs[ui].first = -1;
		else
			GTs[ui].first = (int)tmp_gt;
	}
}

void variant_file::read_big_temp_site(ifstream &tmp_file, string &CHROM, int &POS, int &alleles, vector< pair<int,int> > &GTs)
{
	stringstream chr;
	char tmp_char;
	while(true)
	{
		tmp_file.read(&tmp_char,sizeof(char));
		if (tmp_char == '\n')
			break;
		chr << tmp_char;
	}
	CHROM = chr.str();

	tmp_file.read((char*)&POS,sizeof(POS));

	int8_t tmp_alleles;
	tmp_file.read((char*)&tmp_alleles,sizeof(tmp_alleles));
	alleles = (int)tmp_alleles;

	char in_byte = 0xFF;
	for(unsigned int ui=0; ui<GTs.size(); ui++)
	{
		tmp_file.read(&in_byte,sizeof(in_byte));
		if (in_byte == (char)0xFF)
			GTs[ui].first = -1;
		else
			GTs[ui].first = (int)in_byte;

		tmp_file.read(&in_byte,sizeof(in_byte));
		if (in_byte == (char)0xFF)
			GTs[ui].second = -1;
		else
			GTs[ui].second = (int)in_byte;
	}
}
