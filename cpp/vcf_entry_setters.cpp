/*
 * vcf_entry_setters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"
#include "entry.h"

void vcf_entry::set_ALT(const string &in)
{
	istringstream ss(in);
	string tmpstr;
	ALT.resize(0);
	while(!ss.eof())
	{
		getline(ss, tmpstr, ',');
		add_ALT_allele(tmpstr);
	}
	parsed_ALT = true;
}

void vcf_entry::set_FORMAT(const string &in)
{
	FORMAT.resize(0);
	FORMAT_to_idx.clear();

	if (in.size() > 0)
	{
		istringstream ss(in);
		string tmpstr;

		unsigned int pos=0;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ':');
			add_FORMAT_entry(tmpstr, pos);
			pos++;
		}
	}

	GT_idx = -1;
	GQ_idx = -1;
	DP_idx = -1;
	FT_idx = -1;

	if (FORMAT_to_idx.find("GT") != FORMAT_to_idx.end())
		GT_idx = FORMAT_to_idx["GT"];
	if (FORMAT_to_idx.find("GQ") != FORMAT_to_idx.end())
		GQ_idx = FORMAT_to_idx["GQ"];
	if (FORMAT_to_idx.find("DP") != FORMAT_to_idx.end())
		DP_idx = FORMAT_to_idx["DP"];
	if (FORMAT_to_idx.find("FT") != FORMAT_to_idx.end())
		FT_idx = FORMAT_to_idx["FT"];

	parsed_FORMAT = true;
}

void vcf_entry::add_FORMAT_entry(const string &in, unsigned int pos)
{
	FORMAT.push_back(in);
	FORMAT_to_idx[in] = pos;
}

// The following function reads in a genotype from a '0/1'-like string.
// Should handle haploid types to, but NOT polyploidy.
void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const string &in)
{
	ploidy.resize(N_indv);
	if ((in.size() == 3) && ((in.c_str()[1] == '/') || (in.c_str()[1] == '|')))
	{	// Fast, diploid case...
		ploidy[indv] = 2;
		set_indv_PHASE(indv, in.c_str()[1]);
		set_indv_GENOTYPE_alleles(indv, in.c_str()[0], in.c_str()[2]);
	}
	else
	{	// More complex case...
		size_t pos = in.find_first_of("/|");
		if (pos != string::npos)
		{	// autosome
			ploidy[indv] = 2;
			set_indv_PHASE(indv, in[pos]);
			set_indv_GENOTYPE_alleles(indv, make_pair(in.substr(0,pos), in.substr(pos+1)));
		}
		else
		{	// Male chrX, or chrY
			ploidy[indv] = 1;
			set_indv_PHASE(indv, '|');
			set_indv_GENOTYPE_alleles(indv, make_pair(in.substr(0,pos), "."));
		}

		// Check for polypoidy
		size_t pos2 = in.find_last_of("/|");
		if (pos != pos2)
			LOG.error("Polyploidy found, and not supported by vcftools: " + CHROM + ":" + header::int2str(POS));
	}

	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase)
{
	ploidy.resize(N_indv);
	ploidy[indv] = 2;
	set_indv_GENOTYPE_ids(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase)
{
	ploidy.resize(N_indv);
	ploidy[indv] = 2;
	set_indv_GENOTYPE_alleles(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, const pair<string, string> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (in.first != ".")
		a.first = header::str2int(in.first);

	if (in.second != ".")
		a.second = header::str2int(in.second);

	GENOTYPE[indv] = a;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, char a1, char a2)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (a1 != '.')
		a.first = a1 - '0';

	if (a2 != '.')
		a.second = a2 - '0';

	GENOTYPE[indv] = a;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));
	GENOTYPE[indv] = in;
}

void vcf_entry::set_indv_PHASE(unsigned int indv, char in)
{
	if (PHASE.size() == 0)
		PHASE.resize(N_indv, '/');

	PHASE[indv] = in;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GQUALITY(unsigned int indv, double in)
{
	parsed_GQ[indv] = true;
	if (in == -1)
	{
		if (GQUALITY.size() > 0)
			GQUALITY[indv] = -1;
		return;
	}
	if (GQUALITY.size() == 0)
		GQUALITY.resize(N_indv, -1);

	if (in > 99)
		in = 99;
	GQUALITY[indv] = in;
}

void vcf_entry::set_indv_GFILTER(unsigned int indv, const string &in)
{
	parsed_FT[indv] = true;

	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	GFILTER[indv].resize(0);
	if ((in.size() == 0) || (in == "."))
		return;

	static istringstream ss;
	static string ith_FILTER;
	ss.clear();
	ss.str(in);
	while (!ss.eof())
	{
		getline(ss, ith_FILTER, ';');

		if ((ith_FILTER.size()==0) || (ith_FILTER == "."))
			continue;	// Don't bother storing "unfiltered" state.

		GFILTER[indv].push_back(ith_FILTER);
	}
}

void vcf_entry::set_FILTER(const string &FILTER_str)
{
	FILTER.resize(0);

	if (FILTER_str != ".")
	{
		istringstream ss(FILTER_str);
		string ith_FILTER;
		while (!ss.eof())
		{
			getline(ss, ith_FILTER, ';');
			FILTER.push_back(ith_FILTER);
		}
	}
	sort(FILTER.begin(), FILTER.end());
	parsed_FILTER = true;
}

void vcf_entry::set_INFO(const string &INFO_str)
{
	INFO.resize(0);
	if ((INFO_str.size() > 0) && (INFO_str != "."))
	{
		istringstream ss(INFO_str);
		string tmpstr;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ';');

			istringstream ss2(tmpstr);
			getline(ss2, tmpstr, '=');
			pair<string, string> INFO_entry(tmpstr, ".");

			if (!ss2.eof())
			{	// If there is a value entry, read it now
				getline(ss2, tmpstr);
				INFO_entry.second = tmpstr;
			}
			else	// Otherwise, set it equal to 1
				INFO_entry.second = "1";

			INFO.push_back(INFO_entry);
		}
	}
	parsed_INFO = true;
}
