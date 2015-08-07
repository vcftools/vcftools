/*
	entry_getters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "entry.h"

// Return the CHROMosome name
string entry::get_CHROM() const
{
	return CHROM;
}

// Return the CHROMosome name
void entry::get_CHROM(string &out) const
{
	out = CHROM;
}

int entry::get_POS() const
{
	return POS;
}

string entry::get_ID() const
{
	if (ID.size() == 0)
		return ".";
	return ID;
}

string entry::get_REF() const
{
	if (REF == "")
		return ".";
	else
		return REF;
}

string entry::get_ALT() const
{
	assert(parsed_ALT == true);

	string out;
	if (ALT.empty())
		out = ".";
	else if (ALT.size() == 1 && ALT[0] == "")
		out = ".";
	else
	{
		out = ALT[0];
		for (unsigned int ui=1; ui<ALT.size(); ui++)
			out += "," + ALT[ui];
	}
	return out;
}

bool entry::is_SNP() const
{
	assert(parsed_ALT == true);

	if (REF.size() != 1)
		return false;	// Reference isn't a single base

	if (ALT.empty())
		return false;	// No alternative allele

	for (unsigned int ui=0; ui<ALT.size(); ui++)
		if (ALT[ui].size() != 1)
			return false;	// Alternative allele isn't a single base

	return true;
}

bool entry::is_biallelic_SNP() const
{
	assert(parsed_ALT == true);

	if (REF.size() != 1)
		return false;	// Reference isn't a single base

	if (ALT.size() != 1)
		return false;	// Not biallelic

	if (ALT[0].size() != 1)
		return false;	// Alternative allele isn't a single base

	return true;
}

bool entry::is_diploid() const
{
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			if (ploidy[ui] != 2)
				return false;
		}
	}
	return true;
}

void entry::get_allele(int allele_num, string &out) const
{
	assert(parsed_ALT == true);

	if (allele_num == -2)
		out = "";
	else if (allele_num == 0)
		out = REF;
	else if ((allele_num == -1) || (unsigned(allele_num - 1) >=  ALT.size()))
		out = ".";
	else
		out = ALT[allele_num-1];
}

string entry::get_allele(int allele_num) const
{
	assert(parsed_ALT == true);

	if (allele_num == -2)
		return "";
	else if (allele_num == 0)
		return REF;
	else if ((allele_num < 0) || (unsigned(allele_num - 1) >=  ALT.size()))
		return ".";
	else
		return ALT[allele_num-1];
}

string entry::get_ALT_allele(int allele_num) const
{
	assert(parsed_ALT == true);

	if (allele_num == -2)
		return "";
	else if ((allele_num == -1) || (unsigned(allele_num) >=  ALT.size()))
		return ".";
	return ALT[allele_num];
}

void entry::get_alleles_vector(vector<string> &out) const
{
	assert(parsed_ALT == true);
	out.resize(ALT.size()+1);
	out[0] = REF;
	copy(ALT.begin(), ALT.end(), out.begin()+1);
}

double entry::get_QUAL() const
{
	return QUAL;
}

string entry::get_FILTER() const
{
	assert(parsed_FILTER == true);

	ostringstream out;
	if (FILTER.empty())
		out << ".";
	else
	{
		out << FILTER[0];
		for (unsigned int ui=1; ui<FILTER.size(); ui++)
			out << ";" << FILTER[ui];
	}
	return out.str();
}

void entry::get_FILTER_vector(vector<string> &out) const
{
	assert(parsed_FILTER == true);
	out = FILTER;
}

string entry::get_INFO(const set<string> &INFO_to_keep, bool keep_all_INFO) const
{
	assert(parsed_INFO == true);

	ostringstream sout;
	sout.str("");
	sout.clear();

	bool first=true;
	if ( ( (!INFO.empty()) && (!INFO_to_keep.empty())  ) || keep_all_INFO )
	{
		string key;
		for (unsigned int ui=0; ui<INFO.size();ui++)
		{
			key = INFO[ui].first;
			if ( keep_all_INFO or (INFO_to_keep.find(key) != INFO_to_keep.end() ) )
			{
				if (first != true)
					sout << ";";

				sout << key;
				if (INFO[ui].second != "")
					sout << "=" << INFO[ui].second;
				first = false;
			}
		}
	}

	if (first == true)
	{	// Didn't find any INFO fields to keep
		sout.str(".");
	}
	return sout.str();
}

vector<pair<string, string> > entry::get_INFO_vector(const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	assert(parsed_INFO == true);

	vector<pair<string, string> > out_vector;
	if (keep_all_INFO == true)
		return INFO;

	if ( (!INFO.empty()) && (!INFO_to_keep.empty()) )
	{
		string key;
		for (unsigned int ui=0; ui<INFO.size();ui++)
		{
			key = INFO[ui].first;
			if ( keep_all_INFO or (INFO_to_keep.find(key) != INFO_to_keep.end() ) )
				out_vector.push_back( INFO[ui] );
			else
				N_INFO_removed++;
		}
	}

	return out_vector;
}

string entry::get_INFO_value(const string &key) const
{
	assert(parsed_INFO == true);

	for (unsigned int ui=0; ui<INFO.size(); ui++)
	{
		if (INFO[ui].first == key)
			return INFO[ui].second;
	}
	return "?";
}

vector<string> entry::get_INFO_values(const string &key) const
{
	vector<string> out;
	string tmp;

	tmp = get_INFO_value(key);
	if (tmp != "?")	header::tokenize(tmp, ',', out);

	return out;
}

string entry::get_FORMAT() const
{
	assert(parsed_FORMAT == true);

	string out;
	bool first = true;
	for (unsigned int ui=0; ui<FORMAT.size(); ui++)
	{
		if (first == false)
			out += ":";
		out += FORMAT[ui];
		first = false;
	}
	return out;
}

void entry::get_FORMAT_binary(vector<char> &out) const
{
	assert(parsed_FORMAT_binary == true);
	out = FORMAT_binary;
}

// Return the alleles of a genotype as a pair of strings.
void entry::get_indv_GENOTYPE_strings(unsigned int indv, pair<string, string> &out) const
{
	assert(parsed_GT[indv] == true);

	static string out_allele1, out_allele2;

	get_allele(GENOTYPE[indv].first, out_allele1);
	get_allele(GENOTYPE[indv].second, out_allele2);
	out = make_pair(out_allele1, out_allele2);
}

void entry::get_indv_GENOTYPE_ids(unsigned int indv, pair<int, int> &out) const
{
	assert(parsed_GT[indv] == true);
	out = GENOTYPE[indv];
}

char entry::get_indv_PHASE(unsigned int indv) const
{
	assert(parsed_GT[indv] == true);
	return PHASE[indv];
}

int entry::get_indv_DEPTH(unsigned int indv) const
{
	assert(parsed_DP[indv] == true);
	if (DEPTH.empty())
		return -1;
	return DEPTH[indv];
}

double entry::get_indv_GQUALITY(unsigned int indv) const
{
	assert(parsed_GQ[indv] == true);
	if (GQUALITY.empty())
		return -1;
	return GQUALITY[indv];
}

void entry::get_indv_GFILTER_vector(unsigned int indv, vector<string> &out) const
{
	assert(parsed_FT[indv] == true);
	if (!GFILTER.empty())
		out = GFILTER[indv];
	else
		out.resize(0);
}

void entry::get_indv_GFILTER(unsigned int indv, string &out) const
{
	assert(parsed_FT[indv] == true);

	if ((!GFILTER.empty()) && (GFILTER[indv].size()>0))
	{
		out="";
		for (unsigned int ui=0; ui<GFILTER[indv].size(); ui++)
		{
			if (ui!=0)
				out += ";";
			out += GFILTER[indv][ui];
		}
	}
	else
		out = ".";
}

int entry::get_indv_ploidy(unsigned int indv) const
{
	assert (parsed_GT[indv]==true);
	return ploidy[indv];
}

bool entry::FORMAT_id_exists(const string &FORMAT_id)
{
	assert(parsed_FORMAT == true);
	if (FORMAT_to_idx.find(FORMAT_id) != FORMAT_to_idx.end())
		return true;
	return false;
}

unsigned int entry::get_N_alleles() const
{
	assert(parsed_ALT == true);
	return (ALT.size()+1);
}

unsigned int entry::get_N_chr() const
{
	unsigned int out=0;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == true)
		{
			assert(parsed_GT[ui] == true);
			out += ploidy[ui];
		}
	}
	return out;
}

// Return the frequency (counts) of each allele.
void entry::get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out) const
{
	get_allele_counts(out, N_non_missing_chr_out, include_indv, include_genotype);
}

// Return the frequency (counts) of each allele.
void entry::get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out, const vector<bool> &include_indv, const vector<bool> &include_genotype) const
{
	pair<int,int> genotype;
	vector<int> allele_counts(get_N_alleles(), 0);
	N_non_missing_chr_out = 0;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		//FILTERING BY INDIVIDUAL
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			get_indv_GENOTYPE_ids(ui, genotype);

			if (genotype.first > -1)
			{
				allele_counts[genotype.first]++;
				N_non_missing_chr_out++;
			}
			if (genotype.second > -1)
			{
				allele_counts[genotype.second]++;
				N_non_missing_chr_out++;
			}
		}
	}
	out = allele_counts;
}

void entry::get_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const
{
	out_N_hom1 = 0; out_N_hom2 = 0; out_N_het = 0;
	pair<int, int> genotype;
	if (ALT.size() > 1)
		LOG.error("Tried to return the genotype counts of a non-biallelic SNP", 99);

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			get_indv_GENOTYPE_ids(ui, genotype);
			if ((genotype.first > -1) && (genotype.second > -1))
			{
				if (genotype.first != genotype.second)
					out_N_het++;
				else if (genotype.first == 0)
					out_N_hom1++;
				else if (genotype.first == 1)
					out_N_hom2++;
				else
					LOG.error("Unknown allele in genotype", 98);
			}
		}
	}
}

void entry::get_multiple_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, vector<unsigned int> &out_N_hom, vector<unsigned int> &out_N_het) const
{
	out_N_hom.assign(ALT.size()+1, 0);
	out_N_het.assign(ALT.size()+1, 0);
	pair<int, int> genotype;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			get_indv_GENOTYPE_ids(ui, genotype);

			for (int uj=0; uj<=(int)ALT.size(); uj++)
			{
				if ((genotype.first == uj) && (genotype.second == uj))
					out_N_hom[uj]++;
				else if (((genotype.first == uj) || (genotype.second == uj)) && (genotype.first != -1) && (genotype.second != -1))
					out_N_het[uj]++;
			}
		}
	}
}

// Return the counts of homozygote1, heterozygotes, and homozygote2
void entry::get_genotype_counts(unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const
{
	get_genotype_counts(include_indv, include_genotype, out_N_hom1, out_N_het, out_N_hom2);
}

void entry::get_POS_binary(vector<char> &out) const
{
	out.resize(sizeof(uint32_t));
	uint32_t pos = POS - 1;
	memcpy(&out[0], &pos, sizeof(pos));
}

void entry::get_rlen(vector<char> &out) const
{
	out.resize(sizeof(int32_t));
	int32_t rlen;
	if (REF != "" and REF != "." and REF != " ")
		rlen = (int32_t)REF.length();
	else
		rlen = (int32_t)0;
	memcpy(&out[0], &rlen, sizeof(rlen));
}

void entry::get_QUAL_binary(vector<char> &out) const
{
	out.resize(sizeof(float));
	float qual = (float)QUAL;
	memcpy(&out[0], &qual, sizeof(qual));
}

void entry::get_n_allele_info(vector<char> &out) const
{
	out.resize(sizeof(uint32_t));
	uint32_t n_allele_info = (uint32_t)ALT.size() + 1;
	uint32_t n_info = (uint32_t)(INFO.size()-N_INFO_removed);

	n_allele_info = n_allele_info << 16;
	n_allele_info = n_allele_info | n_info;

	memcpy(&out[0], &n_allele_info, sizeof(n_allele_info));
}

void entry::get_n_fmt_sample(vector<char> &out) const
{
	out.resize(sizeof(uint32_t));
	uint32_t n_fmt_sample = (uint32_t)(FORMAT.size()-N_FORMAT_removed);
	uint32_t n_sample = (uint32_t)N_indv;

	n_fmt_sample = n_fmt_sample << 24;
	n_fmt_sample = n_fmt_sample | n_sample;

	memcpy(&out[0], &n_fmt_sample, sizeof(n_fmt_sample));
}

void entry::get_ID_binary(vector<char> &out)
{
	make_typed_string(out, ID, true );
}

void entry::get_ALLELES_binary(vector<char> &out)
{
	vector<char> tmp;
	out.resize(0);

	make_typed_string(tmp, REF, true );
	out.insert(out.end(), tmp.begin(), tmp.end());

	for (unsigned int ui=0; ui<ALT.size(); ui++)
	{
		tmp.resize(0);
		make_typed_string(tmp, ALT[ui], true );
		out.insert(out.end(), tmp.begin(), tmp.end());
	}
}
