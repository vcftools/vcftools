/*
 * bcf_entry_setters.cpp
 *
 *  Created on: Sep 20, 2012
 *      Author: Anthony Marcketta
 *      ($Revision: 1 $)
 */

#include "bcf_entry.h"

void bcf_entry::set_ALT(const int n_allele)
{
	ALT.resize(n_allele-1);
	unsigned int pos = ALT_pos;
	string allele;
	for (int ui=0; ui<(n_allele-1); ui++)
	{
		allele = get_typed_string( &pos, line );
		std::transform(allele.begin(), allele.end(), allele.begin(), ::toupper);
		ALT[ui] = allele;
	}

	parsed_ALT = true;
}

void bcf_entry::set_ALT(const string &in)
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

void bcf_entry::set_INFO()
{
	int key;
	unsigned int size, type, i = INFO_pos;
	string data_type;
	INFO.resize(N_info);
	bool miss = true;

	for (unsigned int ui=0; ui<N_info; ui++)
	{
		key = get_typed_int(&i, line, type, size);
		get_type(&i, line, type, size);

		pair<string, string> INFO_entry(entry_header.INFO_map[key].ID, ".");
		data_type = entry_header.INFO_map[key].Type_str;
		ostringstream ss(ostringstream::out);

		for (unsigned int uj=0; uj<size; uj++)
		{
			if (uj !=0 && type != 7)
				ss << ",";

			if ( check_missing(i, type, line) )
				ss << ".";

			else if(type==1)
			{
				int8_t tmp;
				memcpy(&tmp, &line[i], sizeof(tmp));
				i += sizeof(tmp);
				ss << (int)tmp;
				miss = false;
			}
			else if(type==2)
			{
				int16_t tmp;
				memcpy(&tmp, &line[i], sizeof(tmp));
				i += sizeof(tmp);
				ss << (int)tmp;
				miss = false;
			}
			else if(type==3)
			{
				int32_t tmp;
				memcpy(&tmp, &line[i], sizeof(tmp));
				i += sizeof(tmp);
				ss << (int)tmp;
				miss = false;
			}
			else if(type==5)
			{
				float tmp;
				memcpy(&tmp, &line[i], sizeof(tmp));
				i += sizeof(tmp);
				ss << tmp;
				miss = false;
			}
			else if(type==7)
			{
				char tmp;
				memcpy(&tmp, &line[i], sizeof(tmp));
				i += sizeof(tmp);
				ss << tmp;
				miss = false;
			}
			else
			{
				LOG.printLOG("Error: Unknown type: " + header::int2str(type) + "\n");
				parsed_INFO = false;
				return;
			}
		}

		string value = ss.str();
		if (miss)
			INFO_entry.second = ".";
		else if (data_type == "Flag")
			INFO_entry.second = "1";
		else if (value != "")
			INFO_entry.second = ss.str();
		else
			INFO_entry.second = ".";
		INFO[ui] = INFO_entry;
	}
	parsed_INFO = true;
}

void bcf_entry::set_FORMAT()
{
	FORMAT.resize(0);
	FORMAT_to_idx.clear();
	unsigned int l_pos = L_shared + 2*sizeof(uint32_t);
	unsigned int fmt_key, type, size, skip;

	for (unsigned int ui=0; ui<N_format; ui++)
	{
		fmt_key = get_typed_int(&l_pos, line, type, size);
		get_type(&l_pos, line, type, size);

		if ( (type == 1) or (type == 7) )
			skip = sizeof(int8_t)*size*N_indv;
		else if (type == 2)
			skip = sizeof(int16_t)*size*N_indv;
		else if ( (type == 3) or (type == 5) )
			skip = sizeof(int32_t)*size*N_indv;
		else
		{
			LOG.printLOG("Error: Unknown type: " + header::int2str(type) + "\n");
			parsed_FORMAT = false;
			return;
		}
		add_FORMAT_entry( entry_header.FORMAT_map[fmt_key].ID, fmt_key, ui, l_pos, type, size);
		l_pos += skip;
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

void bcf_entry::add_FORMAT_entry(const string &in, const unsigned int &fmt_key, const unsigned int &pos, const unsigned int &line_pos, const unsigned int &type, const unsigned int &size)
{
	FORMAT.push_back(in);
	FORMAT_to_idx[in] = pos;
	FORMAT_positions[pos] = line_pos;
	FORMAT_types[pos] = type;
	FORMAT_sizes[pos] = size;
	FORMAT_keys[pos] = fmt_key;

	if ( (type==1) or (type==7) )
		FORMAT_skip[pos] = ( sizeof(int8_t)*size );
	else if ( type==2 )
		FORMAT_skip[pos] = ( sizeof(int16_t)*size );
	else if ( (type==3) or (type==5) )
		FORMAT_skip[pos] = ( sizeof(int32_t)*size );
}

void bcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const unsigned int &pos, const unsigned int &size)
{
	int8_t tmp, tmp2;
	unsigned int cur_pos = pos;
	char phased[2] = {'/', '|'};
	ploidy.resize(N_indv);
	ploidy[indv] = 0;

	for (unsigned int ui=0; ui<size; ui++)
	{
		tmp = *reinterpret_cast<const int8_t*>(&line[cur_pos]);
		if ( tmp == (int8_t)0x81 )
			break;
		ploidy[indv]++;
		cur_pos += sizeof(int8_t);
	}

	if (ploidy[indv] == 0)
	{
		set_indv_GENOTYPE_alleles(indv, make_pair(-2, -2));
	}
	else if (ploidy[indv] == 1)
	{
		set_indv_PHASE(indv, '|');
		tmp = *reinterpret_cast<const int8_t*>(&line[pos]);
		if (tmp == (int8_t)0x80)
			tmp = -1;
		else
			tmp = (tmp >> 1) - 1;
		set_indv_GENOTYPE_alleles(indv, make_pair(tmp, -2));
	}
	else if (ploidy[indv] == 2)
	{
		tmp = *reinterpret_cast<const int8_t*>(&line[pos]);
		tmp2 = *reinterpret_cast<const int8_t*>(&line[pos+sizeof(int8_t)]);

		if (tmp == (int8_t)0x80)
			tmp = -1;
		else
			tmp = (tmp >> 1) - 1;

		if (tmp2 == (int8_t)0x80)
		{
			tmp2 = -1;
			set_indv_PHASE(indv, '/');
		}
		else
		{
			char phase = phased[ tmp2 & (int8_t)1 ];
			tmp2 = (tmp2 >> 1) - 1;
			set_indv_PHASE(indv, phase);
		}
		set_indv_GENOTYPE_alleles(indv, make_pair((int)tmp, (int)tmp2));
	}
	else if (ploidy[indv] > 2)
		LOG.error("Polyploidy found, and is not supported by vcftools: " + CHROM + ":" + header::int2str(POS));
	parsed_GT[indv] = true;
}

void bcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase)
{
	set_indv_GENOTYPE_ids(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void bcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase)
{
	pair<int, int> a(-1,-1);
	if (genotype.first != ".")
		a.first = header::str2int(genotype.first);
	if (genotype.second != ".")
		a.second = header::str2int(genotype.second);

	set_indv_GENOTYPE_alleles(indv, a);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void bcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, const pair<int, int> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (in.first == 0x81)
		a.first = -2;
	else if (in.first != 0x80)
		a.first = in.first;

	if (in.second == 0x81)
		a.second = -2;
	else if (in.second != 0x80)
		a.second = in.second;

	GENOTYPE[indv] = in;
	parsed_GT[indv] = true;
}

void bcf_entry::set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-2,-2));
	GENOTYPE[indv] = in;
}

void bcf_entry::set_indv_PHASE(unsigned int indv, char in)
{
	if (PHASE.size() == 0)
		PHASE.resize(N_indv, '/');

	PHASE[indv] = in;
	parsed_GT[indv] = true;
}

void bcf_entry::set_indv_GQUALITY(unsigned int indv, const vector<char> &in)
{
	float tmp;
	memcpy(&tmp, &in[0], sizeof(tmp));

	parsed_GQ[indv] = true;
	if (tmp == 0x7F800001)
	{
		if (GQUALITY.size() > 0)
			GQUALITY[indv] = -1;
		return;
	}
	if (GQUALITY.size() == 0)
		GQUALITY.resize(N_indv, -1);

	if (tmp > 99.0)
		tmp = 99;
	GQUALITY[indv] = tmp;
}

void bcf_entry::set_indv_GQUALITY(unsigned int indv, const float &in)
{
	parsed_GQ[indv] = true;
	if ( (in == -1) or (in == 0x7F800001) )
	{
		if (GQUALITY.size() > 0)
			GQUALITY[indv] = -1;
		return;
	}
	if (GQUALITY.size() == 0)
		GQUALITY.resize(N_indv, -1);

	if (in > 99)
		GQUALITY[indv] = 99;
	else
		GQUALITY[indv] = in;
}

void bcf_entry::set_indv_GFILTER(unsigned int indv, const vector<char> &in)
{
	parsed_FT[indv] = true;

	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	GFILTER[indv].resize(0);
	if (in.empty())
		return;
	else if ((in.size() == 1) and (in[0] == '\0') )
		return;

	ostringstream ss;
	string ith_FILTER;
	ss.clear();
	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		if (in[ui] == ';')
		{
			ith_FILTER = ss.str();
			ss.clear();

			if ((ith_FILTER.size()==0) || (ith_FILTER == "."))
				continue;	// Don't bother storing "unfiltered" state.

			GFILTER[indv].push_back(ith_FILTER);
		}
		else
			ss << in[ui];
	}
	ith_FILTER = ss.str();
	ss.clear();

	if ((ith_FILTER.size()!=0) || (ith_FILTER != "."))
		GFILTER[indv].push_back(ith_FILTER);
}

void bcf_entry::set_indv_GFILTER(unsigned int indv, const string &in)
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

void bcf_entry::set_FILTER()
{
	FILTER_str = get_int_vector( &FILTER_pos, line );
	FILTER.resize(0);

	if ( (FILTER_str.size() == 1) and entry_header.FILTER_map[ FILTER_str[0] ].ID == "")
	{
		FILTER.push_back("PASS");
		return;
	}

	for (unsigned int ui=0; ui<FILTER_str.size(); ui++)
	{
		if ((int8_t)FILTER_str[ui] != (int8_t)0x80)
			FILTER.push_back( entry_header.FILTER_map[ FILTER_str[ui] ].ID );
		else
			FILTER.push_back( "." );
	}
	sort(FILTER.begin(), FILTER.end());
	parsed_FILTER = true;
}
