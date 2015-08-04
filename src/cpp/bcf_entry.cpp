/*
 * bcf_entry.cpp
 *
 *  Created on: Sep 20, 2012
 *      Author: Anthony Marcketta
 *      ($Revision: 1 $)
 */

#include "bcf_entry.h"

bcf_entry::bcf_entry(header &header_obj, vector<bool> &include_individual)
{
	N_indv = header_obj.N_indv;
	include_indv = include_individual;
	include_genotype = vector<bool>(N_indv, true);
	basic_parsed = false; fully_parsed = false;
	parsed_ALT = false; parsed_FILTER = false;
	parsed_INFO = false; parsed_FORMAT = false;
	CHROM = ""; POS = -1; REF = ""; QUAL = -1;
	N_INFO_removed = 0; N_FORMAT_removed = 0;
	passed_filters = true; parsed_FORMAT_binary = false;
	parsed_GT = vector<bool>(N_indv, false); parsed_GQ = vector<bool>(N_indv, false);
	parsed_DP = vector<bool>(N_indv, false); parsed_FT = vector<bool>(N_indv, false);
	GT_idx = -1; GQ_idx = -1; DP_idx = -1; FT_idx = -1;
	N_info = 0; N_format = 0; L_shared = 0; L_indiv = 0; line_pos = 0;
	N_allele = 0; INFO_pos = 0; FILTER_pos = 0; ALT_pos = 0; FORMAT_pos = 0;
	FORMAT_positions.resize(0); FORMAT_types.resize(0); FORMAT_sizes.resize(0); FORMAT_skip.resize(0); FORMAT_keys.resize(0);
	line.clear();

	entry_header = header_obj;
}

bcf_entry::~bcf_entry() {}

void bcf_entry::reset(const vector<char> &data_line)
{
	basic_parsed = false;
	fully_parsed = false;
	parsed_ALT = false;
	parsed_FILTER = false;
	parsed_INFO = false;
	parsed_FORMAT = false;
	parsed_FORMAT_binary = false;
	passed_filters = true;

	line = data_line;

	fill(parsed_GT.begin(), parsed_GT.end(), false);
	fill(parsed_GQ.begin(), parsed_GQ.end(), false);
	fill(parsed_DP.begin(), parsed_DP.end(), false);
	fill(parsed_FT.begin(), parsed_FT.end(), false);
	fill(include_genotype.begin(), include_genotype.end(), true);

	INFO_pos = 0; FILTER_pos = 0; ALT_pos = 0; FORMAT_pos = 0;
	FORMAT_positions.clear(); FORMAT_types.clear(); FORMAT_sizes.clear(); FORMAT_skip.clear(); FORMAT_keys.clear();
	N_INFO_removed = 0; N_FORMAT_removed = 0;
}

void bcf_entry::parse_basic_entry(bool parse_ALT, bool parse_FILTER, bool parse_INFO)
{
	if (line.empty())
	{
		if (parse_ALT)
			set_ALT("");
		return;
	}

	if (!basic_parsed)
	{
		uint32_t n_allele_info, n_fmt_sample;
		uint32_t chrom, pos, rlen;
		uint32_t shared, indiv;
		float qual;

		line_pos = 0;

		get_number(shared, &line_pos, line);
		get_number(indiv, &line_pos, line);
		L_shared = shared;
		L_indiv = indiv;

		get_number(chrom, &line_pos, line);
		get_number(pos, &line_pos, line);
		get_number(rlen, &line_pos, line);

		qual = *reinterpret_cast<const float*>(&line[line_pos]);
		line_pos += sizeof(qual);

		get_number(n_allele_info, &line_pos, line);
		get_number(n_fmt_sample, &line_pos, line);

		N_format = n_fmt_sample >> 24;
		CHROM = entry_header.CONTIG_map[chrom].ID;
		POS = pos + 1;
		ID = get_typed_string( &line_pos, line );
		REF = get_typed_string( &line_pos, line );
		QUAL = qual;

		N_allele = n_allele_info >> 16;
		N_info = n_allele_info & (uint32_t)65535;

		ALT_pos = line_pos;
		for (unsigned int ui=1; ui<N_allele; ui++)
			skip_section(&line_pos, line);

		FILTER_pos = line_pos;
		skip_section(&line_pos, line);

		INFO_pos = line_pos;
		std::transform(REF.begin(), REF.end(), REF.begin(), ::toupper);

		basic_parsed = true;
	}
	if (parse_ALT && !parsed_ALT)
		set_ALT(N_allele);
	if (parse_FILTER && !parsed_FILTER)
		set_FILTER();
	if (parse_INFO && !parsed_INFO)
		set_INFO();
}

void bcf_entry::parse_full_entry(bool parse_FORMAT)
{
	if (fully_parsed)
		return;

	if (basic_parsed == false)
		parse_basic_entry();

	FORMAT_pos = L_shared + 2*sizeof(uint32_t);
	FORMAT_positions.resize(N_format);
	FORMAT_types.resize(N_format);
	FORMAT_sizes.resize(N_format);
	FORMAT_skip.resize(N_format);
	FORMAT_keys.resize(N_format);

	if (parse_FORMAT)
		set_FORMAT();

	fully_parsed = true;
}

// Filter specific genotypes by quality
void bcf_entry::filter_genotypes_by_quality(double min_genotype_quality)
{
	if (fully_parsed == false)
		parse_full_entry();

	if (GQ_idx != -1)
	{	// Have quality info
		double quality;
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_GQ[ui] == false)
				parse_genotype_entry(ui, false, true);
			quality = get_indv_GQUALITY(ui);
			if (quality < min_genotype_quality)
				include_genotype[ui] = false;
		}
	}
}

// Set the include_genotype flag on the basis of depth
void bcf_entry::filter_genotypes_by_depth(int min_depth, int max_depth)
{
	if (fully_parsed == false)
		parse_full_entry();

	if (DP_idx != -1)
	{	// Have depth info
		int depth;
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_DP[ui] == false)
				parse_genotype_entry(ui, false, false, true);
			depth = get_indv_DEPTH(ui);
			if ((depth < min_depth) || (depth > max_depth))
				include_genotype[ui] = false;
		}
	}
}

void bcf_entry::filter_genotypes_by_filter_status(const set<string> &filter_flags_to_remove, bool remove_all)
{
	if (fully_parsed == false)
		parse_full_entry();

	vector<string> GFILTERs;
	if (FT_idx != -1)
	{	// Have GFilter info
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_FT[ui] == false)
				parse_genotype_entry(ui, false, false, false, true);
			get_indv_GFILTER_vector(ui, GFILTERs);

			if (remove_all == true)
			{	// If removing all filters, only keep things with label PASS
				if (!GFILTERs.empty())
					if ((GFILTERs[0] != "PASS") && (GFILTERs[0] != "."))
						include_genotype[ui] = false;
			}
			else
			{	// Only removing specific filters
				for (unsigned int uj=0; uj<GFILTERs.size(); uj++)
					if (filter_flags_to_remove.find(GFILTERs[uj]) != filter_flags_to_remove.end())
							include_genotype[ui] = false;
			}
		}
	}
}

void bcf_entry::parse_genotype_entry(unsigned int indv, bool GT, bool GQ, bool DP, bool FT)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT();

	unsigned int l_pos, indv_size, type, size, ui;
	static vector<int> ids;
	    ids.resize(0);

	if ( GT && !parsed_GT[indv] && GT_idx != -1 )
		ids.push_back(GT_idx);
	if (GQ && !parsed_GQ[indv] && GQ_idx != -1)
		ids.push_back(GQ_idx);
	if (DP && !parsed_DP[indv] && DP_idx != -1)
		ids.push_back(DP_idx);
	if (FT && !parsed_FT[indv] && FT_idx != -1)
		ids.push_back(FT_idx);

	for(unsigned int i=0; i<ids.size(); i++)
	{
		ui = ids[i];
		type = FORMAT_types[ui];
		size = FORMAT_sizes[ui];
		indv_size = FORMAT_skip[ui];
		l_pos = FORMAT_positions[ui] + indv*indv_size;

		if ((int)ui == GT_idx)
			set_indv_GENOTYPE_and_PHASE(indv, l_pos, size);
		else if ((int)ui == GQ_idx)
		{
			if (size>1)
				LOG.error("Error: Only expect single value for QUALITY.\n");

			float tmp;
			if (type==5)
				tmp = *reinterpret_cast<const float*>(&line[l_pos]);
			else if (type==1)
			{
				int8_t tmp2 = *reinterpret_cast<const int8_t*>(&line[l_pos]);
				tmp = (float)tmp2;
			}
			else if (type==2)
			{
				int16_t tmp2 = *reinterpret_cast<const int16_t*>(&line[l_pos]);
				tmp = (float)tmp2;
			}
			else if (type==3)
			{
				int32_t tmp2 = *reinterpret_cast<const int32_t*>(&line[l_pos]);
				tmp = (float)tmp2;
			}
			else
				LOG.error("Error: Invalid type for QUALITY.\n");

			set_indv_GQUALITY(indv, tmp);
		}
		else if ((int)ui == DP_idx)
		{
			if (size>1)
				LOG.error("Error: Only expect single value for DEPTH.\n");

			int tmp = -1;

			if (type==1)
			{
				if ( !check_missing(l_pos, 1, line) )
					tmp = *reinterpret_cast<const int8_t*>(&line[l_pos]);
			}
			else if (type==2)
			{
				if ( !check_missing(l_pos, 2, line) )
					tmp = *reinterpret_cast<const int16_t*>(&line[l_pos]);
			}
			else if (type==3)
			{
				if ( !check_missing(l_pos, 3, line) )
					tmp = *reinterpret_cast<const int32_t*>(&line[l_pos]);
			}
			else if (type==5)
			{
				float tmp2 = -1;
				if ( !check_missing(l_pos, 5, line) )
					tmp2 = *reinterpret_cast<const float*>(&line[l_pos]);
				tmp = (int)tmp2;
			}
			else
				LOG.error("Error: Invalid type for DEPTH.\n");

			set_indv_DEPTH(indv, tmp);
		}
		else if ((int)ui == FT_idx)
		{
			if (type == 7)
			{
				vector<char> tmp;
				tmp.resize( size*sizeof(char) );
				memcpy(&tmp[0], &line[l_pos], size*sizeof(char));
				set_indv_GFILTER(indv, tmp);
			}
			else
				LOG.one_off_warning("Warning: FT values must be encoded in string format.\n");
		}
	}
	// Set missing return values if requested a value, but couldn't find it
	if (GT && (parsed_GT[indv] == false))
	{
		set_indv_GENOTYPE_and_PHASE(indv, make_pair(-1,-1), '/');
	}
	if (GQ && (parsed_GQ[indv] == false))
	{
		set_indv_GQUALITY(indv, -1);
	}
	if (DP && (parsed_DP[indv] == false))
	{
		set_indv_DEPTH(indv, -1);
	}
	if (FT && (parsed_FT[indv] == false))
	{
		set_indv_GFILTER(indv, "");
	}
}

void bcf_entry::parse_genotype_entries(bool GT, bool GQ, bool DP, bool FT)
{
	for (unsigned int ui=0; ui<N_indv; ui++)
		parse_genotype_entry(ui, GT, GQ, DP, FT);
}

void bcf_entry::read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out)
{
	read_indv_generic_entry(indv, FORMAT_to_idx[FORMAT_id], out);
}

void bcf_entry::read_indv_generic_entry(unsigned int indv, const int &idx, string &out)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT();

	if(idx == GT_idx && !parsed_GT[indv])
		parse_genotype_entry(indv, true);

	out = ".";
	outstream.str("");

	string tmpstr;
	unsigned int l_pos, type, size;
	bool miss, format_miss;

	l_pos = FORMAT_positions[idx] + FORMAT_skip[idx]*indv;
	type = FORMAT_types[idx];
	size = FORMAT_sizes[idx];

	format_miss = true;
	if(type == 1)
	{
		int8_t tmp;
		if (idx == GT_idx)
		{
			pair<int, int> genotype;
			char phase;

			get_indv_GENOTYPE_ids(indv, genotype);
			phase = get_indv_PHASE(indv);
			if ((genotype.first == -2) && (genotype.second == -2))
				outstream << ".";
			else if ((genotype.first == -1) && (genotype.second == -2))
				outstream << ".";
			else if ((genotype.first > -1) && (genotype.second == -2))
				outstream << genotype.first;
			else if ((genotype.first > -1) && (genotype.second > -1))
				outstream << genotype.first << phase << genotype.second;
			else
				outstream << genotype.first << phase << genotype.second;

			out = outstream.str();
		}
		else
		{
			format_miss = true;
			for (unsigned int uj=0; uj<size; uj++)
			{
				if (check_end(l_pos, type, line))
					break;

				miss = check_missing(l_pos, type, line);
				if (uj != 0)
					outstream << ",";

				if (miss)
					outstream << ".";
				else
				{
					tmp = *reinterpret_cast<const int8_t*>(&line[l_pos]);
					outstream << int(tmp);
				}
				l_pos += sizeof(int8_t);
				format_miss = format_miss && miss;
			}

			tmpstr = outstream.str();
			if ( (tmpstr.length() > 0) and !format_miss)
				out = tmpstr;
		}
	}
	else if (type == 2)
	{
		int16_t tmp;
		format_miss = true;

		for (unsigned int uj=0; uj<size; uj++)
		{
			if (check_end(l_pos, type, line))
				break;

			miss = check_missing(l_pos, type, line);
			if (uj != 0)
				outstream << ",";

			if (miss)
				outstream << ".";
			else
			{
				tmp = *reinterpret_cast<const int16_t*>(&line[l_pos]);
				outstream << int(tmp);
			}
			l_pos += sizeof(int16_t);
			format_miss = format_miss && miss;
		}

		tmpstr = outstream.str();
		if ( (tmpstr.length() > 0) and !format_miss )
			out = tmpstr;
	}
	else if (type == 3)
	{
		int32_t tmp;
		format_miss = true;

		for (unsigned int uj=0; uj<size; uj++)
		{
			if (check_end(l_pos, type, line))
				break;

			miss = check_missing(l_pos, type, line);
			if (uj != 0)
				outstream << ",";

			if (miss)
				outstream << ".";
			else
			{
				tmp = *reinterpret_cast<const int32_t*>(&line[l_pos]);
				outstream << int(tmp);
			}
			l_pos += sizeof(int32_t);
			format_miss = format_miss && miss;
		}

		tmpstr = outstream.str();
		if ( (tmpstr.length() > 0) and !format_miss )
			out = tmpstr;
	}
	else if (type == 5)
	{
		float tmp;
		format_miss = true;

		for (unsigned int uj=0; uj<size; uj++)
		{
			if (check_end(l_pos, type, line))
				break;

			miss = check_missing(l_pos, type, line);
			if (uj != 0)
				outstream << ",";

			if (miss)
				outstream << ".";
			else
			{
				tmp = *reinterpret_cast<const float*>(&line[l_pos]);
				outstream << float(tmp);
			}
			l_pos += sizeof(float);
			format_miss = format_miss && miss;
		}

		tmpstr = outstream.str();
		if ( (tmpstr.length() > 0) and !format_miss )
			out = tmpstr;
	}
	else if (type == 7)
	{
		stringstream str_stream;
		string tmp_string;
		char tmp = '.';

		for (unsigned int uj=0; uj<size; uj++)
		{
			tmp = *reinterpret_cast<const char*>(&line[l_pos]);
			l_pos += sizeof(char);
			str_stream << tmp;
		}
		tmp_string = str_stream.str();
		tmp_string.erase( remove( tmp_string.begin(), tmp_string.end(), ' ' ), tmp_string.end() );

		if (tmp_string != "")
			out = tmp;
		else
			out = ".";
	}
}

void bcf_entry::read_all_entries(string &out)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT();

	ostringstream outstream;
	string tmpstr;
	outstream.str("");
	tmpstream.str("");
	bool format_miss, indv_miss;

	for(unsigned int ui=0; ui<N_indv; ui++)
	{
		if(include_indv[ui] == false)
			continue;
		outstream << "\t";
		indv_miss = true;

		for(unsigned int uj=0; uj<N_format; uj++)
		{
			if (uj != 0)
				tmpstream << ":";

			if ( ((int)uj == GT_idx) && (include_genotype[ui] == false) )
			{
				tmpstr = ".";
				for (int uk=1; uk<ploidy[ui]; uk++)
					tmpstr = tmpstr + "/.";
				format_miss = false;
			}
			else if ((int)uj == GT_idx)
			{
				read_indv_generic_entry(ui, uj, tmpstr);
				format_miss = false;
			}
			else
			{
				read_indv_generic_entry(ui, uj, tmpstr);
				format_miss = (tmpstr == ".");
			}
			indv_miss = indv_miss && format_miss;
			tmpstream << tmpstr;

			if (!format_miss)
			{
				tmpstr = tmpstream.str();
				outstream << tmpstr;
				tmpstream.str("");
			}
		}
		tmpstream.str("");
	}
	out = outstream.str();
}

// Output BCF entry to output stream in VCF format
void bcf_entry::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	if (fully_parsed == false)
		parse_full_entry();

	out << get_CHROM() << '\t' << POS << '\t' << get_ID() << '\t' << REF << '\t' << get_ALT();
	out << '\t' << header::double2str(QUAL);
	out << '\t' << get_FILTER();
	out << '\t' << get_INFO(INFO_to_keep, keep_all_INFO);

	if (FORMAT.size() > 0)
	{
		string indv_entries;
		out << '\t' << get_FORMAT();

		read_all_entries(indv_entries);
		out << indv_entries;
	}
	out << '\n';	// endl flushes the buffer, which is slow. This (should be) quicker.
}

// Output BCF entry to output stream in BCF format
void bcf_entry::print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	vector<char> out_vector, tmp_vector;
	vector<pair< string, string > > tmp_info;
	int index;

	out_vector.resize(INFO_pos);
	memcpy(&out_vector[0], &line[0], INFO_pos);

	if (keep_all_INFO)
	{
		unsigned int curr_size = out_vector.size();
		out_vector.resize(curr_size + (FORMAT_pos - INFO_pos) );
		memcpy(&out_vector[curr_size], &line[INFO_pos], (FORMAT_pos - INFO_pos));
	}
	else
	{
		int map_type, number;
		tmp_info = get_INFO_vector(INFO_to_keep, keep_all_INFO);
		N_INFO_removed = INFO.size()-tmp_info.size();

		get_n_allele_info(tmp_vector);
		memcpy(&out_vector[6*sizeof(int32_t)], &tmp_vector[0], sizeof(char));

		for(unsigned int ui=0; ui<tmp_info.size(); ui++)
		{
			tmp_vector.resize(0);
			index = entry_header.INFO_reverse_map[ tmp_info[ui].first ];
			make_typed_int(tmp_vector, index, true);
			out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());

			tmp_vector.resize(0);
			map_type = entry_header.INFO_map[ index ].Type;
			number = entry_header.INFO_map[ index ].N_entries;

			if (map_type == Integer)
				make_typed_int_vector(tmp_vector, tmp_info[ui].second, number );
			else if (map_type == Float)
				make_typed_float_vector(tmp_vector, tmp_info[ui].second, number );
			else if ( (map_type == Character) or (map_type == String) )
				make_typed_string(tmp_vector, tmp_info[ui].second, true );
			else if (map_type == Flag)
				make_typed_int(tmp_vector, 1, true );
			else
				LOG.error("Invalid type in INFO definition", 0);

			out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());
		}
	}

	uint32_t l_shared = (uint32_t)out_vector.size() - (uint32_t)(2*sizeof(uint32_t));
	memcpy(&out_vector[0], &l_shared, sizeof(uint32_t));

	unsigned int size, l_pos, type;
	for (unsigned int ui=0; ui<N_format; ui++)
	{
		type = FORMAT_types[ui];
		size = FORMAT_sizes[ui];

		tmp_vector.resize(0);
		make_typed_int(tmp_vector, FORMAT_keys[ui], true);
		out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end() );

		tmp_vector.resize(0);
		make_type_size(tmp_vector, type, size);
		out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end() );

		for (unsigned int uj=0; uj<N_indv; uj++)
		{
			if (include_indv[uj] == false)
				continue;

			tmp_vector.resize(FORMAT_skip[ui]);
			l_pos = FORMAT_positions[ui] + FORMAT_skip[ui]*uj;

			if ( ((int)uj == GT_idx) and (include_genotype[uj] == false) )
				for (unsigned int ploidy = 0; ploidy < FORMAT_skip[ui]; ploidy++)
					tmp_vector[ploidy] = (int8_t)0x81;
			else
				memcpy(&tmp_vector[0], &line[l_pos], FORMAT_skip[ui]);

			out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end() );
		}
	}
	uint32_t l_indv = (uint32_t)out_vector.size() - l_shared - (uint32_t)(2*sizeof(uint32_t));
	memcpy(&out_vector[sizeof(uint32_t)], &l_indv, sizeof(uint32_t));

	bgzf_write(out, &out_vector[0], out_vector.size());
}
