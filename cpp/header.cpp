/*
 * header.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: amarcketta
 */

#include "header.h"

header::header()
{
	has_contigs = false; has_file_format = false; has_genotypes = false;
	has_header = false; has_idx = false;
	contig_index = 0; N_indv = 0;
}

void header::parse_meta(const string &line, unsigned int &line_index)
{
	lines.push_back(line);
	if (line.find("##fileformat=")!=string::npos)
	{
		has_file_format = true;
		size_t found = line.find_first_of("=");
		string version = line.substr(found+1);
		if ((version != "VCFv4.0") && (version != "VCFv4.1") && (version != "VCFv4.2"))
			LOG.error("VCF version must be v4.0, v4.1 or v4.2:\nYou are using version " + version);
	}
	else if (line.find("##INFO=")!=string::npos)
	{	// Found an INFO descriptor
		line_index += add_INFO_descriptor(line.substr(8, line.size()-8), line_index);
	}
	else if (line.find("##FILTER=")!=string::npos)
	{	// Found a FILTER descriptor
		line_index += add_FILTER_descriptor(line.substr(10, line.size()-8), line_index);
	}
	else if (line.find("##FORMAT=")!=string::npos)
	{	// Found a genotype filter descriptor
		line_index += add_FORMAT_descriptor(line.substr(10, line.size()-8), line_index);
	}
	else if (line.find("##contig=")!=string::npos)
	{	// Found a contig descriptor
		add_CONTIG_descriptor(line.substr(10, line.size()-8), contig_index);
		contig_index++;
		has_contigs = true;
	}
	else
	{
		Field_description I;
		size_t found = line.find_first_of("=");
		I.Field = line.substr(0,found);
		I.Other = line.substr(found+1);
		parsed_lines.push_back(I);
	}
}

void header::parse_header(const string &line)
{
	// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	(FORMAT	NA00001 NA00002 ... )
	if (has_header == true)
		LOG.warning("Multiple Header lines.");

	has_header = true;
	istringstream header(line);
	int count = 0;
	string tmp_str;
	unsigned int N_header_indv = 0;
	has_genotypes = false;
	while (!header.eof())
	{
		getline(header, tmp_str, '\t');
		switch (count)
		{
			case 0: if (tmp_str != "#CHROM") LOG.warning("First Header entry should be #CHROM: " + tmp_str); break;
			case 1: if (tmp_str != "POS") LOG.warning("Second Header entry should be POS: " + tmp_str); break;
			case 2: if (tmp_str != "ID") LOG.warning("Third Header entry should be ID: " + tmp_str); break;
			case 3: if (tmp_str != "REF") LOG.warning("Fourth Header entry should be REF: " + tmp_str); break;
			case 4: if (tmp_str != "ALT") LOG.warning("Fifth Header entry should be ALT: " + tmp_str); break;
			case 5: if (tmp_str != "QUAL") LOG.warning("Sixth Header entry should be QUAL: " + tmp_str); break;
			case 6: if (tmp_str != "FILTER") LOG.warning("Seventh Header entry should be FILTER: " + tmp_str); break;
			case 7: if (tmp_str != "INFO") LOG.warning("Eighth Header entry should be INFO: " + tmp_str); break;
			case 8:
				if (tmp_str != "FORMAT")
					LOG.warning("Ninth Header entry should be FORMAT: " + tmp_str);
				else
					has_genotypes = true;
				break;
			default:
			{
				if (count <= 8)
					LOG.error("Incorrectly formatted header.");
				indv.push_back(tmp_str);
				N_header_indv++;
			}
			break;
		}
		count++;
	}
	N_indv = N_header_indv;

	if ((has_genotypes == true ) && (N_indv == 0))
		LOG.warning("FORMAT field without genotypes?");
}

int header::add_INFO_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end);

	Field_description I;
	I.Field = "INFO";
	vector<string> tokens;
	tokenize(details, ',', tokens);

	if (tokens.size() < 4)
		LOG.error("Expected at least 4 parts in INFO definition: " + in);

	vector<string> entry;
	for (unsigned int ui=0; ui<tokens.size(); ui++)
	{
		tokenize(tokens[ui], '=', entry);

		if (entry[0] == "ID") I.ID = entry[1];
		else if (entry[0] == "Number")
		{
			if ((entry[1] == "A") || (entry[1] == "G"))
				I.N_entries = -1;
			else
				I.N_entries =  str2int(entry[1]);
			I.N_entries_str = entry[1];
		}
		else if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") { I.Type_str = "Integer"; I.Type = Integer; }
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) {I.Type_str = "Float"; I.Type = Float;}
			else if (entry[1] == "Character") {I.Type_str = "Character"; I.Type = Character;}
			else if (entry[1] == "String") {I.Type_str = "String"; I.Type = String;}
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				I.Type_str = "Flag";
				if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
			}
			else LOG.error("Unknown Type in INFO meta-information: " + in);
		}
		else if (entry[0] == "Description") I.Description = entry[1];
		else if (entry[0] == "Source") I.Source = entry[1];
		else if (entry[0] == "Version") I.Version = entry[1];
		else if (entry[0] == "IDX")
		{
			has_idx = true;
			I.idx = str2int(entry[1]);
		}
		else
		{
			if (I.Other != "")
				I.Other += ",";
			I.Other += tokens[ui];
		}
	}

	if (I.ID == "") LOG.error("ID required in INFO field description: " + in);
	if (I.N_entries_str == "") LOG.error("Number required in INFO field description: " + in);
	if (I.Type_str == "") LOG.error("Type required in INFO field description: " + in);
	if (I.Description == "") LOG.error("Description required in INFO field description: " + in);
	if ((I.idx == -1) && (has_idx)) LOG.error("Missing index in INFO field description: " + in);
	parsed_lines.push_back(I);

	if (I.idx != -1)
	{
		INFO_map[I.idx] = I;
		INFO_reverse_map[I.ID] = I.idx;
		return 0;
	}
	else
	{
		if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
		{
			INFO_map[ FORMAT_reverse_map[ I.ID ] ] = I;
			INFO_reverse_map[I.ID] = FORMAT_reverse_map[I.ID];
			return 0;
		}
		else if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
		{
			INFO_map[ FILTER_reverse_map[ I.ID ] ] = I;
			INFO_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
			return 0;
		}
		else if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
			return 0;
		else
		{
			INFO_map[index] = I;
			INFO_reverse_map[I.ID] = index;
			return 1;
		}
	}
}

int header::add_FORMAT_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end);

	vector<string> tokens;
	tokenize(details, ',', tokens);
	Field_description I;
	I.Field = "FORMAT";

	if (tokens.size() < 4)
		LOG.error("Expected at least 4 parts in FORMAT definition: " + in);

	vector<string> entry;
	for (unsigned int ui=0; ui<tokens.size(); ui++)
	{
		tokenize(tokens[ui], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else if (entry[0] == "Number")
		{
			if ((entry[1] == "A") || (entry[1] == "G"))
				I.N_entries = -1;
			else
				I.N_entries = str2int(entry[1]);
			I.N_entries_str = entry[1];
		}
		else if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") {I.Type_str = "Integer"; I.Type = Integer;}
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) {I.Type_str = "Float"; I.Type = Float;}
			else if (entry[1] == "Character") {I.Type_str = "Character"; I.Type = Character;}
			else if (entry[1] == "String") {I.Type_str = "String"; I.Type = String;}
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				I.Type_str = "Flag";
				if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
			}
			else LOG.error("Unknown Type in FORMAT meta-information: " + in);
		}
		else if (entry[0] == "Description")	I.Description = entry[1];
		else if (entry[0] == "IDX")
		{
			has_idx = true;
			I.idx = str2int(entry[1]);
		}
		else
		{
			if (I.Other != "")
				I.Other += ",";
			I.Other += tokens[ui];
		}
	}

	if (I.ID == "") LOG.error("ID required in FORMAT field description: " + in);
	if (I.N_entries_str == "") LOG.error("Number required in FORMAT field description: " + in);
	if (I.Type_str == "") LOG.error("Type required in FORMAT field description: " + in);
	if (I.Description == "") LOG.error("Description required in FORMAT field description: " + in);
	if ((I.idx == -1) && (has_idx)) LOG.error("Missing index in FORMAT field description: " + in);
	parsed_lines.push_back(I);

	if (I.idx != -1)
	{
		FORMAT_map[I.idx] = I;
		FORMAT_reverse_map[I.ID] = I.idx;
		return 0;
	}
	else
	{
		if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
		{
			FORMAT_map[ FILTER_reverse_map[ I.ID ] ] = I;
			FORMAT_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
			return 0;
		}
		else if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
		{
			FORMAT_map[ INFO_reverse_map[ I.ID ] ] = I;
			FORMAT_reverse_map[I.ID] = INFO_reverse_map[ I.ID ];
			return 0;
		}
		else if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
			return 0;
		else
		{
			FORMAT_map[ index ] = I;
			FORMAT_reverse_map[I.ID] = index;
			return 1;
		}
	}
}

void header::add_CONTIG_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end);

	vector<string> tokens;
	tokenize(details, ',', tokens);
	Field_description I;
	I.Field = "contig";
	bool id_found = false;
	vector<string> entry;

	for (unsigned int ui=0; ui<tokens.size(); ui++)
	{
		tokenize(tokens[ui], '=', entry);
		if (entry[0] == "ID")
		{
			I.ID = entry[1];
			id_found = true;
		}
		else if (entry[0] == "length") I.Length = entry[1];
		else if (entry[0] == "assembly") I.Assembly = entry[1];
		else
		{
			if (I.Other != "")
				I.Other += ",";
			I.Other += tokens[ui];
		}
	}
	if (id_found == false)
		LOG.warning("CONTIG declaration found without ID: "+ in + "\n");

	parsed_lines.push_back(I);
	CONTIG_map[index] = I;
	CONTIG_reverse_map[I.ID] = index;
}

int header::add_FILTER_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end);
	vector<string> tokens;
	tokenize(details, ',', tokens);
	if (tokens.size() < 2)
		LOG.error("Expected at least 2 parts in FILTER definition: " + in);

	string Description;
	Field_description I;
	I.Field = "FILTER";
	vector<string> entry;

	for (unsigned int ui=0; ui<tokens.size(); ui++)
	{
		tokenize(tokens[ui], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else if (entry[0] == "Description")
		{
			I.Description = entry[1];
		}
		else if (entry[0] == "IDX")
		{
			has_idx = true;
			I.idx = str2int(entry[1]);
		}
		else
		{
			if (I.Other != "")
				I.Other += ",";
			I.Other += tokens[ui];
		}
	}

	if (I.ID == "") LOG.error("ID required in FILTER field description: " + in);
	if (I.Description == "") LOG.error("Description required in FILTER field description: " + in);
	if ((I.idx == -1) && (has_idx)) LOG.error("Missing index in FILTER field description: " + in);
	if (I.ID != "PASS") parsed_lines.push_back(I);

	if (I.idx != -1)
	{
		FILTER_map[I.idx] = I;
		FILTER_reverse_map[I.ID] = I.idx;
		return 0;
	}
	else
	{
		if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
		{
			FILTER_map[ INFO_reverse_map[ I.ID ] ] = I;
			FILTER_reverse_map[I.ID] = INFO_reverse_map[ I.ID ];
			return 0;
		}
		else if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
		{
			FILTER_map[ FORMAT_reverse_map[ I.ID ] ] = I;
			FILTER_reverse_map[I.ID] = FORMAT_reverse_map[ I.ID ];
			return 0;
		}
		else if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
			return 0;
		else
		{
			FILTER_map[index] = I;
			FILTER_reverse_map[I.ID] = index;
			return 1;
		}
	}
}

void header::reprint()
{
	lines.resize(0);
	Field_description I;
	for (unsigned int ui=0; ui<parsed_lines.size(); ui++)
	{
		I = parsed_lines[ui];
		stringstream new_line;

		if (I.Field[0] == '#')
			new_line << I.Field << "=" << I.Other;
		else
		{
			new_line << "##" << I.Field << "=<";
			if (I.ID != "") new_line << "ID=" << I.ID;
			if (I.N_entries_str != "") new_line << ",Number=" << I.N_entries_str;
			if (I.Type_str != "") new_line << ",Type=" << I.Type_str;
			if (I.Description != "") new_line << ",Description=" << I.Description;
			if (I.Source != "") new_line << ",Source=" << I.Source;
			if (I.Version != "") new_line << ",Version=" << I.Version;
			if (I.Length != "") new_line << ",Length=" << I.Length;
			if (I.Assembly != "") new_line << ",Assembly=" << I.Assembly;
			if (I.Other != "") new_line << "," << I.Other;
			new_line << ">";
		}
		lines.push_back(new_line.str());
	}
}

void header::reparse()
{
	unsigned int index = 0;
	has_idx = false; contig_index = 0;
	vector<string> old_lines(lines.size(),"");
	copy(lines.begin(), lines.end(), old_lines.begin());
	lines.resize(0);

	INFO_map.clear(); INFO_reverse_map.clear();
	FILTER_map.clear(); FILTER_reverse_map.clear();
	FORMAT_map.clear(); FORMAT_reverse_map.clear();
	CONTIG_map.clear(); CONTIG_reverse_map.clear();

	index += add_FILTER_descriptor("ID=PASS,Description=PASS", index);
	for (unsigned int ui=0; ui<old_lines.size(); ui++)
		parse_meta(old_lines[ui],index);
}

void header::tokenize(const string &in, char token, vector<string> &out)
{
	out.resize(0);
	istringstream ss(in);
	string tmp;
	while( getline(ss, tmp, token) )
	{
		out.push_back(tmp);
	}
}

void header::split(const string &text, char sep, vector<string> &tokens)
{
	int start = 0, end = 0, idx = 0, max = tokens.size();
	while ((end = text.find(sep, start)) != string::npos)
	{
		if (idx < max)
			tokens[idx] = text.substr(start, end - start);
		else
			tokens.push_back(text.substr(start, end - start));
		start = end + 1;
		idx++;
	}
	if (idx < max)
		tokens[idx] = text.substr(start);
	else
		tokens.push_back(text.substr(start));
}

string header::int2str(const int in, const int missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}

int header::str2int(const string &in, const int missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atoi(in.c_str());
}

double header::str2double(const string &in, const double missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atof(in.c_str());
}

string header::double2str(const double in, const double missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}
