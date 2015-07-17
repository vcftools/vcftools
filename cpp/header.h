/*
 * header.h
 *
 *  Created on: Apr 29, 2013
 *      Author: amarcketta
 */

#ifndef HEADER_H_
#define HEADER_H_

#include <cstring>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "output_log.h"

using namespace std;
extern output_log LOG;

enum Type_enum {Integer=0, Float=1, Character=2, String=3, Flag=4};

class Field_description
{
public:
	string Field;
	string ID;
	int idx;
	int N_entries;
	string N_entries_str;
	string Type_str;
	Type_enum Type;
	string Description;
	string Length;
	string Assembly;
	string Source;
	string Version;
	string Other;

	Field_description() : Field(""), ID(""), idx(-1), N_entries(0), N_entries_str(""), Type_str(""), Type(Integer), Description(""), Length(""), Assembly(""), Source(""), Version(""), Other("") {};
	~Field_description() {};
};

class header
{
public:
	unsigned int contig_index;
	bool has_contigs;
	bool has_genotypes;
	bool has_header;
	bool has_file_format;
	bool has_idx;
	vector<string> indv;
	vector<string> lines;
	vector<Field_description> parsed_lines;
	unsigned int N_indv;

	map<int, Field_description> INFO_map;
	map<int, Field_description> FILTER_map;
	map<int, Field_description> FORMAT_map;
	map<int, Field_description> CONTIG_map;
	map<string, int> CONTIG_reverse_map;
	map<string, int> FILTER_reverse_map;
	map<string, int> INFO_reverse_map;
	map<string, int> FORMAT_reverse_map;

	header();
	~header() {};

	void reprint();
	void reparse();
	void parse_meta(const string &line, unsigned int &line_index);
	void parse_header(const string &line);

	int add_INFO_descriptor(const string &in, int index);
	int add_FILTER_descriptor(const string &in, int index);
	int add_FORMAT_descriptor(const string &in, int index);
	void add_CONTIG_descriptor(const string &in, int index);

	static void tokenize(const string &in, char token, vector<string> &out);
	static void split(const string &in, char token, vector<string> &out);
	static int str2int(const string &in, const int missing_value=-1);
	static string int2str(const int in, const int missing_value=-1);
	static double str2double(const string &in, const double missing_value=-1.0);
	static string double2str(const double in, const double missing_value=-1.0);
};

#endif /* HEADER_H_ */
