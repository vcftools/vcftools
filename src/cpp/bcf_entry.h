/*
 * bcf_entry.h
 *
 *  Created on: Sep 20, 2012
 *      Author: Anthony Marcketta
 *      ($Revision: 1 $)
 */

#include "output_log.h"
#include "entry.h"
#include "header.h"

extern output_log LOG;

class bcf_entry : public entry {
public:
	bcf_entry(header &header_obj, vector<bool> &include_individual);
	~bcf_entry();

	unsigned int N_info;
	unsigned int N_format;
	unsigned int N_allele;
	unsigned int L_shared;
	unsigned int L_indiv;
	unsigned int line_pos;
	ostringstream outstream;
	ostringstream tmpstream;

	void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false);
	void parse_full_entry(bool parse_FORMAT=true);
	void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false);
	void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false);

	void set_ALT(const int n_allele);
	void set_ALT(const string &in);
	void set_FILTER();
	void set_FORMAT();
	void set_INFO();
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const unsigned int &pos, const unsigned int &size);
	void set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in);
	void set_indv_GQUALITY(unsigned int indv, const vector<char> &in);
	void set_indv_GQUALITY(unsigned int indv, const float &in);
	void set_indv_GFILTER(unsigned int indv, const string &in);
	void set_indv_GFILTER(unsigned int indv, const vector<char> &in);
	void set_indv_PHASE(unsigned int indv, char in);
	void set_indv_GENOTYPE_alleles(unsigned int indv, const pair<int, int> &in);
	void reset(const vector<char> &data_line);

	void add_FORMAT_entry(const string &in, const unsigned int &fmt_key, const unsigned int &pos, const unsigned int &line_pos, const unsigned int &type, const unsigned int &size);
	void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out);
	void read_indv_generic_entry(unsigned int indv, const int &idx, string &out);
	void read_all_entries(string &out);

	void filter_genotypes_by_quality(double min_genotype_quality);
	void filter_genotypes_by_depth(int min_depth, int max_depth);
	void filter_genotypes_by_filter_status(const set<string> &filter_flags_to_remove, bool remove_all = false);

	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);

private:
	vector<int> FILTER_str;
	unsigned int INFO_pos, FILTER_pos, ALT_pos, FORMAT_pos;
};
