/*
 * entry.h
 *
 *  Created on: Dec 12, 2012
 *      Author: amarcketta
 */
#ifndef ENTRY_H_
#define ENTRY_H_

#include <set>
#include <deque>
#include <vector>
#include <cassert>

#include "header.h"
#include "bgzf.h"
#include "output_log.h"
#include "parameters.h"

using namespace std;
extern output_log LOG;

class entry
{
public:
	virtual ~entry() {};
	unsigned int N_indv;
	bool passed_filters;
	header entry_header;
	vector<bool> include_indv;
	vector<bool> include_genotype;

	virtual void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false) = 0;
	virtual void parse_full_entry(bool parse_FORMAT=true) = 0;
	virtual void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false) = 0;
	virtual void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false) = 0;

	virtual void reset(const vector<char> &data_line) = 0;
	int apply_filters(const parameters &params);

	void filter_sites(const set<string> &snps_to_keep, const string &snps_to_keep_file, const string &snps_to_exclude_file, bool keep_then_exclude = false);
	void filter_sites_to_keep(const set<string> &snps_to_keep, const string &snps_to_keep_file);
	void filter_sites_to_exclude(const string &snps_to_exclude_file);
	void filter_sites_by_position(const string &chr, int start_pos, int end_pos);
	void filter_sites_by_positions(const string &positions_file, const string &exclude_positions_file);
	void filter_sites_by_overlap_positions(const string &positions_overlap_file, const string &exclude_positions_overlap_file);
	void filter_sites_by_chromosome(const set<string> &chrs_to_keep, const set<string> &chrs_to_exclude);
	void filter_sites_by_quality(double min_quality);
	void filter_sites_by_mean_depth(double min_mean_depth, double max_mean_depth);
	void filter_sites_by_frequency_and_call_rate(double min_maf, double max_maf, double min_non_ref_af, double max_non_ref_af, double min_non_ref_af_any, double max_non_ref_af_any, double min_site_call_rate);
	void filter_sites_by_allele_type(bool keep_only_indels, bool remove_indels);
	void filter_sites_by_allele_count(double min_mac, double max_mac, double min_non_ref_ac, double max_non_ref_ac, double min_non_ref_ac_any, double max_non_ref_ac_any, double max_missing_call_count);
	void filter_sites_by_number_of_alleles(int min_alleles, int max_alleles);
	void filter_sites_by_HWE_pvalue(double min_HWE_pvalue);
	void filter_sites_by_BED_file(const string &bed_file, bool BED_exclude = false);
	void filter_sites_by_mask(const string &mask_file, bool invert_mask = false, int min_kept_mask_value=0);
	void filter_sites_by_filter_status(const set<string> &filter_flags_to_remove, const set<string> &filter_flags_to_keep, bool remove_all = false);
	void filter_sites_by_phase();
	void filter_sites_by_thinning(int min_SNP_distance);
	void filter_sites_by_INFO(const set<string> &flags_to_remove, const set<string> &flags_to_keep);

	void filter_genotypes_by_quality_value(double min_genotype_quality);
	void filter_genotypes_by_depth_range(int min_depth, int max_depth);
	void filter_genotypes_by_filter_flag(const set<string> &filter_flags_to_remove, bool remove_all = false);

	string get_CHROM() const;
	void get_CHROM(string &out) const;
	int get_POS() const;
	string get_ID() const;
	string get_REF() const;
	string get_ALT() const;
	string get_ALT_allele(int allele_num) const;
	void get_allele(int allele_num, string &out) const;
	string get_allele(int allele_num) const;
	void get_alleles_vector(vector<string> &out) const;
	string get_FILTER() const;
	void get_FILTER_vector(vector<string> &out) const;
	double get_QUAL() const;
	string get_INFO(const set<string> &INFO_to_keep, bool keep_all_INFO=false) const;
	string get_INFO_value(const string &key) const;
	vector<string> get_INFO_values(const string &key) const;
	string get_FORMAT() const;

	void get_indv_GENOTYPE_ids(unsigned int indv, pair<int, int> &out) const;
	void get_indv_GENOTYPE_strings(unsigned int indv, pair<string, string> &out) const;
	char get_indv_PHASE(unsigned int indv) const;
	double get_indv_GQUALITY(unsigned int indv) const;
	int get_indv_DEPTH(unsigned int indv) const;
	void get_indv_GFILTER(unsigned int indv, string &out) const;
	void get_indv_GFILTER_vector(unsigned int indv, vector<string> &out) const;
	int get_indv_ploidy(unsigned int indv) const;

	bool is_SNP() const;
	bool is_biallelic_SNP() const;
	bool is_diploid() const;
	virtual void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out) = 0;
	bool FORMAT_id_exists(const string &FORMAT_id);

	void get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out) const;
	void get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out, const vector<bool> &include_indv, const vector<bool> &include_genotype) const;
	void get_genotype_counts(unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const;
	void get_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const;
	void get_multiple_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, vector<unsigned int> &out_N_hom, vector<unsigned int> &out_N_het) const;
	unsigned int get_N_alleles() const;
	unsigned int get_N_chr() const;

	void get_POS_binary(vector<char> &out) const;
	void get_ID_binary(vector<char> &out);
	void get_rlen(vector<char> &out) const;
	void get_QUAL_binary(vector<char> &out) const;
	void get_n_allele_info(vector<char> &out) const;
	void get_n_fmt_sample(vector<char> &out) const;
	void get_ALLELES_binary(vector<char> &out);
	vector<pair<string, string> > get_INFO_vector(const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void get_FORMAT_binary(vector<char> &out) const;

	string get_typed_string( unsigned int * line_position, const vector<char>& line );
	void get_type(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size);
	vector<int> get_int_vector(unsigned int * line_position, const vector<char>& line);
	int get_typed_int(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size);
	void get_number(uint32_t &out, unsigned int * line_position, const vector<char>& line);

	void make_typed_string(vector<char> &out, const string &in, bool typed);
	void make_typed_int(vector<char> &out, const int &in, bool typed);
	void make_int(vector<char> &out, const int &in, int type);
	void make_typed_int_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_int_vector(vector<char> &out, const string &in, int number = -1);
	void make_typed_int_vector(vector<char> &out, const vector<int> &in);
	void make_typed_float_vector(vector<char> &out, const string &in, int number = -1);
	void make_typed_float_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_string_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_GT_vector(vector<char> &out, vector<string> &in);
	void make_type_size(vector<char> &out, const unsigned int &type, const unsigned int &size);
	void encode_genotype(vector<char> &out, string &in, int exp_size);

	void skip_section(unsigned int *line_position, const vector<char> &line);
	bool check_missing(unsigned int line_position, const unsigned int type, const vector<char> &line);
	bool check_end(unsigned int line_position, const unsigned int type, const vector<char> &line);
	void add_ALT_allele(const string &in);

	virtual void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false) = 0;
	virtual void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO=false) = 0;

	virtual void filter_genotypes_by_depth(int min_depth, int max_depth) = 0;
	virtual void filter_genotypes_by_quality(double min_genotype_quality) = 0;
	virtual void filter_genotypes_by_filter_status(const set<string> &filter_flags_to_remove, bool remove_all = false) = 0;

	static void SNPHWE(int obs_hets, int obs_hom1, int obs_hom2, double &p_hwe, double &p_lo, double &p_hi);

	static set<string> local_snps_to_keep;
	static set<string> snps_to_exclude;
	static vector< set<int > > keep_positions;
	static vector< set<int > > exclude_positions;
	static map<string,int> chr_to_idx;
	static vector< deque<pair<int,int> > > lims;
	static ifstream mask;
	static string mask_chr;
	static string mask_line;
	static int mask_pos;
	static string thin_chrom;
	static int thin_pos;

protected:
	istringstream data_stream;
	vector<char> line;

	bool basic_parsed;
	bool fully_parsed;
	bool parsed_ALT;
	bool parsed_FILTER;
	bool parsed_INFO;
	bool parsed_FORMAT;
	bool parsed_FORMAT_binary;

	string CHROM;
	int POS;
	string ID;
	string REF;
	vector<string> ALT;
	double QUAL;
	vector<string> FILTER;
	vector<pair<string, string> > INFO;
	vector<string> FORMAT;
	vector<char> FORMAT_binary;
	int N_INFO_removed;
	int N_FORMAT_removed;

	vector< pair<int,int> > GENOTYPE;
	vector<int> ploidy;
	vector<char>  PHASE;
	vector<double> GQUALITY;
	vector<int>   DEPTH;
	vector< vector<string> > GFILTER;

	vector<bool> parsed_GT;
	vector<bool> parsed_GQ;
	vector<bool> parsed_DP;
	vector<bool> parsed_FT;

	map<string, unsigned int> FORMAT_to_idx;
	int GT_idx;
	int GQ_idx;
	int DP_idx;
	int FT_idx;

	void set_indv_DEPTH(unsigned int indv, int in);
	vector<unsigned int> FORMAT_positions, FORMAT_types, FORMAT_sizes, FORMAT_skip, FORMAT_keys;
};

#endif /* ENTRY_H_ */
