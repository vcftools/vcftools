/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */

// Class for reading in, checking and storing user parameters
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include <unistd.h>

#include "output_log.h"

extern output_log LOG;

using namespace std;

const string VCFTOOLS_VERSION="v0.1.13";
static const uint8_t bgzf_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0"; //just compare the first 16 chars? though
static const uint8_t gzip_magic[2] = {0x1f,0x8b};

class parameters
{
public:
	bool stream_in;
	bool bcf_format;
	bool BED_exclude;
	string BED_file;
	set<string> chrs_to_exclude;
	set<string> chrs_to_keep;
	string chrom_map_file;
	string contigs_file;
	bool derived;
	bool diff_discordance_matrix;
	string diff_file;
	bool diff_file_bcf;
	bool diff_file_compressed;
	bool diff_indv;
	bool diff_indv_discordance;
	string diff_indv_map_file;
	bool diff_site;
	bool diff_site_discordance;
	bool diff_switch_error;
	int end_pos;
	string exclude_positions_file;
	string exclude_positions_overlap_file;
	string FORMAT_id_to_extract;
	set<string> geno_filter_flags_to_exclude;
	string geno_rsq_position_list;
	string hap_rsq_position_list;
	string hapcount_BED;
	vector<string> weir_fst_populations;
	int fst_window_size;
	int fst_window_step;
	vector<string> indv_exclude_files;
	vector<string> indv_keep_files;
	set<string> indv_to_exclude;
	set<string> indv_to_keep;
	vector<string> INFO_to_extract;
	bool invert_mask;
	bool keep_only_indels;
	int ld_bp_window_size;
	int ld_snp_window_size;
	int ld_bp_window_min;
	int ld_snp_window_min;
	int min_mac;
	double min_maf;
	string mask_file;
	int max_alleles;
	int max_genotype_depth;
	int max_mac;
	double max_maf;
	double max_mean_depth;
	int max_missing_call_count;
	int max_non_ref_ac;
	double max_non_ref_af;
	int max_non_ref_ac_any;
	double max_non_ref_af_any;
	int max_N_indv;
	string mendel_ped_file;
	int min_alleles;
	int min_genotype_depth;
	double min_genotype_quality;
	double min_HWE_pvalue;
	int min_interSNP_distance;
	int min_kept_mask_value;
	double min_mean_depth;
	int min_non_ref_ac;
	double min_non_ref_af;
	int min_non_ref_ac_any;
	double min_non_ref_af_any;
	double min_quality;
	double min_r2;
	double min_site_call_rate;
	int num_outputs;
	bool output_012_matrix;
	bool output_as_IMPUTE;
	bool output_as_ldhat_phased;
	bool output_as_ldhat_unphased;
	bool output_BEAGLE_genotype_likelihoods_GL;
	bool output_BEAGLE_genotype_likelihoods_PL;
	bool output_counts;
	bool output_filter_summary;
	bool output_freq;
	bool output_geno_depth;
	bool output_geno_chisq;
	bool output_geno_rsq;
	bool output_hap_rsq;
	bool output_het;
	bool output_HWE;
	bool output_indel_hist;
	bool output_indv_burden;
	bool output_indv_depth;
	bool output_indv_freq_burden;
	bool output_indv_freq_burden2;
	bool output_indv_missingness;
	bool output_interchromosomal_hap_rsq;
	bool output_interchromosomal_geno_rsq;
	bool output_kept_sites;
	bool output_LROH;
	int output_N_PCA_SNP_loadings;
	bool output_PCA;
	string output_prefix;
	bool output_relatedness_Yang;
	bool output_relatedness_Manichaikul;
	bool output_removed_sites;
	bool output_singletons;
	bool output_site_depth;
	bool output_site_mean_depth;
	bool output_site_missingness;
	bool output_site_pi;
	bool output_site_quality;
	int output_SNP_density_bin_size;
	int output_Tajima_D_bin_size;
	int output_TsTv_bin_size;
	bool output_TsTv_by_count;
	bool output_TsTv_by_qual;
	bool output_TsTv_summary;
	bool phased_only;
	bool PCA_no_normalisation;
	int pi_window_size;
	int pi_window_step;
	bool plink_output;
	bool plink_tped_output;
	string positions_file;
	string positions_overlap_file;
	bool recode;
	bool recode_bcf;
	set<string> recode_INFO_to_keep;
	bool recode_all_INFO;
	bool remove_all_filtered_genotypes;
	bool remove_all_filtered_sites;
	bool remove_indels;
	set<string> site_filter_flags_to_exclude;
	set<string> site_filter_flags_to_keep;
	set<string> site_INFO_flags_to_keep;
	set<string> site_INFO_flags_to_remove;
	string snps_to_exclude_file;
	string snps_to_keep_file;
	set<string> snps_to_keep;
	int start_pos;
	bool stream_err;
	bool stream_out;
	bool suppress_allele_output;
	string temp_dir;
	string vcf_filename;
	bool vcf_format;
	bool vcf_compressed;

	parameters(int argc, char *argv[]);
	~parameters(){};

	void read_parameters();
	void print_help();
	void print_params();

private:
	void check_parameters();
	static void error(string err_msg, int code);

	vector<string> argv;

	string get_arg(unsigned int i);
};


#endif /* PARAMETERS_H_ */
