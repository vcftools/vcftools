/*
 * variant_file.h
 *
 *  Created on: Dec 12, 2012
 *      Author: amarcketta
 */

#ifndef VARIANT_FILE_H_
#define VARIANT_FILE_H_

#include <algorithm>
#include <bitset>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <set>
#include <sstream>
#include <map>
#include <numeric>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>
#include <zlib.h>

#include "parameters.h"
#include "entry.h"
#include "gamma.h"
#include "vcf_entry.h"
#include "bcf_entry.h"
#include "header.h"

#ifdef VCFTOOLS_PCA
	#include "dgeev.h"
#endif

extern output_log LOG;

using namespace std;

class variant_file
{
public:
	string filename;
	bool compressed;
	istream *file_in;
	ifstream file_tmp;
	unsigned int gzMAX_LINE_LEN;
	gzFile gzfile_in;

	header meta_data;
	vector<bool> include_indv;
	unsigned int N_entries;
	unsigned int N_kept_entries;

	int N_kept_individuals() const;
	int N_kept_sites() const;
	int N_total_sites() const;

	virtual void open() = 0;
	virtual void open_gz() = 0;
	virtual void close() = 0;
	virtual bool eof() = 0;

	virtual void get_entry(vector<char> &out) = 0;
	virtual entry* get_entry_object() = 0;
	void ByteSwap(unsigned char *b, int n) const;
	static inline bool is_big_endian() { long one= 1; return !(*((char *)(&one))); };

	void apply_filters(const parameters &params);
	void filter_individuals(const set<string> &indv_to_keep, const set<string> &indv_to_exclude, const vector<string> &indv_to_keep_filename, const vector<string> &indv_to_exclude_filename, bool keep_then_exclude=true);
	void filter_individuals_by_keep_list(const set<string> &indv_to_keep, const vector<string> &indv_to_keep_filenames);
	void filter_individuals_by_exclude_list(const set<string> &indv_to_exclude, const vector<string> &indv_to_exclude_filenames);
	void filter_individuals_randomly(int max_N_indv);

	void output_frequency(const parameters &params, bool output_counts=false);
	void output_individuals_by_mean_depth(const parameters &params);
	void output_site_depth(const parameters &params, bool output_mean=true);
	void output_genotype_depth(const parameters &params);
	void output_het(const parameters &params);
	void output_hwe(const parameters &params);
	void output_SNP_density(const parameters &params);
	void output_indv_missingness(const parameters &params);
	void output_indv_burden(const parameters &params);
	void output_indv_freq_burden(const parameters &params, int double_count_hom_alt=0);
	void output_site_missingness(const parameters &params);
	void output_haplotype_r2(const parameters &params);
	void output_genotype_r2(const parameters &params);
	void output_genotype_chisq(const parameters &params, double min_pval);
	void output_interchromosomal_genotype_r2(const parameters &params);
	void output_interchromosomal_haplotype_r2(const parameters & params);
	void output_haplotype_r2_of_SNP_list_vs_all_others(const parameters &params);
	void output_haplotype_count(const parameters &params);
	void output_genotype_r2_of_SNP_list_vs_all_others(const parameters &params);
	void output_singletons(const parameters &params);
	void output_TsTv(const parameters &params);
	void output_TsTv_by_count(const parameters &params);
	void output_TsTv_by_quality(const parameters &params);
	void output_TsTv_summary(const parameters &params);
	void output_per_site_nucleotide_diversity(const parameters &params);
	void output_windowed_nucleotide_diversity(const parameters &params);
	void output_Tajima_D(const parameters &params);
	void output_site_quality(const parameters &params);
	void output_FILTER_summary(const parameters &params);
	void output_kept_sites(const parameters &params);
	void output_removed_sites(const parameters &params);
	void output_LROH(const parameters &params);
	void output_indv_relatedness_Yang(const parameters &params);
	void output_indv_relatedness_Manichaikul(const parameters &params);
	void output_PCA(const parameters &params);
	void output_PCA_SNP_loadings(const parameters &params);
	void output_indel_hist(const parameters &params);

	void output_as_012_matrix(const parameters &params);
	void output_as_plink(const parameters &params);
	void output_as_plink_tped(const parameters &params);
	void output_BEAGLE_genotype_likelihoods(const parameters &params, int GL_or_PL=0);
	void output_as_IMPUTE(const parameters &params);
	void output_as_LDhat_phased(const parameters &params);
	void output_as_LDhat_unphased(const parameters &params);
	void output_FORMAT_information(const parameters &params);

	void output_weir_and_cockerham_fst(const parameters &params);
	void output_windowed_weir_and_cockerham_fst(const parameters &params);

	void output_sites_in_files(const parameters &params, variant_file &diff_vcf_file);
	void output_indv_in_files(const parameters &params, variant_file &diff_vcf_file);
	void output_discordance_by_site(const parameters &params, variant_file &diff_vcf_file);
	void output_discordance_matrix(const parameters &params, variant_file &diff_vcf_file);
	void output_discordance_by_indv(const parameters &params, variant_file &diff_vcf_file);
	void output_switch_error(const parameters &params, variant_file &diff_vcf_file);
	void output_INFO_for_each_site(const parameters &params);
	void output_mendel_inconsistencies(const parameters &params);
	void write_stats(const parameters &params);

	virtual void print(const parameters &params) = 0;
	virtual void print_bcf(const parameters &params) = 0;

	void calc_hap_r2(vector<pair<int,int> > &GT1, vector<pair<int,int> > &GT2, double &r2, double &D, double &Dprime, int &chr_count);
	void calc_geno_r2(vector<pair<int,int> > &GT1, vector<pair<int,int> > &GT2, double &r2, int &indv_count);
	void calc_r2_em(entry *e, entry *e2, double &r2, int &indv_count);
	void calc_geno_chisq(vector<pair<int,int> > &GT1, vector<pair<int,int> > &GT2, int &N0, int &N1, double &chisq, double &dof, double &pval, int &indv_count);
	void read_temp_site(ifstream &tmp_file, string &CHROM, int &POS, vector< pair<int,int> > &GTs);
	void read_big_temp_site(ifstream &tmp_file, string &CHROM, int &POS, int &alleles, vector< pair<int,int> > &GTs);
	void return_indv_union(variant_file &file2, map<string, pair< int, int> > &combined_individuals, const string &indv_ID_map_file="");

	void get_contigs(const std::string &contigs_file, vector<string> &contig_vector);
	virtual ~variant_file();
};

#endif /* VARIANT_FILE_H_ */
