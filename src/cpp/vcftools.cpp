/*
 * vcftools.cpp
 */

#include "vcftools.h"

output_log LOG;

int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);

	// The following turns off sync between C and C++ streams.
	// Apparently it's faster to turn sync off, and as I don't use C streams, it's okay to turn off.
	ios_base::sync_with_stdio(false);

	parameters params(argc, argv);
	params.print_help();
	params.read_parameters();

	LOG.open(params.stream_out, params.stream_err, params.output_prefix);

	LOG.printLOG("\nVCFtools - " + VCFTOOLS_VERSION + "\n");
	LOG.printLOG("(C) Adam Auton and Anthony Marcketta 2009\n\n");

	params.print_params();

	variant_file *vf;
	if (!params.bcf_format)
		vf = new vcf_file(params);
	else
		vf = new bcf_file(params);

	vf->apply_filters(params);

	LOG.printLOG("After filtering, kept " + output_log::int2str(vf->N_kept_individuals()) + " out of " + output_log::int2str(vf->meta_data.N_indv) + " Individuals\n");

	if (params.diff_file != "")
	{
		variant_file *variant_diff;
		if (params.diff_file_bcf)
			variant_diff = new bcf_file(params, true);
		else
			variant_diff = new vcf_file(params, true);

		variant_diff->apply_filters(params);
		if (params.diff_indv == true) vf->output_indv_in_files(params, *variant_diff);
		else if (params.diff_site_discordance == true) vf->output_discordance_by_site(params, *variant_diff);
		else if (params.diff_discordance_matrix == true) vf->output_discordance_matrix(params, *variant_diff);
		else if (params.diff_indv_discordance == true) vf->output_discordance_by_indv(params, *variant_diff);
		else if (params.diff_switch_error == true) vf->output_switch_error(params, *variant_diff);
		else if (params.diff_site == true)	vf->output_sites_in_files(params, *variant_diff);
		else LOG.warning("Diff file provided, but no additional option.\n");

		delete variant_diff;
	}
	if (params.num_outputs == 0) vf->write_stats(params);
	if (!params.INFO_to_extract.empty()) vf->output_INFO_for_each_site(params);
	if (params.FORMAT_id_to_extract != "") vf->output_FORMAT_information(params);
	if (params.output_indv_burden == true) vf->output_indv_burden(params);
	if (params.output_indv_depth == true) vf->output_individuals_by_mean_depth(params);
	if (params.output_indv_freq_burden == true) vf->output_indv_freq_burden(params);
	if (params.output_indv_freq_burden2 == true) vf->output_indv_freq_burden(params, 1);
	if (params.output_geno_depth == true) vf->output_genotype_depth(params);
	if (params.output_site_depth == true) vf->output_site_depth(params, false);
	if (params.output_site_mean_depth == true) vf->output_site_depth(params, true);
	if (params.output_freq == true) vf->output_frequency(params, false);
	if (params.output_counts == true) vf->output_frequency(params, true);
	if (params.plink_output == true) vf->output_as_plink(params);
	if (params.plink_tped_output == true) vf->output_as_plink_tped(params);
	if (params.output_HWE == true) vf->output_hwe(params);
	if (params.output_SNP_density_bin_size > 0) vf->output_SNP_density(params);
	if (params.output_indv_missingness == true) vf->output_indv_missingness(params);
	if (params.output_site_missingness == true) vf->output_site_missingness(params);
	if (params.output_geno_chisq == true) vf->output_genotype_chisq(params, -1.0);
	if (params.output_geno_rsq == true) vf->output_genotype_r2(params);
	if (params.output_interchromosomal_hap_rsq == true) vf->output_interchromosomal_haplotype_r2(params);
	if (params.output_interchromosomal_geno_rsq == true) vf->output_interchromosomal_genotype_r2(params);
	if (params.output_hap_rsq == true) vf->output_haplotype_r2(params);
	if (params.hap_rsq_position_list != "") vf->output_haplotype_r2_of_SNP_list_vs_all_others(params);
	if (params.geno_rsq_position_list != "") vf->output_genotype_r2_of_SNP_list_vs_all_others(params);
	if (params.output_het == true) vf->output_het(params);
	if (params.hapcount_BED != "") vf->output_haplotype_count(params);
	if (params.output_site_quality == true) vf->output_site_quality(params);
	if (params.output_012_matrix == true) vf->output_as_012_matrix(params);
	if (params.output_as_IMPUTE == true) vf->output_as_IMPUTE(params);
	if (params.output_BEAGLE_genotype_likelihoods_GL == true) vf->output_BEAGLE_genotype_likelihoods(params, 0);
	if (params.output_BEAGLE_genotype_likelihoods_PL == true) vf->output_BEAGLE_genotype_likelihoods(params, 1);
	if (params.output_as_ldhat_unphased == true) vf->output_as_LDhat_unphased(params);
	if (params.output_as_ldhat_phased == true) vf->output_as_LDhat_phased(params);
	if (params.output_singletons == true) vf->output_singletons(params);
	if (params.output_site_pi == true) vf->output_per_site_nucleotide_diversity(params);
	if (params.pi_window_size > 0) vf->output_windowed_nucleotide_diversity(params);
	if (params.output_Tajima_D_bin_size > 0) vf->output_Tajima_D(params);
	if (params.output_TsTv_bin_size > 0) vf->output_TsTv(params);
	if (params.output_TsTv_by_count) vf->output_TsTv_by_count(params);
	if (params.output_TsTv_by_qual) vf->output_TsTv_by_quality(params);
	if (params.output_TsTv_summary) vf->output_TsTv_summary(params);
	if (params.recode == true) vf->print(params);
	if (params.recode_bcf == true) vf->print_bcf(params);
	if (params.output_filter_summary == true) vf->output_FILTER_summary(params);
	if (params.output_kept_sites == true) vf->output_kept_sites(params);
	if (params.output_removed_sites == true) vf->output_removed_sites(params);

	if (params.output_LROH == true) vf->output_LROH(params);
	if (params.output_relatedness_Yang == true) vf->output_indv_relatedness_Yang(params);
	if (params.output_relatedness_Manichaikul == true) vf->output_indv_relatedness_Manichaikul(params);
	if (params.output_PCA == true) vf->output_PCA(params);
	if (params.output_N_PCA_SNP_loadings > 0) vf->output_PCA_SNP_loadings(params);
	if (params.mendel_ped_file != "") vf->output_mendel_inconsistencies(params);

	if (params.fst_window_size <= 0 && params.weir_fst_populations.size() > 0) vf->output_weir_and_cockerham_fst(params);
	else if (params.weir_fst_populations.size() > 0) vf->output_windowed_weir_and_cockerham_fst(params);

	if (params.output_indel_hist == true) vf->output_indel_hist(params);

	LOG.printLOG("After filtering, kept " + header::int2str(vf->N_kept_sites()) + " out of a possible " + header::int2str(vf->N_total_sites()) + " Sites\n");
	if (vf->N_total_sites() <= 0)
		LOG.warning("File does not contain any sites");
	else if (vf->N_kept_sites() <= 0)
		LOG.warning("No data left for analysis!");

	time(&end);
	double running_time = difftime(end,start);
	LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
	LOG.close();
	delete vf;
	return 0;
}



