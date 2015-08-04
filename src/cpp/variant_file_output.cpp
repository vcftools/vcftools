/*
 * variant_file_output.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */
#include "variant_file.h"

void variant_file::output_frequency(const parameters &params, bool output_counts)
{
	// Output statistics of frequency at each site
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Frequency Statistics.");

	LOG.printLOG("Outputting Frequency Statistics...\n");
	string output_file = params.output_prefix + ".frq";
	if (output_counts)
		output_file += ".count";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Frequency output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	if (params.suppress_allele_output == false)
	{
		out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:";
		if (output_counts)
			out << "COUNT}\n";
		else
			out << "FREQ}\n";
	}
	else
	{
		if (output_counts)
			out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\n";
		else
			out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n";
	}

	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	vector<char> variant_line;
	entry *e = get_entry_object();
	unsigned int aa_idx = 0;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		if (params.derived)
			e->parse_basic_entry(true, false, true);
		else
			e->parse_basic_entry(true);
		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();

		if (params.derived)
		{
			aa_idx = 0;
			string AA = e->get_INFO_value("AA");
			std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper);	// Comment this out if only want high quality sites.
			if ((AA == "?") || (AA == "."))
			{
				LOG.one_off_warning("\tWarning: Cannot output derived allele frequencies without Ancestral Alleles (AA)");
				continue;
			}
			else
			{
				bool found = false;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (AA == e->get_allele(ui))
					{
						aa_idx = ui; found = true;
						break;
					}
				}
				if (found == false)
				{
					LOG.one_off_warning("\tWarning: Ancestral allele does not match any SNP allele.");
					continue;
				}
			}
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << N_alleles << "\t" << N_non_missing_chr;
		if (output_counts)
		{
			if (params.suppress_allele_output == false)
			{
				out << "\t" << e->get_allele(aa_idx) << ":" << allele_counts[aa_idx];
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
						out << "\t" << e->get_allele(ui) << ":" << allele_counts[ui];
				}
				out << "\n";
			}
			else
			{
				out << "\t" << allele_counts[aa_idx];
				for (unsigned ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
						out << "\t" << allele_counts[ui];
				}
				out << "\n";
			}
		}
		else
		{
			double freq;
			if (params.suppress_allele_output == false)
			{
				freq = allele_counts[aa_idx] / (double)N_non_missing_chr;
				out << "\t" << e->get_allele(aa_idx) << ":" << freq;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
					{
						freq = allele_counts[ui] / (double)N_non_missing_chr;
						out << "\t" << e->get_allele(ui) << ":" << freq;
					}
				}
				out << "\n";
			}
			else
			{
				freq = allele_counts[aa_idx] / (double)N_non_missing_chr;
				out << "\t" << freq;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
					{
						freq = allele_counts[ui] / (double)N_non_missing_chr;
						out << "\t" << freq;
					}
				}
				out << "\n";
			}
		}
	}
	delete e;
}

void variant_file::output_het(const parameters &params)
{
	// Output statistics on Heterozygosity for each individual
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Heterozygosity Statistics.");
	// Following the calculations in PLINK....
	// Note this assumes Biallelic SNPs.
	LOG.printLOG("Outputting Individual Heterozygosity\n");
	streambuf * buf;
	ofstream temp_out;
	string output_file = params.output_prefix + ".het";
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Heterozygosity output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF" << endl;

	// P(Homo) = F + (1-F)P(Homo by chance)
	// P(Homo by chance) = p^2+q^2 for a biallelic locus.
	// For an individual with N genotyped loci, we
	//   1. count the total observed number of loci which are homozygous (O),
	//   2. calculate the total expected number of loci homozygous by chance (E)
	// Then, using the method of moments, we have
	//    O = NF + (1-F)E
	// Which rearranges to give
	//    F = (O-E)/(N-E)

	double freq;
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	vector<int> N_sites_included(meta_data.N_indv, 0);
	vector<int> N_obs_hom(meta_data.N_indv, 0);
	vector<double> N_expected_hom(meta_data.N_indv, 0.0);
	pair<int, int> alleles;

	vector<char> variant_line;
	entry *e = get_entry_object();
	while (!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);
		e->parse_basic_entry(true);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tIndividual Heterozygosity: Only using biallelic SNPs.");
			continue;
		}

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tIndividual Heterozygosity: Only using fully diploid SNPs.");
			continue;
		}

		// Frequency of non-reference allele
		e->get_allele_counts(allele_counts, N_non_missing_chr);

		if (N_non_missing_chr > 0)
			freq = allele_counts[1] / double(N_non_missing_chr);
		else
			freq = -1;

		if ((freq <= numeric_limits<double>::epsilon())  || (1.0 - freq <= numeric_limits<double>::epsilon()))
			continue;

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, alleles);
				if ((alleles.first > -1) && (alleles.second > -1))
				{
					N_sites_included[ui]++;
					if (alleles.first == alleles.second)
						N_obs_hom[ui]++;

					N_expected_hom[ui] += 1.0 - (2.0 * freq * (1.0 - freq) * (N_non_missing_chr / (N_non_missing_chr - 1.0)));
				}
			}
		}
	}

	out.setf(ios::fixed,ios::floatfield);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (N_sites_included[ui] > 0)
		{
			double F = (N_obs_hom[ui] - N_expected_hom[ui]) / double(N_sites_included[ui] - N_expected_hom[ui]);
			out << meta_data.indv[ui] << "\t" << N_obs_hom[ui] << "\t";
			out.precision(1);
			out << N_expected_hom[ui] << "\t";
			out.precision(5);
			out << N_sites_included[ui] << "\t" << F << endl;
		}
	}
	delete e;
}

void variant_file::output_hwe(const parameters &params)
{
	// Output HWE statistics for each site as described in Wigginton, Cutler, and Abecasis (2005)
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output HWE Statistics.");
	// Note this assumes Biallelic SNPs.
	LOG.printLOG("Outputting HWE statistics (but only for biallelic loci)\n");

	string output_file = params.output_prefix + ".hwe";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS" << endl;

	/* PLINK code:
	// b11 = Nhom1, b12 = Nhet, b22 = Nhom2
	double tot = b11 + b12 + b22;
	double exp_11 = freq * freq * tot;
	double exp_12 = 2 * freq * (1-freq) * tot;
	double exp_22 = (1-freq) * (1-freq) * tot;

	double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
			    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12
			    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;

	p = chiprobP(chisq,1);
	*/
	double freq;
	unsigned int b11, b12, b22;
	double exp_11, exp_12, exp_22;
	double chisq;
	double tot;
	double p_hwe, p_lo, p_hi;
	unsigned int precision = out.precision();
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	vector<char> variant_line;
	entry *e = get_entry_object();

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tHWE: Only using biallelic SNPs.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tHWE: Only using fully diploid SNPs.");
			continue;	// Isn't diploid
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		freq = allele_counts[0] / (double)N_non_missing_chr;
		e->get_genotype_counts(b11, b12, b22);
		tot = b11 + b12 + b22;
		exp_11 = freq * freq * tot;
		exp_12 = 2.0 * freq * (1.0-freq) * tot;
		exp_22 = (1.0-freq) * (1.0-freq) * tot;

		chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
				+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
				+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22;

		entry::SNPHWE(b12, b11, b22, p_hwe, p_lo, p_hi);
		out << e->get_CHROM() << "\t" << e->get_POS();
		out << "\t" << b11 << "/" << b12 << "/" << b22;
		out.precision(2);
		out << fixed << "\t" << exp_11 << "/" << exp_12 << "/" << exp_22;
		out.precision(precision);
		out << scientific;
		out << "\t" << chisq << "\t" << p_hwe << "\t" << p_lo << "\t" << p_hi << endl;
	}
	delete e;
}


void variant_file::output_indv_burden(const parameters &params)
{
	// Output the burden within each individual of variants at each frequency.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Burden Statistics.");
	LOG.printLOG("Outputting variant burden by individual\n");
	string output_file = params.output_prefix + ".iburden";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Burden Output File: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	vector< int > hom_ref_burden(meta_data.N_indv, 0);
	vector< int > het_burden(meta_data.N_indv, 0);
	vector< int > hom_alt_burden(meta_data.N_indv, 0);
	vector< int > missing_burden(meta_data.N_indv, 0);

	if (params.derived)
		out << "INDV\tN_HOM_ANC\tN_HET\tN_HOM_DER\tN_MISS" << endl;
	else
		out << "INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS" << endl;

	unsigned int N_alleles;
	vector<char> variant_line;
	entry *e = get_entry_object();
	int aa_idx = 0;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		if (params.derived)
			e->parse_basic_entry(true, false, true);
		else
			e->parse_basic_entry(true);

		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tWarning: Only using fully diploid sites.");
			continue;
		}

		if (params.derived)
		{
			aa_idx = 0;
			string AA = e->get_INFO_value("AA");
			std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper);	// Comment this out if only want high quality sites.
			if ((AA == "?") || (AA == "."))
			{
				LOG.one_off_warning("\tWarning: Cannot find Ancestral Alleles (AA)");
				continue;
			}
			else
			{
				bool found = false;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (AA == e->get_allele(ui))
					{
						aa_idx = ui;
						found = true;
						break;
					}
				}
				if (found == false)
				{
					LOG.one_off_warning("\tWarning: Ancestral allele does not match any SNP allele.");
					continue;
				}
			}
		}

		pair<int, int> geno;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (e->include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, geno);

				if ((geno.first == aa_idx) && (geno.second == aa_idx))
					hom_ref_burden[ui]++;
				else if ((geno.first >= 0) && (geno.second >= 0) && (geno.first != geno.second))
					het_burden[ui]++;
				else if ((geno.first >= 0) && (geno.second >= 0) && (geno.first == geno.second))
					hom_alt_burden[ui]++;
				else
					missing_burden[ui]++;
			}
		}
	}
	delete e;

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		out << meta_data.indv[ui];
		out << "\t" << hom_ref_burden[ui];
		out << "\t" << het_burden[ui];
		out << "\t" << hom_alt_burden[ui];
		out << "\t" << missing_burden[ui] << endl;
	}
}


void variant_file::output_indv_freq_burden(const parameters &params, int double_count_hom_alt)
{
	// Output the burden within each individual of variants at each frequency.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Frequency Burden Statistics.");
	LOG.printLOG("Outputting frequency burden by individual\n");
	string output_file = params.output_prefix + ".ifreqburden";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Frequency Burden Output File: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	int N = N_kept_individuals();
	int max_chr_count = N * 2;	// Assuming diploidy...
	vector< vector<int> > burden_matrix(N, vector<int>(max_chr_count+1, 0));

	out << "INDV";
	for (int i=0; i<=max_chr_count; i++)
		out << "\t" << LOG.int2str(i);
	out << endl;

	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	vector<char> variant_line;
	entry *e = get_entry_object();
	int aa_idx = 0;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		if (params.derived)
			e->parse_basic_entry(true, false, true);
		else
			e->parse_basic_entry(true);

		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();

		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tWarning: Only using fully diploid sites.");
			continue;
		}

		if (params.derived)
		{
			aa_idx = 0;
			string AA = e->get_INFO_value("AA");
			std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper);	// Comment this out if only want high quality sites.
			if ((AA == "?") || (AA == "."))
			{
				LOG.one_off_warning("\tWarning: Cannot find Ancestral Alleles (AA)");
				continue;
			}
			else
			{
				bool found = false;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (AA == e->get_allele(ui))
					{
						aa_idx = ui;
						found = true;
						break;
					}
				}
				if (found == false)
				{
					LOG.one_off_warning("\tWarning: Ancestral allele does not match any SNP allele.");
					continue;
				}
			}
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		pair<int, int> geno;

		int indv_count = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (e->include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, geno);

				if ((geno.first != aa_idx) && (geno.first >= 0))
					burden_matrix[indv_count][allele_counts[geno.first]]++;
				if ((double_count_hom_alt == 0) || (geno.first != geno.second))
				{	// Count the second allele if required
					if ((geno.second != aa_idx) && (geno.second >= 0))
						burden_matrix[indv_count][allele_counts[geno.second]]++;
				}
			}
			indv_count++;
		}
	}
	delete e;

	int indv_count = 0;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		out << meta_data.indv[indv_count];
		for (int i=0; i<=max_chr_count; i++)
			out << "\t" << LOG.int2str(burden_matrix[indv_count][i]);
		out << endl;
		indv_count++;
	}
}

void variant_file::output_individuals_by_mean_depth(const parameters &params)
{
	// Output information regarding the mean depth for each individual
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Individuals by Mean Depth Statistics.");

	LOG.printLOG("Outputting Mean Depth by Individual\n");
	string output_file = params.output_prefix + ".idepth";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Individual Depth Output File: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	out << "INDV\tN_SITES\tMEAN_DEPTH" << endl;
	vector<double> depth_sum(meta_data.N_indv, 0.0);
	vector<int> count(meta_data.N_indv, 0);
	int depth;
	vector<char> variant_line;
	entry *e = get_entry_object();
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (e->include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->parse_genotype_entry(ui, false, false, true);
				depth = e->get_indv_DEPTH(ui);
				if (depth >= 0)
				{
					depth_sum[ui] += depth;
					count[ui]++;
				}
			}
		}
	}

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		double mean_depth = depth_sum[ui] / count[ui];
		out << meta_data.indv[ui] << "\t" << count[ui] << "\t" << mean_depth << endl;
	}
	delete e;
}

void variant_file::output_SNP_density(const parameters &params)
{
	// Output SNP density (technically variant density)
	int bin_size = params.output_SNP_density_bin_size;
	if (bin_size <= 0)
		return;
	LOG.printLOG("Outputting SNP density\n");

	string output_file = params.output_prefix + ".snpden";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open SNP Density Output File: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	string CHROM; int POS;
	vector<char> variant_line;
	entry *e = get_entry_object();
	map<string, vector<int> > bins;
	vector<string> chrs;
	unsigned int idx;
	double C = 1.0 / double(bin_size);
	int prev_pos = -1; string prev_chrom = "";
	string alt;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		CHROM = e->get_CHROM();
		POS = e->get_POS();
		alt = e->get_ALT();

		if (alt != "." && (POS != prev_pos || CHROM != prev_chrom))
		{
			idx = (unsigned int)(POS * C);

			if (idx>=bins[CHROM].size())
				bins[CHROM].resize(idx+1,0);

			bins[CHROM][idx]++;
		}
		if (CHROM != prev_chrom)
			chrs.push_back(CHROM);

		prev_pos = POS;
		prev_chrom = CHROM;
	}

	out << "CHROM\tBIN_START\tSNP_COUNT\tVARIANTS/KB" << endl;
	int bin_tot;
	C = 1000.0 / bin_size;
	for (unsigned int ui=0; ui<chrs.size(); ui++)
	{
		bool output = false;
		CHROM = chrs[ui];
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			bin_tot = bins[CHROM][s];
			if (bin_tot > 0)
				output = true;
			if (output == true)
				out << CHROM << "\t" << s*bin_size << "\t" << bin_tot << "\t" << bin_tot * C << endl;
		}
	}
	delete e;
}

void variant_file::output_indv_missingness(const parameters &params)
{
	// Output missingness by individual
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Missingness Statistics.");

	LOG.printLOG("Outputting Individual Missingness\n");
	string output_file = params.output_prefix + ".imiss";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Individual Missingness Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	out << "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS" << endl;
	unsigned int ui;
	vector<unsigned int> indv_N_missing(meta_data.N_indv, 0), indv_N_tot(meta_data.N_indv, 0);
	vector<unsigned int> indv_N_geno_filtered(meta_data.N_indv, 0);
	pair<int, int> alleles;
	vector<char> variant_line;
	entry *e = get_entry_object();

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();
		for (ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (e->include_genotype[ui] == false)
			{
				indv_N_geno_filtered[ui]++;
				continue;
			}

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);
			if (alleles.first == -1)
				indv_N_missing[ui]++;
			indv_N_tot[ui]++;
		}
	}

	for (ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		out << meta_data.indv[ui] << "\t" << indv_N_tot[ui] << "\t";
		out << indv_N_geno_filtered[ui] << "\t" << indv_N_missing[ui] << "\t";
		out << indv_N_missing[ui] / double(indv_N_tot[ui]) << endl;
	}

	delete e;
}

void variant_file::output_site_missingness(const parameters &params)
{
	// Output missingness by site
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Missingness Statistics.");

	LOG.printLOG("Outputting Site Missingness\n");

	string output_file = params.output_prefix + ".lmiss";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Site Missingness Output File: " + output_file, 4);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	unsigned int ui;
	unsigned int site_N_missing, site_N_tot, site_N_geno_filtered;
	pair<int, int> alleles;
	vector<char> variant_line;
	entry *e = get_entry_object();

	out << "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS" << endl;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();

		site_N_missing = 0;
		site_N_tot = 0;
		site_N_geno_filtered = 0;
		for (ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (e->include_genotype[ui] == false)
			{
				site_N_geno_filtered++;
				continue;
			}

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);
			if (alleles.first == -1)
				site_N_missing++;
			if (alleles.second == -1)
				site_N_missing++;
			site_N_tot+=2;

			if ((alleles.second == -1) && (e->get_indv_PHASE(ui) == '|'))
			{	// Phased missing genotypes indicate haploid genome
				site_N_tot--;
				site_N_missing--;
			}
			else if (alleles.second == -2)
			{
				site_N_tot--;
				site_N_missing--;
			}
		}
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << site_N_tot << "\t" << site_N_geno_filtered << "\t";
		out << site_N_missing << "\t" << double(site_N_missing) / double(site_N_tot) << endl;
	}
	delete e;
}

void variant_file::calc_hap_r2(vector<pair<int, int> > &GT1, vector<pair<int,int> > &GT2, double &r2, double &D, double &Dprime, int &chr_count)
{
	double x11=0, x12=0, x21=0, x22=0;
	double X=0, X2=0, Y=0, Y2=0, XY=0;
	double sx, sy;
	double rel_x11, p1, p2, q1, q2, Dmax;
	double var1, var2, cov12;
	chr_count = 0;
	int allele1, allele2;

	for (unsigned int ui=0; ui<GT1.size(); ui++)
	{
		for (unsigned int c=0; c<2; c++)
		{
			if (c==0)
			{
				allele1 = GT1[ui].first;
				allele2 = GT2[ui].first;
			}
			else
			{
				allele1 = GT1[ui].second;
				allele2 = GT2[ui].second;
			}

			if ((allele1 < 0) || (allele2 < 0))
				continue;

			if (allele1 == 0 && allele2 == 0){
			  x11++;
			} else if (allele1 == 0 && allele2 != 0){
			  x12++;
			} else if (allele1 != 0 && allele2 == 0){
			  x21++;
			} else { // (allele1 !=0 && allele2 != 0)
			  x22++;
			}

			sx=0, sy=0;
			if (allele1 == 0)
				sx += 1;

			if (allele2 == 0)
				sy += 1;

			X += sx; Y += sy;
			XY += sx*sy;
			sx *= sx; sy *= sy;
			X2 += sx;
			Y2 += sy;

			chr_count++;
		}
	}

	rel_x11 = x11/double(chr_count);
	p1 = (x11 + x12)/double(chr_count);
	p2 = (x21 + x22)/double(chr_count);
	q1 = (x11 + x21)/double(chr_count);
	q2 = (x12 + x22)/double(chr_count);
	D = rel_x11 - p1*q1;
	if (D < 0)
		Dmax = min(p1*q1,p2*q2);
	else
		Dmax = min(p1*q2,p2*q1);
	Dprime = D/Dmax;

	X /= chr_count; X2 /= chr_count;
	Y /= chr_count;	Y2 /= chr_count;
	XY /= chr_count;
	var1 = X2 - X*X;
	var2 = Y2 - Y*Y;
	cov12 = XY - X*Y;

	r2 = cov12 * cov12 / (var1 * var2);
}

// Calculate r2 for either haplotypes or genotypes using the em algorithm...
void variant_file::calc_r2_em(entry *e, entry *e2, double &r2, int &indv_count)
{
	r2 = 0;	indv_count = 0;
	pair<int, int> geno1, geno2;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if ((include_indv[ui] == false) || (e->include_genotype[ui] == false) || (e2->include_genotype[ui] == false))
			continue;

		e->get_indv_GENOTYPE_ids(ui, geno1);
		e2->parse_genotype_entry(ui, true);
		e2->get_indv_GENOTYPE_ids(ui, geno2);
		// TODO... not yet implemented...!
		LOG.error("Not yet implmented!\n");
	}
}

void variant_file::calc_geno_r2(vector<pair<int,int> > &GT1, vector<pair<int,int> > &GT2, double &r2, int &indv_count)
{
	double X=0, X2=0, Y=0, Y2=0, XY=0;
	double sx, sy;
	indv_count = 0;
	pair<int, int> geno1, geno2;

	for (unsigned int ui=0; ui<GT1.size(); ui++)
	{
		geno1 = GT1[ui];
		geno2 = GT2[ui];

		if ((geno1.first < 0) || (geno1.second < 0))
			continue;

		if ((geno2.first < 0) || (geno2.second < 0))
			continue;

		sx=0, sy=0;
		if (geno1.first == geno1.second)
		{
			if (geno1.first == 0)
				sx = 2;
		}
		else
			sx = 1;

		if (geno2.first == geno2.second)
		{
			if (geno2.first == 0)
				sy = 2;
		}
		else
			sy = 1;

		X += sx; Y += sy;
		XY += sx*sy;
		sx *= sx; sy *= sy;
		X2 += sx; Y2 += sy;

		indv_count++;
	}

	X /= indv_count; X2 /= indv_count;
	Y /= indv_count; Y2 /= indv_count;
	XY /= indv_count;

	double var1 = X2 - X*X;
	double var2 = Y2 - Y*Y;
	double cov12 = XY - X*Y;

	r2 = cov12 * cov12 / (var1 * var2);
}

void variant_file::calc_geno_chisq(vector<pair<int,int> > &GT1, vector<pair<int,int> > &GT2, int &N0, int &N1, double &chisq, double &dof, double &pval, int &indv_count)
{
	int N_genotypes0 = N0 * (N0+1) / 2;
	int N_genotypes1 = N1 * (N1+1) / 2;

	vector<vector<double> > observed(N_genotypes0, vector<double>(N_genotypes1,0));
	indv_count = 0;
	pair<int, int> geno1, geno2;
	for (unsigned int ui=0; ui<GT1.size(); ui++)
	{
		geno1 = GT1[ui];
		geno2 = GT2[ui];
		if ((geno1.first < 0) || (geno1.second < 0))
			continue;

		if ((geno2.first < 0) || (geno2.second < 0))
			continue;

		map<pair<int, int>, int> idx_lookup1;
		int count = 0;
		for (int uj=0; uj<N0; uj++)
		{
			for (int uk=uj; uk<N0; uk++)
			{
				idx_lookup1[make_pair(uj,uk)] = count;
				idx_lookup1[make_pair(uk,uj)] = count;
				count++;
			}
		}
		map<pair<int, int>, int> idx_lookup2;
		count = 0;
		for (int uj=0; uj<N1; uj++)
		{
			for (int uk=uj; uk<N1; uk++)
			{
				idx_lookup2[make_pair(uj,uk)] = count;
				idx_lookup2[make_pair(uk,uj)] = count;
				count++;
			}
		}

		int idx1 = idx_lookup1[geno1];
		int idx2 = idx_lookup2[geno2];

		observed[idx1][idx2]++;
		indv_count++;
	}

	vector<vector<double> > expected(N_genotypes0, vector<double>(N_genotypes1,0));
	vector<double> row_tot(N_genotypes0, 0);
	vector<double> col_tot(N_genotypes1, 0);
	double tot=0;
	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			row_tot[ui] += observed[ui][uj];
			col_tot[uj] += observed[ui][uj];
			tot += observed[ui][uj];
		}
	}

	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			expected[ui][uj] = row_tot[ui] * col_tot[uj] / tot;
		}
	}

	chisq = 0;
	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			if ((row_tot[ui] > 0) && (col_tot[uj] > 0))	// Don't use incomplete cases
				chisq += pow(observed[ui][uj] - expected[ui][uj], 2) / expected[ui][uj];
		}
	}
	int n_col=0, n_row=0;
	for (int ui=0; ui<N_genotypes0; ui++)
		if (row_tot[ui] > 0)
			n_row++;
	for (int ui=0; ui<N_genotypes1; ui++)
		if (col_tot[ui] > 0)
			n_col++;

	dof = (n_row-1) * (n_col-1);
	pval = 1.0-gammp(dof/2, chisq/2);
}

// Count the number of haplotypes within user-defined bins
void variant_file::output_haplotype_count(const parameters &params)
{
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Haplotype Counts.");

	LOG.printLOG("Outputting Haplotype Counts\n");

	ifstream BED(params.hapcount_BED.c_str());
	if (!BED.is_open())
		LOG.error("Could not open BED file: " + params.hapcount_BED);

	string line; stringstream ss;
	string CHROM; int POS1, POS2; int idx;
	unsigned int N_chr=0;
	BED.ignore(numeric_limits<streamsize>::max(), '\n');;
	vector< vector< pair<int, int> > > bin_positions;
	map<string, int> chr_to_idx;
	while (!BED.eof())
	{
		getline(BED, line);
		if ((line[0] == '#') || (line.size() == 0))
			continue;

		ss.clear();
		ss.str(line);
		ss >> CHROM >> POS1 >> POS2;

		if (chr_to_idx.find(CHROM) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[CHROM] = (N_chr-1);
			bin_positions.resize(N_chr);
		}

		idx = chr_to_idx[CHROM];
		bin_positions[idx].push_back(make_pair(POS1, POS2));
	}
	BED.close();

	for (unsigned int ui=0; ui<bin_positions.size(); ui++)
	{
		sort(bin_positions[ui].begin(), bin_positions[ui].end());
		for (unsigned int uj=1; uj<bin_positions[ui].size(); uj++)
		{	// Very naive check to ensure BED file is non-overlapping
			if (bin_positions[ui][uj-1].second > bin_positions[ui][uj].first)
				LOG.error("BED file must be non-overlapping.\n", 33);
		}
	}

	vector< vector<int> > haplotypes(2*meta_data.N_indv);
	vector<char> variant_line;
	entry *e = get_entry_object();
	string haplotype;
	pair<int, int> geno;
	vector<unsigned int> min_ui(N_chr, 0);
	bool have_data=false;

	string output_file = params.output_prefix + ".hapcount";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Haplotype Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "#CHROM\tBIN_START\tBIN_END\tN_SNP\tN_UNIQ_HAPS\tN_GROUPS\t{MULTIPLICITY:FREQ}" << endl;
	int bin_idx=0, prev_bin_idx=-1;
	int prev_idx = -1;
	vector<int> haplotype_count;
	vector<int> SNP_count;
	vector< map<int, int> > haplotype_frequencies;
	string prev_CHROM="";

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		CHROM = e->get_CHROM();
		POS1 = e->get_POS();

		if (chr_to_idx.find(CHROM) == chr_to_idx.end())
			continue;

		idx = chr_to_idx[CHROM];
		if (idx != prev_idx)
		{	// Moved to a new chromosome, so output last chromosome
			map<vector<int>, int > haplotype_set;
			if (have_data == true)
			{	// Process any remaining data
				for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
				{
					if (e->include_indv[ui] == false)
						continue;

					haplotype_set[haplotypes[(2*ui)]]++;
					haplotype_set[haplotypes[(2*ui)+1]]++;
					SNP_count[prev_bin_idx] = haplotypes[2*ui].size();
					haplotypes[(2*ui)].resize(0);
					haplotypes[(2*ui)+1].resize(0);
				}
				haplotype_count[prev_bin_idx] = haplotype_set.size();
				for (map<vector<int>, int >::iterator it=haplotype_set.begin(); it != haplotype_set.end(); ++it)
					haplotype_frequencies[prev_bin_idx][it->second]++;
			}
			have_data = false;

			for (unsigned int ui=0; ui<haplotype_count.size(); ui++)
			{
				out << prev_CHROM << "\t" << bin_positions[prev_idx][ui].first << "\t" << bin_positions[prev_idx][ui].second;
				out << "\t" << SNP_count[ui] << "\t" << haplotype_count[ui];
				out << "\t" << haplotype_frequencies[ui].size();
				for (map<int, int>::iterator it = haplotype_frequencies[ui].begin(); it != haplotype_frequencies[ui].end(); ++it)
					out << "\t" << it->second << ":" << it->first;
				out << endl;
			}

			// Set up for new chromosome
			unsigned int N_bins = bin_positions[idx].size();
			haplotype_count.clear(); haplotype_count.resize(N_bins, 0);
			SNP_count.clear(); SNP_count.resize(N_bins, 0);
			haplotype_frequencies.clear(); haplotype_frequencies.resize(N_bins);
			bin_idx=0, prev_bin_idx=-1;
			prev_idx = idx;
			prev_CHROM = CHROM;
		}

		bool found=false;
		unsigned int max_ui = bin_positions[idx].size();
		for (unsigned int ui=min_ui[idx]; ui<max_ui; ui++)
		{	// No need to start this loop at zero every time...
			if ((POS1 > bin_positions[idx][ui].first) && (POS1 <= bin_positions[idx][ui].second))
			{	// We're in a BED bin, so add to haplotypes
				found=true;
				prev_bin_idx = bin_idx;
				bin_idx = ui;
				break;
			}
			else if (POS1 > bin_positions[idx][ui].second)
				min_ui[idx] = ui+1;
		}

		if ((found == false) || (prev_bin_idx != bin_idx))
		{	// Changed bin, so update haplotype count in previous bin, and reset for next bin
			if (have_data == true)
			{
				map<vector<int>, int > haplotype_set;
				for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
				{
					if (e->include_indv[ui] == false)
						continue;

					haplotype_set[haplotypes[(2*ui)]]++;
					haplotype_set[haplotypes[(2*ui)+1]]++;
					SNP_count[prev_bin_idx] = haplotypes[2*ui].size();
					haplotypes[(2*ui)].resize(0);
					haplotypes[(2*ui)+1].resize(0);
				}
				haplotype_count[prev_bin_idx] = haplotype_set.size();
				for (map<vector<int>, int >::iterator it=haplotype_set.begin(); it != haplotype_set.end(); ++it)
					haplotype_frequencies[prev_bin_idx][it->second]++;
			}
			have_data = false;
		}

		if (found == true)
		{	// Inside a bin, so append to haplotypes
			have_data = true;
			e->parse_genotype_entries(true);

			if (e->is_diploid() == false)
			{
				LOG.one_off_warning("\tWarning: Only using fully diploid sites.");
				continue;
			}

			for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
			{
				if (e->include_indv[ui] == false)
					continue;

				geno.first = -1; geno.second = -1;
				if (e->include_genotype[ui] == true)
					e->get_indv_GENOTYPE_ids(ui, geno);

				haplotypes[(2*ui)].push_back(geno.first);
				haplotypes[(2*ui)+1].push_back(geno.second);
			}
		}
	}
	delete e;

	if (idx == prev_idx)
	{	// Output any remaining data from last chromosome
		if (have_data == true)
		{	// Process any remaining data
			map<vector<int>, int > haplotype_set;
			for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
			{
				if (e->include_indv[ui] == false)
					continue;

				haplotype_set[haplotypes[(2*ui)]]++;
				haplotype_set[haplotypes[(2*ui)+1]]++;
				SNP_count[prev_bin_idx] = haplotypes[2*ui].size();
			}
			haplotype_count[prev_bin_idx] = haplotype_set.size();
			for (map<vector<int>, int >::iterator it=haplotype_set.begin(); it != haplotype_set.end(); ++it)
				haplotype_frequencies[prev_bin_idx][it->second]++;
		}

		for (unsigned int ui=0; ui<haplotype_count.size(); ui++)
		{
			out << CHROM << "\t" << bin_positions[prev_idx][ui].first << "\t" << bin_positions[prev_idx][ui].second << "\t";
			out << SNP_count[ui] << "\t" << haplotype_count[ui];
			out << "\t" << haplotype_frequencies[ui].size();
			for (map<int, int>::iterator it = haplotype_frequencies[ui].begin(); it != haplotype_frequencies[ui].end(); ++it)
				out << "\t" << it->second << ":" << it->first;
			out << endl;
		}
	}
}

void variant_file::output_haplotype_r2(const parameters &params)
{
	// Output pairwise LD statistics, using traditional r^2. Requires phased haplotypes.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	int snp_window_size = params.ld_snp_window_size;
	int snp_window_min = params.ld_snp_window_min;
	int bp_window_size = params.ld_bp_window_size;
	int bp_window_min = params.ld_bp_window_min;
	double min_r2 = params.min_r2;

	LOG.printLOG("Outputting Pairwise LD (phased bi-allelic only)\n");
	string output_file = params.output_prefix + ".hap.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime" << endl;

	double r2, D, Dprime;
	int chr_count, site_count = 0;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	pair<int, int> geno;
	vector<char> variant_line;
	string CHROM,CHROM2;
	int POS,POS2;
	vector<char> out_line, tmp_int;
	entry *e = get_entry_object();

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	out_line.reserve(meta_data.N_indv+10);
	int indv_miss = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);
		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();

		site_count++;
		string chrom_str = CHROM+"\n";

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (e->get_indv_PHASE(ui) != '|')
			{
				remove(tmpname);
				LOG.error("Require phased haplotypes for r^2 calculation (use --phased)\n");
			}
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}

			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) > 2)
			{
				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Cannot use polyploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		fd.write(&out_line[0],out_line.size());
	}
	fd.close();

	if (N_kept_entries <= 1)
	{
		remove(tmpname);
		LOG.error("Insufficient sites remained after filtering");
	}

	ifstream tmp_file(tmpname, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;

	unsigned int uj = 0;
	for(unsigned int ui=0; ui<site_count-1; ui++)
	{
		tmp_file.seekg(file_pos);
		GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
		read_temp_site(tmp_file, CHROM, POS, GTs);
		file_pos = tmp_file.tellg();

		for(uj=ui+1; uj<site_count; uj++)
		{
			if (int(uj-ui) > snp_window_size)
				break;
			GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if(uj < (ui+skip))
				continue;
			if (CHROM != CHROM2)
				continue;
			if (POS2 < POS)
				LOG.one_off_warning("Warning: Input is unsorted, results may not be complete.");
			if ((POS2 - POS) < bp_window_min)
				continue;
			if ((POS2 - POS) > bp_window_size)
				break;

			calc_hap_r2(GTs, GTs2, r2, D, Dprime, chr_count);
			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << POS2 << "\t" << chr_count << "\t" << r2 << "\t" << D << "\t" << Dprime << "\t" << endl;
		}
	}
	tmp_file.close();
	remove(tmpname);
	delete e;
}

void variant_file::output_genotype_r2(const parameters &params)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	int snp_window_size = params.ld_snp_window_size;
	int snp_window_min = params.ld_snp_window_min;
	int bp_window_size = params.ld_bp_window_size;
	int bp_window_min = params.ld_bp_window_min;
	double min_r2 = params.min_r2;

	LOG.printLOG("Outputting Pairwise LD (bi-allelic only)\n");
	string output_file = params.output_prefix + ".geno.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR\tPOS1\tPOS2\tN_INDV\tR^2" << endl;

	double r2;
	int indv_count;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	vector<char> variant_line;
	entry *e = get_entry_object();
	int count = 0;
	string CHROM, CHROM2;
	int POS, POS2;
	pair<int, int> geno;
	vector<char> out_line, tmp_int;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	out_line.reserve(meta_data.N_indv+10);
	int indv_miss = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tgenoLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		count++;
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();
		string chrom_str = CHROM+"\n";

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}

			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) != 2)
			{
				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Only using diploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		fd.write(&out_line[0],out_line.size());
	}
	fd.close();

	if (N_kept_entries <= 1)
	{
		remove(tmpname);
		LOG.error("Insufficient sites remained after filtering");
	}

	ifstream tmp_file(tmpname, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;

	unsigned int uj = 0;
	for(unsigned int ui=0; ui<count-1; ui++)
	{
		tmp_file.seekg(file_pos);
		GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
		read_temp_site(tmp_file, CHROM, POS, GTs);
		file_pos = tmp_file.tellg();

		for(uj=ui+1; uj<count; uj++)
		{
			if (int(uj-ui) > snp_window_size)
				break;

			GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if(uj < (ui+skip))
				continue;
			if (CHROM != CHROM2)
				continue;
			if (POS2 < POS)
				LOG.one_off_warning("Warning: Input is unsorted, results may not be complete.");
			if ((POS2 - POS) < bp_window_min)
				continue;
			if ((POS2 - POS) > bp_window_size)
				break;

			calc_geno_r2(GTs, GTs2, r2, indv_count);
			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << POS2 << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	remove(tmpname);
	delete e;
}

void variant_file::output_genotype_chisq(const parameters &params, double min_pval)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	int snp_window_size = params.ld_snp_window_size;
	int snp_window_min = params.ld_snp_window_min;
	int bp_window_size = params.ld_bp_window_size;
	int bp_window_min = params.ld_bp_window_min;

	LOG.printLOG("Outputting Pairwise LD\n");
	string output_file = params.output_prefix + ".geno.chisq";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR\tPOS1\tPOS2\tN_INDV\tCHI^2\tDOF\tPVAL" << endl;

	double chisq, dof, pval;
	int indv_count;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	vector<char> variant_line;
	entry *e = get_entry_object();
	string CHROM, CHROM2;
	int POS, POS2;
	pair<int,int> geno;
	int8_t tmp_alleles;
	int alleles, alleles2;
	vector<char> out_line, tmp_int;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	out_line.reserve(2*meta_data.N_indv+11);
	int indv_miss = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;

		N_kept_entries++;
		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();
		tmp_alleles = (int8_t)e->get_N_alleles();
		string chrom_str = CHROM+"\n";

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));
		out_line.push_back(tmp_alleles);

		int8_t out_byte = 0x00;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}
			if (e->include_genotype[ui] == false)
			{
				out_byte = 0xFF;
				out_line.push_back(out_byte);
				out_line.push_back(out_byte);
				continue;
			}

			if (e->get_indv_ploidy(ui) != 2)
			{
				out_byte = 0xFF;
				out_line.push_back(out_byte);
				out_line.push_back(out_byte);
				LOG.one_off_warning("\tgenoLD: Only using diploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte = 0xFF;
			else
				out_byte = (int8_t)geno.first;
			out_line.push_back(out_byte);

			if (geno.second == -1)
				out_byte = 0xFF;
			else
				out_byte = (int8_t)geno.second;
			out_line.push_back(out_byte);
		}
		fd.write(&out_line[0],out_line.size());
	}
	fd.close();

	if (N_kept_entries <= 1)
	{
		remove(tmpname);
		LOG.error("Insufficient sites remained after filtering");
	}

	ifstream tmp_file(tmpname, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;

	unsigned int uj = 0;
	for(unsigned int ui=0; ui<N_kept_entries-1; ui++)
	{
		tmp_file.seekg(file_pos);
		GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
		read_big_temp_site(tmp_file, CHROM, POS, alleles, GTs);
		file_pos = tmp_file.tellg();

		for(uj=ui+1; uj<N_kept_entries; uj++)
		{
			if (int(uj-ui) > snp_window_size)
				break;

			GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
			read_big_temp_site(tmp_file, CHROM2, POS2, alleles2, GTs2);

			if(uj < (ui+skip))
				continue;
			if (CHROM != CHROM2)
				continue;
			if (POS2 < POS)
				LOG.one_off_warning("Warning: Input is unsorted, results may not be complete.");
			if ((POS2 - POS) < bp_window_min)
				continue;
			if ((POS2 - POS) > bp_window_size)
				break;
			calc_geno_chisq(GTs, GTs2, alleles, alleles2, chisq, dof, pval, indv_count);

			if (min_pval > 0)
				if ((pval < min_pval) | (pval != pval))
					continue;
			out << CHROM << "\t" << POS << "\t" << POS2 << "\t" << indv_count << "\t" << chisq << "\t" << dof << "\t" << pval << endl;
		}
	}
	remove(tmpname);
	delete e;
}

void variant_file::output_interchromosomal_genotype_r2(const parameters &params)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	double min_r2 = params.min_r2;

	LOG.printLOG("Outputting Interchromosomal Pairwise Genotype LD (bi-allelic only)\n");
	string output_file = params.output_prefix + ".interchrom.geno.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
			buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2" << endl;

	int indv_count;
	double r2;
	vector<char> variant_line;
	entry *e = get_entry_object();
	int count = 0;
	string CHROM, CHROM2;
	int POS, POS2;
	pair<int,int> geno;
	vector<char> out_line, tmp_int;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	out_line.reserve(meta_data.N_indv+10);
	int indv_miss = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		count++;
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();
		string chrom_str = CHROM+"\n";

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}
			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) != 2)
			{
				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Only using diploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		fd.write(&out_line[0],out_line.size());
	}
	fd.close();

	if (N_kept_entries <= 1)
	{
		remove(tmpname);
		LOG.error("Insufficient sites remained after filtering");
	}

	ifstream tmp_file(tmpname, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;
	unsigned int uj=0;

	for(unsigned int ui=0; ui<count-1; ui++)
	{
		tmp_file.seekg(file_pos);
		GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
		read_temp_site(tmp_file, CHROM, POS, GTs);
		file_pos = tmp_file.tellg();

		for(uj=ui+1; uj<count; uj++)
		{
			GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if (CHROM == CHROM2)
				continue;

			calc_geno_r2(GTs, GTs2, r2, indv_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << CHROM2 << "\t" << POS2 << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	remove(tmpname);
	delete e;
}

void variant_file::output_interchromosomal_haplotype_r2(const parameters &params)
{
	double min_r2 = params.min_r2;
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	LOG.printLOG("Outputting Interchromosomal Pairwise LD (bi-allelic only)\n");
	string output_file = params.output_prefix + ".interchrom.hap.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2" << endl;

	double D, Dprime;
	int chr_count, site_count = 0;
	double r2;
	entry *e;
	e = get_entry_object();
	pair<int, int> geno;
	vector<char> variant_line;
	string CHROM,CHROM2;
	int POS,POS2;
	vector<char> out_line, tmp_int;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	out_line.reserve(meta_data.N_indv+10);
	int indv_miss = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();
		site_count++;
		string chrom_str = CHROM+"\n";

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (e->get_indv_PHASE(ui) != '|')
			{
				remove(tmpname);
				LOG.error("Require phased haplotypes for r^2 calculation (use --phased)\n");
			}
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}
			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) > 2)
			{
				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Cannot use polyploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		fd.write(&out_line[0],out_line.size());
	}
	fd.close();

	if (N_kept_entries <= 1)
	{
		remove(tmpname);
		LOG.error("Insufficient sites remained after filtering");
	}

	ifstream tmp_file(tmpname, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;

	unsigned int uj = 0;
	for(unsigned int ui=0; ui<site_count-1; ui++)
	{
		tmp_file.seekg(file_pos);
		GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
		read_temp_site(tmp_file, CHROM, POS, GTs);
		file_pos = tmp_file.tellg();

		for(uj=ui+1; uj<site_count; uj++)
		{
			GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if (CHROM == CHROM2)
				continue;

			calc_hap_r2(GTs, GTs2, r2, D, Dprime, chr_count);
			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << CHROM2 << "\t" << POS2 << "\t" << chr_count << "\t" << r2 << endl;
		}
	}
	remove(tmpname);
	delete e;
}

void variant_file::output_haplotype_r2_of_SNP_list_vs_all_others(const parameters &params)
{
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
			LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	LOG.printLOG("Outputting haplotype pairwise LD (bi-allelic only) for a set of SNPs versus all others.\n");

	int snp_window_size = params.ld_snp_window_size;
	int snp_window_min = params.ld_snp_window_min;
	int bp_window_size = params.ld_bp_window_size;
	int bp_window_min = params.ld_bp_window_min;

	string positions_file = params.hap_rsq_position_list;
	double min_r2 = params.min_r2;
	vector< set<int > > keep_positions;
	vector<int> list_positions;
	map<string, int> chr_to_idx;
	string line;
	stringstream ss;
	pair<int, int> geno;
	string CHROM, CHROM2;
	int POS, POS2, idx;
	unsigned int N_chr=0;
	vector<char> out_line, tmp_int;

	ifstream BED(positions_file.c_str());
	if (!BED.is_open())
		LOG.error("Could not open Positions file: " + positions_file);

	BED.ignore(numeric_limits<streamsize>::max(), '\n');
	int nlist = 0;
	while (!BED.eof())
	{
		getline(BED, line);
		if (line[0] == '#' || line == "")
			continue;

		ss.clear();
		ss.str(line);
		ss >> CHROM >> POS;

		if (chr_to_idx.find(CHROM) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[CHROM] = (N_chr-1);
			keep_positions.resize(N_chr);
		}

		idx = chr_to_idx[CHROM];
		keep_positions[idx].insert(POS);
		nlist += 1;
	}
	BED.close();

	if (nlist == 0)
		LOG.error("No sites found in positions file.\n",0);

	LOG.printLOG("\tRead in "+header::int2str(nlist)+" site(s) for LD analysis.\n");

	string output_file = params.output_prefix + ".list.hap.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2" << endl;

	double D, Dprime;
	int chr_count, site_count = 0;
	double r2;
	vector<char> variant_line;
	entry *e = get_entry_object();

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	char tmpname2[new_tmp.size()];
	strcpy(tmpname2, new_tmp.c_str());
	ret = mktemp(tmpname2);
	ofstream fd_POS(tmpname2, std::ios::out | std::ios::binary);
	if (!fd_POS.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	nlist = 0;
	int indv_miss = 0;
	out_line.reserve(meta_data.N_indv+10);
	while (!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();

		bool check_pos = false;
		if ( chr_to_idx.find(CHROM) != chr_to_idx.end() )
			if ( keep_positions[ chr_to_idx[CHROM] ].find(POS) != keep_positions[ chr_to_idx[CHROM] ].end() )
				check_pos = true;

		string chrom_str = CHROM+"\n";

		if (check_pos)
		{
			nlist++;
			list_positions.push_back(site_count);
		}
		else
			site_count++;

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (e->get_indv_PHASE(ui) != '|')
			{
				remove(tmpname);
				remove(tmpname2);
				LOG.error("Require phased haplotypes for r^2 calculation (use --phased)\n");
			}
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}
			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) > 2)
			{

				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Cannot use polyploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		if (check_pos)
			fd_POS.write(&out_line[0],out_line.size());
		else
			fd.write(&out_line[0],out_line.size());
	}
	fd.close();
	fd_POS.close();

	ifstream tmp_file(tmpname, ios::binary);
	ifstream tmp_file2(tmpname2, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;
	streampos file_pos2 = 0;

	GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
	GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));

	for(unsigned int ui=0; ui<nlist; ui++)
	{
		tmp_file2.seekg(file_pos2);
		read_temp_site(tmp_file2, CHROM, POS, GTs);
		file_pos2 = tmp_file2.tellg();

		tmp_file.seekg(0);
		for(unsigned int uj=0; uj<site_count; uj++)
		{
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if (CHROM != CHROM2)
				continue;
			if ( abs(POS2 - POS) < bp_window_min)
				continue;
			if ( abs(POS2 - POS) > bp_window_size)
				continue;

			int list_pos = list_positions[ui];
			if ( fabs(list_pos - uj) < snp_window_min)
				continue;
			if ( fabs(list_pos - uj) > snp_window_size)
				continue;

			calc_hap_r2(GTs, GTs2, r2, D, Dprime, chr_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << CHROM2 << "\t" << POS2 << "\t" << chr_count << "\t" << r2 << endl;
		}
	}
	remove(tmpname);
	remove(tmpname2);
	delete e;
}

void variant_file::output_genotype_r2_of_SNP_list_vs_all_others(const parameters &params)
{
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
			LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	LOG.printLOG("Outputting genotype pairwise LD (bi-allelic only) for a set of SNPs versus all others.\n");

	int snp_window_size = params.ld_snp_window_size;
	int snp_window_min = params.ld_snp_window_min;
	int bp_window_size = params.ld_bp_window_size;
	int bp_window_min = params.ld_bp_window_min;

	vector< set<int > > keep_positions;
	vector<int> list_positions;
	map<string, int> chr_to_idx;
	string line;
	stringstream ss;
	string CHROM, CHROM2;
	int POS, POS2, idx;
	pair<int,int> geno;
	unsigned int N_chr=0;
	double min_r2 = params.min_r2;
	vector<char> out_line, tmp_int;

	ifstream BED(params.geno_rsq_position_list.c_str());
	if (!BED.is_open())
		LOG.error("Could not open Positions file: " + params.geno_rsq_position_list);

	BED.ignore(numeric_limits<streamsize>::max(), '\n');
	int nlist = 0;
	while (!BED.eof())
	{
		getline(BED, line);
		if (line[0] == '#' || line == "")
			continue;

		ss.clear();
		ss.str(line);
		ss >> CHROM >> POS;

		if (chr_to_idx.find(CHROM) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[CHROM] = (N_chr-1);
			keep_positions.resize(N_chr);
		}

		idx = chr_to_idx[CHROM];
		keep_positions[idx].insert(POS);
		nlist++;
	}
	BED.close();

	if (nlist == 0)
		LOG.error("No sites found in positions file.\n",0);

	LOG.printLOG("\tRead in "+header::int2str(nlist)+" site(s) for LD analysis.\n");

	string output_file = params.output_prefix + ".list.geno.ld";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LD Output File: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2" << endl;

	int indv_count, site_count = 0;
	double r2;
	vector<char> variant_line;
	entry *e = get_entry_object();

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	char *ret = mktemp(tmpname);
	ofstream fd(tmpname, std::ios::out | std::ios::binary);
	if (!fd.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	char tmpname2[new_tmp.size()];
	strcpy(tmpname2, new_tmp.c_str());
	ret = mktemp(tmpname2);
	ofstream fd_POS(tmpname2, std::ios::out | std::ios::binary);
	if (!fd_POS.is_open())
		LOG.error(" Could not open temporary file.\n", 12);

	nlist = 0;
	int indv_miss = 0;
	out_line.reserve(meta_data.N_indv+10);
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);
		CHROM = e->get_CHROM();
		POS = e->get_POS();

		bool check_pos = false;
		if ( chr_to_idx.find(CHROM) != chr_to_idx.end() )
			if ( keep_positions[ chr_to_idx[CHROM] ].find(POS) != keep_positions[ chr_to_idx[CHROM] ].end() )
				check_pos = true;

		string chrom_str = CHROM+"\n";

		if (check_pos)
		{
			nlist++;
			list_positions.push_back(site_count);
		}
		else
			site_count++;

		out_line.resize(0);
		copy(chrom_str.begin(), chrom_str.end(), back_inserter(out_line));

		tmp_int.resize(0);
		e->make_int(tmp_int, POS, 3);
		copy(tmp_int.begin(), tmp_int.end(), back_inserter(out_line));

		char out_byte;
		indv_miss = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			out_byte = 0x00;
			if (include_indv[ui] == false)
			{
				indv_miss++;
				continue;
			}
			if (e->include_genotype[ui] == false)
			{
				out_line.push_back(0x22);
				continue;
			}

			if (e->get_indv_ploidy(ui) != 2)
			{
				out_line.push_back(0x22);
				LOG.one_off_warning("\tLD: Only using diploid individuals.");
				continue;
			}

			e->get_indv_GENOTYPE_ids(ui, geno);
			if (geno.first == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.first;
			out_byte = out_byte << 4;
			if (geno.second == -1)
				out_byte |= 0x02;
			else
				out_byte |= (char)geno.second;
			out_line.push_back(out_byte);
		}
		if (check_pos)
			fd_POS.write(&out_line[0],out_line.size());
		else
			fd.write(&out_line[0],out_line.size());
	}
	fd.close();
	fd_POS.close();

	ifstream tmp_file(tmpname, ios::binary);
	ifstream tmp_file2(tmpname2, ios::binary);
	vector<pair<int,int> > GTs, GTs2;
	streampos file_pos = 0;
	streampos file_pos2 = 0;

	GTs.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
	GTs2.resize(meta_data.N_indv-indv_miss, make_pair(-1,-1));
	for(unsigned int ui=0; ui<nlist; ui++)
	{
		tmp_file2.seekg(file_pos2);
		read_temp_site(tmp_file2, CHROM, POS, GTs);
		file_pos2 = tmp_file2.tellg();

		tmp_file.seekg(0);
		for(unsigned int uj=0; uj<site_count; uj++)
		{
			read_temp_site(tmp_file, CHROM2, POS2, GTs2);

			if (CHROM != CHROM2)
				continue;
			if ( abs(POS2 - POS) < bp_window_min)
				continue;
			if ( abs(POS2 - POS) > bp_window_size)
				continue;

			int list_pos = list_positions[ui];
			if ( fabs(list_pos - uj) < snp_window_min)
				continue;
			if ( fabs(list_pos - uj) > snp_window_size)
				continue;

			calc_geno_r2(GTs, GTs2, r2, indv_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << CHROM << "\t" << POS << "\t" << CHROM2 << "\t" << POS2 << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	remove(tmpname);
	remove(tmpname2);
	delete e;
}

void variant_file::output_singletons(const parameters &params)
{
	// Locate and output singletons (and private doubletons)
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Singletons.");

	LOG.printLOG("Outputting Singleton Locations\n");
	string output_file = params.output_prefix + ".singletons";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Singleton output file: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	out << "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV" << endl;

	int a;
	vector<int> allele_counts;
	unsigned int N_non_missing_chr, N_alleles, ui;
	pair<int, int> geno;
	string allele;
	vector<char> variant_line;
	entry *e = get_entry_object();
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		N_alleles = e->get_N_alleles();

		for (a=0; a<(signed)N_alleles; a++)
		{
			if (allele_counts[a] == 1)
			{	// Singleton
				for (ui=0; ui<meta_data.N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e->get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) || (geno.second == a))
					{
						e->get_allele(a, allele);
						out << e->get_CHROM() << "\t" << e->get_POS() << "\tS\t" << allele << "\t" << meta_data.indv[ui] << endl;
						ui=meta_data.N_indv;
						break;
					}
				}
			}
			else if (allele_counts[a] == 2)
			{	// Possible doubleton
				for (ui=0; ui<meta_data.N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e->get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) && (geno.second == a))
					{
						e->get_allele(a, allele);
						out << e->get_CHROM() << "\t" << e->get_POS() << "\tD\t" << allele << "\t" << meta_data.indv[ui] << endl;
						ui=meta_data.N_indv;
						break;
					}
				}
			}
		}
	}
	delete e;
}

void variant_file::output_genotype_depth(const parameters &params)
{
	// Output genotype depth in tab-delimited format.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Genotype Depth Statistics.");

	LOG.printLOG("Outputting Depth for Each Genotype\n");
	string output_file = params.output_prefix + ".gdepth";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Genotype Depth Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	out << "CHROM\tPOS";
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		out << "\t" << meta_data.indv[ui];
	}
	out << endl;

	vector<char> variant_line;
	entry *e = get_entry_object();
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();

		out << e->get_CHROM() << "\t" << e->get_POS();

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->parse_genotype_entry(ui, false, false, true);
				out << "\t" << e->get_indv_DEPTH(ui);
			}
			else
				out << "\t-1";
		}
		out << endl;
	}
	delete e;
}

void variant_file::output_FILTER_summary(const parameters &params)
{
	// Output a summary of sites in various FILTER categories.
	LOG.printLOG("Outputting Filter Summary (for bi-allelic loci only)\n");

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0;
	model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2;
	model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4;
	model_to_idx["GT"] = 5;
	string FILTER;
	vector<char> variant_line;
	entry *e = get_entry_object();

	map<string, pair<int, int> > FILTER_to_TsTv;
	map<string, int > FILTER_to_Nsites;
	map<string, int >::iterator FILTER_to_Nsites_it;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true, true);

		string model = e->get_REF() + e->get_ALT_allele(0);
		sort(model.begin(), model.end());

		FILTER = e->get_FILTER();
		FILTER_to_Nsites[FILTER]++;
		if (model_to_idx.find(model) != model_to_idx.end())
		{
			switch (model_to_idx[model])
			{
			case 1:
			case 4:
				FILTER_to_TsTv[FILTER].first++;
				break;
			case 0:
			case 2:
			case 3:
			case 5:
				FILTER_to_TsTv[FILTER].second++;
				break;
			default:
				// Don't count this snp towards Ts/Tv
				break;
			}
		}
	}

	vector<pair<int, string > > count_to_FILTER;
	for ( FILTER_to_Nsites_it=FILTER_to_Nsites.begin() ; FILTER_to_Nsites_it != FILTER_to_Nsites.end(); ++FILTER_to_Nsites_it )
	{
		FILTER = (*FILTER_to_Nsites_it).first;
		int Nsites = (*FILTER_to_Nsites_it).second;

		count_to_FILTER.push_back(make_pair(Nsites, FILTER));
	}

	sort(count_to_FILTER.begin(), count_to_FILTER.end());

	string output_file = params.output_prefix + ".FILTER.summary";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Filter Summary Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv" << endl;

	for (int i=count_to_FILTER.size()-1; i > -1; i--)
	{
		FILTER = count_to_FILTER[i].second;
		int Ts = FILTER_to_TsTv[FILTER].first;
		int Tv = FILTER_to_TsTv[FILTER].second;
		int Nsites = FILTER_to_Nsites[FILTER];
		out << FILTER << "\t" << Nsites << "\t";
		out << Ts << "\t" << Tv << "\t" << double(Ts)/Tv << endl;
	}
	delete e;
}

void variant_file::output_TsTv(const parameters &params)
{
	// Output Ts/Tv ratios in bins of a given size.
	int bin_size = params.output_TsTv_bin_size;
	LOG.printLOG("Outputting Ts/Tv in bins of " + header::int2str(bin_size) + "bp\n");

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0;
	model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2;
	model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4;
	model_to_idx["GT"] = 5;

	map<string, int> max_pos;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object();
	map<string, vector<int> > Ts_counts;
	map<string, vector<int> > Tv_counts;
	vector<string> chrs;
	string prev_chr = "";

	vector<int> model_counts(6,0);
	double C = 1.0 / double(bin_size);
	unsigned int idx;
	string model;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (!e->is_biallelic_SNP())
			continue;

		model = e->get_REF() + e->get_ALT_allele(0);
		sort(model.begin(), model.end());

		CHROM = e->get_CHROM();
		idx = (unsigned int)(e->get_POS() * C);

		if(idx>=Ts_counts[CHROM].size())
			Ts_counts[CHROM].resize(idx+1,0);
		if(idx>=Tv_counts[CHROM].size())
			Tv_counts[CHROM].resize(idx+1,0);
		if(CHROM != prev_chr)
		{
			chrs.push_back(CHROM);
			prev_chr = CHROM;
		}

		if (model_to_idx.find(model) != model_to_idx.end())
		{
			model_counts[model_to_idx[model]]++;
			switch (model_to_idx[model])
			{
			case 1:
			case 4:
				Ts_counts[CHROM][idx]++;
				break;
			case 0:
			case 2:
			case 3:
			case 5:
				Tv_counts[CHROM][idx]++;
				break;
			default:
				LOG.error("Unknown idx\n");
				break;
			}
		}
		else
			LOG.warning("Unknown model type. Not a SNP? " + CHROM + ":" + header::int2str(e->get_POS()) +"\n");
	}

	string output_file = params.output_prefix + ".TsTv";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open TsTv Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tBinStart\tSNP_count\tTs/Tv" << endl;
	double ratio;
	for(unsigned int ui=0; ui<chrs.size(); ui++)
	{
		CHROM = chrs[ui];
		for (unsigned int s=0; s<Ts_counts[CHROM].size(); s++)
		{
			ratio = 0.0;
			if (Tv_counts[CHROM][s] != 0)
				ratio = double(Ts_counts[CHROM][s]) / Tv_counts[CHROM][s];
			out << CHROM << "\t" << s*bin_size << "\t" << Ts_counts[CHROM][s]+Tv_counts[CHROM][s] << "\t" << ratio << endl;
		}
	}
	unsigned int Ts = model_counts[1] + model_counts[4];
	unsigned int Tv = model_counts[0] + model_counts[2] + model_counts[3] + model_counts[5];

	LOG.printLOG("Ts/Tv ratio: " + output_log::dbl2str(double(Ts)/Tv, 4) + "\n");
	delete e;
}

void variant_file::output_TsTv_summary(const parameters &params)
{
	// Output Ts/Tv summary.
	LOG.printLOG("Outputting Ts/Tv summary\n");

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0; model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2; model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4; model_to_idx["GT"] = 5;

	vector<char> variant_line;
	entry *e = get_entry_object();

	vector<unsigned int> model_counts(6,0);
	string model;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (!e->is_biallelic_SNP())
			continue;

		model = e->get_REF() + e->get_ALT_allele(0);
		sort(model.begin(), model.end());

		if (model_to_idx.find(model) != model_to_idx.end())
			model_counts[model_to_idx[model]]++;
		else
			LOG.warning("Unknown model type. Not a SNP? " + e->get_CHROM() + ":" + header::int2str(e->get_POS()) +"\n");
	}

	string output_file = params.output_prefix + ".TsTv.summary";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open TsTv Summary Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "MODEL\tCOUNT" << endl;
	out << "AC\t" << model_counts[0] << endl;
	out << "AG\t" << model_counts[1] << endl;
	out << "AT\t" << model_counts[2] << endl;
	out << "CG\t" << model_counts[3] << endl;
	out << "CT\t" << model_counts[4] << endl;
	out << "GT\t" << model_counts[5] << endl;
	unsigned int Ts = model_counts[1] + model_counts[4];
	unsigned int Tv = model_counts[0] + model_counts[2] + model_counts[3] + model_counts[5];
	out << "Ts\t" << Ts << endl;
	out << "Tv\t" << Tv << endl;

	LOG.printLOG("Ts/Tv ratio: " + output_log::dbl2str(double(Ts)/Tv, 4) + "\n");
	delete e;
}

void variant_file::output_TsTv_by_count(const parameters &params)
{
	// Output Ts/Tv ratios in bins of a given size.
	LOG.printLOG("Outputting Ts/Tv by Alternative Allele Count\n");
	vector<unsigned int> Ts_counts, Tv_counts;
	unsigned int N_kept_indv = N_kept_individuals();
	Ts_counts.resize(2*N_kept_indv);
	Tv_counts.resize(2*N_kept_indv);

	string model;
	vector<char> variant_line;
	entry *e = get_entry_object();
	map<string, unsigned int> model_to_Ts_or_Tv;
	model_to_Ts_or_Tv["AC"] = 1;	model_to_Ts_or_Tv["CA"] = 1;
	model_to_Ts_or_Tv["AG"] = 0;	// Ts
	model_to_Ts_or_Tv["GA"] = 0;	// Ts
	model_to_Ts_or_Tv["AT"] = 1;	model_to_Ts_or_Tv["TA"] = 1;
	model_to_Ts_or_Tv["CG"] = 1;	model_to_Ts_or_Tv["GC"] = 1;
	model_to_Ts_or_Tv["CT"] = 0;	// Ts
	model_to_Ts_or_Tv["TC"] = 0;	// Ts
	model_to_Ts_or_Tv["GT"] = 1;	model_to_Ts_or_Tv["TG"] = 1;
	unsigned int idx;
	vector<int> allele_counts;
	unsigned int allele_count;
	unsigned int N_included_indv;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (!e->is_biallelic_SNP())
			continue;

		e->parse_genotype_entries(true);
		e->get_allele_counts(allele_counts, N_included_indv);
		allele_count = allele_counts[1];

		model = e->get_REF() + e->get_ALT_allele(0);
		if (model_to_Ts_or_Tv.find(model) != model_to_Ts_or_Tv.end())
		{
				idx = model_to_Ts_or_Tv[model];
				if (idx == 0) // Ts
					Ts_counts[allele_count]++;
				else if (idx == 1) // Tv;
					Tv_counts[allele_count]++;
				else
					LOG.error("Unknown model type\n");
		}
		else
			LOG.warning("Unknown model type. Not a SNP? " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()) +"\n");
	}

	string output_file = params.output_prefix + ".TsTv.count";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open TsTv by Count Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	double ratio;
	out << "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv" << endl;
	for (unsigned int ui=0; ui<2*N_kept_indv; ui++)
	{
		ratio = double(Ts_counts[ui]) / Tv_counts[ui];
		out << ui << "\t" << Ts_counts[ui] << "\t" << Tv_counts[ui] << "\t" << ratio << endl;
	}
	delete e;
}

void variant_file::output_TsTv_by_quality(const parameters &params)
{
	// Output Ts/Tv ratios in bins of a given size.
	LOG.printLOG("Outputting Ts/Tv By Quality\n");
	map<double, pair<unsigned int, unsigned int> > TsTv_counts;
	double max_qual = -numeric_limits<double>::max(), min_qual=numeric_limits<double>::max();

	string model;
	vector<char> variant_line;
	entry *e = get_entry_object();
	map<string, unsigned int> model_to_Ts_or_Tv;
	model_to_Ts_or_Tv["AC"] = 1;	model_to_Ts_or_Tv["CA"] = 1;
	model_to_Ts_or_Tv["AG"] = 0;	// Ts
	model_to_Ts_or_Tv["GA"] = 0;	// Ts
	model_to_Ts_or_Tv["AT"] = 1;	model_to_Ts_or_Tv["TA"] = 1;
	model_to_Ts_or_Tv["CG"] = 1;	model_to_Ts_or_Tv["GC"] = 1;
	model_to_Ts_or_Tv["CT"] = 0;	// Ts
	model_to_Ts_or_Tv["TC"] = 0;	// Ts
	model_to_Ts_or_Tv["GT"] = 1;	model_to_Ts_or_Tv["TG"] = 1;
	unsigned int idx;
	double QUAL;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);

		if (!e->is_biallelic_SNP())
			continue;

		QUAL = e->get_QUAL();
		if (QUAL > max_qual)
			max_qual = QUAL;
		if (QUAL < min_qual)
			min_qual = QUAL;

		model = e->get_REF() + e->get_ALT_allele(0);;
		if (model_to_Ts_or_Tv.find(model) != model_to_Ts_or_Tv.end())
		{
			idx = model_to_Ts_or_Tv[model];
			if (idx == 0) // Ts
				TsTv_counts[QUAL].first++;
			else if (idx == 1) // Tv;
				TsTv_counts[QUAL].second++;
			else
				LOG.error("Unknown model type\n");
		}
		else
			LOG.warning("Unknown model type. Not a SNP? " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()) +"\n");
	}

	string output_file = params.output_prefix + ".TsTv.qual";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open TsTv by Quality Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	out << "QUAL_THRESHOLD";
	out << "\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD";
	out << "\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD" << endl;

	unsigned int N_TsTv = TsTv_counts.size();

	vector<double> Ts_sum_below(N_TsTv+1, 0.0), Tv_sum_below(N_TsTv+1, 0.0);
	vector<double> QUAL_vector(N_TsTv+1, 0.0);
	QUAL_vector[0] = min_qual;
	QUAL_vector[N_TsTv] = max_qual;
	idx = 1;
	for (map<double, pair<unsigned int, unsigned int> >::iterator it=TsTv_counts.begin(); it != TsTv_counts.end(); ++it)
	{
		QUAL = (it->first);
		double Ts = (it->second).first;
		double Tv = (it->second).second;
		Ts_sum_below[idx] = Ts_sum_below[idx-1]+Ts;
		Tv_sum_below[idx] = Tv_sum_below[idx-1]+Tv;
		QUAL_vector[idx-1] = QUAL;
		idx++;
	}
	QUAL_vector[N_TsTv] = max_qual;

	vector<double> Ts_sum_above(N_TsTv+1, 0.0), Tv_sum_above(N_TsTv+1, 0.0);
	idx = N_TsTv;
	for (map<double, pair<unsigned int, unsigned int> >::reverse_iterator it=TsTv_counts.rbegin(); it != TsTv_counts.rend(); ++it)
	{
		QUAL = (it->first);
		double Ts = (it->second).first;
		double Tv = (it->second).second;
		Ts_sum_above[idx] = Ts_sum_above[idx+1]+Ts;
		Tv_sum_above[idx] = Tv_sum_above[idx+1]+Tv;
		idx--;
	}

	double Ts_sum, Tv_sum, ratio;
	for (unsigned int ui=1; ui<(N_TsTv+1); ui++)
	{
		QUAL = QUAL_vector[ui-1];
		out << QUAL;
		Ts_sum = Ts_sum_below[ui-1]; Tv_sum = Tv_sum_below[ui-1];
		ratio = Ts_sum / Tv_sum;
		out << "\t" << Ts_sum << "\t" << Tv_sum << "\t" << ratio;
		Ts_sum = Ts_sum_above[ui+1]; Tv_sum = Tv_sum_above[ui+1];
		ratio = Ts_sum / Tv_sum;
		out << "\t" << Ts_sum << "\t" << Tv_sum << "\t" << ratio;
		out << endl;
	}
	delete e;
}

void variant_file::output_site_quality(const parameters &params)
{
	// Output per-site quality information.
	LOG.printLOG("Outputting Quality for Each Site\n");
	string output_file = params.output_prefix + ".lqual";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open TsTv by Count Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS\tQUAL" << endl;
	vector<char> variant_line;
	entry *e = get_entry_object();
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e->get_QUAL() << endl;
	}
	delete e;
}

void variant_file::output_site_depth(const parameters &params, bool output_mean)
{
	// Output per-site depth information
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Site Depth Statistics.");

	LOG.printLOG("Outputting Depth for Each Site\n");
	string output_file = params.output_prefix + ".ldepth";
	if (output_mean)
		output_file += ".mean";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Site Depth Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS\t";
	if (output_mean)
		out << "MEAN_DEPTH\tVAR_DEPTH" << endl;
	else
		out << "SUM_DEPTH\tSUMSQ_DEPTH" << endl;

	int depth;
	vector<char> variant_line;
	entry *e = get_entry_object();

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t";

		unsigned int sum=0;
		unsigned int sumsq=0;
		unsigned int n=0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (e->include_genotype[ui] == false)
				continue;

			e->parse_genotype_entry(ui, false, false, true);
			depth = e->get_indv_DEPTH(ui);
			if (depth >= 0)
			{
				sum += depth;
				sumsq += (depth*depth);
				n++;
			}
		}

		if (output_mean)
		{
			double mean = double(sum) / n;
			double var = ((double(sumsq) / n) - (mean*mean)) * double(n) / double(n-1);
			out << mean << "\t" << var << endl;
		}
		else
			out << sum << "\t" << sumsq << endl;
	}
	delete e;
}

void variant_file::output_weir_and_cockerham_fst(const parameters &params)
{	// Implements the bi-allelic version of Weir and Cockerham's Fst
	if (params.weir_fst_populations.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting Weir and Cockerham Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = params.weir_fst_populations.size();
	indvs_in_pops.resize(N_pops, vector<bool>(meta_data.N_indv, false));
	vector<bool> all_indv(meta_data.N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[meta_data.indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(params.weir_fst_populations[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + params.weir_fst_populations[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	string output_file = params.output_prefix + ".weir.fst";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Fst Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST" << endl;

	entry *e = get_entry_object();
	vector<char> variant_line;

	double sum1=0.0, sum2 = 0.0;
	double sum3=0.0, count = 0.0;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		unsigned int N_alleles = e->get_N_alleles();
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tFst: Only using diploid sites.");
			continue;
		}

		vector<unsigned int> N_hom, N_het;
		vector<double> n(N_pops, 0.0);
		vector<vector<double> > p(N_pops, vector<double>(N_alleles,0.0));

		double nbar = 0.0;
		vector<double> pbar(N_alleles, 0.0);
		vector<double> hbar(N_alleles, 0.0);
		vector<double> ssqr(N_alleles, 0.0);
		double sum_nsqr = 0.0;
		double n_sum = 0.0;

		for (unsigned int i=0; i<N_pops; i++)
		{
			e->get_multiple_genotype_counts(indvs_in_pops[i], e->include_genotype, N_hom, N_het);

			for (unsigned int j=0; j<N_alleles; j++)
			{
				n[i] += N_hom[j] + 0.5*N_het[j];
				p[i][j] = N_het[j] + 2*N_hom[j];

				nbar += n[i];
				pbar[j] += p[i][j];
				hbar[j] += N_het[j];
			}
			for (unsigned int j=0; j<N_alleles; j++)
				p[i][j] /= (2.0*n[i]);	// diploid

			sum_nsqr += (n[i] * n[i]);
		}
		n_sum = accumulate(n.begin(),n.end(),0);
		nbar = n_sum / N_pops;

		for (unsigned int j=0; j<N_alleles; j++)
		{
			pbar[j] /= (n_sum*2.0); //diploid
			hbar[j] /= n_sum;
		}

		for (unsigned int j=0; j<N_alleles; j++)
		{
			for (unsigned int i=0; i<N_pops; i++)
				ssqr[j] += (n[i]*(p[i][j] - pbar[j])*(p[i][j] - pbar[j]));
			ssqr[j] /= ((N_pops-1.0)*nbar);
		}
		double nc = (n_sum - (sum_nsqr / n_sum)) / (N_pops - 1.0);

		vector<double> snp_Fst(N_alleles, 0.0);
		vector<double> a(N_alleles, 0.0);
		vector<double> b(N_alleles, 0.0);
		vector<double> c(N_alleles, 0.0);
		double r = double(N_pops);
		double sum_a = 0.0;
		double sum_all = 0.0;

		for(unsigned int j=0; j<N_alleles; j++)
		{
			a[j] = (ssqr[j] - ( pbar[j]*(1.0-pbar[j]) - (((r-1.0)*ssqr[j])/r) - (hbar[j]/4.0) )/(nbar-1.0))*nbar/nc;
			b[j] = (pbar[j]*(1.0-pbar[j]) - (ssqr[j]*(r-1.0)/r) - hbar[j]*( ((2.0*nbar)-1.0) / (4.0*nbar) ))*nbar / (nbar-1.0) ;
			c[j] = hbar[j] / 2.0;
			snp_Fst[j] = a[j]/(a[j]+b[j]+c[j]);

			if ((!std::isnan(a[j])) && (!std::isnan(b[j])) && (!std::isnan(c[j])))
			{
				sum_a += a[j];
				sum_all += (a[j]+b[j]+c[j]);
			}
		}
		double fst = sum_a/sum_all;
		if (!std::isnan(fst))
		{
			sum1 += sum_a;
			sum2 += sum_all;
			sum3 += fst;
			count++;
		}
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << fst << endl;
	}

	double weighted_Fst = sum1 / sum2;
	double mean_Fst = sum3 / count;

	LOG.printLOG("Weir and Cockerham mean Fst estimate: " + output_log::dbl2str(mean_Fst, 5) + "\n");
	LOG.printLOG("Weir and Cockerham weighted Fst estimate: " + output_log::dbl2str(weighted_Fst, 5) + "\n");
	delete e;
}

void variant_file::output_windowed_weir_and_cockerham_fst(const parameters &params)
{
	int fst_window_size = params.fst_window_size;
	int fst_window_step = params.fst_window_step;
	vector<string> indv_files = params.weir_fst_populations;

	if ((fst_window_step <= 0) || (fst_window_step > fst_window_size))
		fst_window_step = fst_window_size;

	if (indv_files.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting Windowed Weir and Cockerham Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = indv_files.size();
	indvs_in_pops.resize(N_pops, vector<bool>(meta_data.N_indv, false));
	vector<bool> all_indv(meta_data.N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[meta_data.indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(indv_files[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + indv_files[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	string CHROM; string last_chr = "";
	vector<string> chrs;
	vector<char> variant_line;
	entry *e = get_entry_object();

	// Calculate number of bins for each chromosome and allocate memory for them.
	// Each bin is a vector with four entries:
	// N_variant_sites: Number of sites in a window that have VCF entries
	// N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
	// N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
	// N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference allele
	const vector< double > empty_vector(4, 0);	// sum1, sum2, sum3, count
	map<string, vector< vector< double > > > bins;
	double sum1=0.0, sum2 = 0.0;
	double sum3=0.0, count = 0.0;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		unsigned int N_alleles = e->get_N_alleles();

		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tFst: Only using diploid sites.");
			continue;
		}

		vector<unsigned int> N_hom, N_het;
		vector<double> n(N_pops, 0.0);
		vector<vector<double> > p(N_pops, vector<double>(N_alleles,0.0));

		double nbar = 0.0;
		vector<double> pbar(N_alleles, 0.0);
		vector<double> hbar(N_alleles, 0.0);
		vector<double> ssqr(N_alleles, 0.0);
		double sum_nsqr = 0.0;
		double n_sum = 0.0;

		for (unsigned int i=0; i<N_pops; i++)
		{
			e->get_multiple_genotype_counts(indvs_in_pops[i], e->include_genotype, N_hom, N_het);

			for (unsigned int j=0; j<N_alleles; j++)
			{
				n[i] += N_hom[j] + 0.5*N_het[j];
				p[i][j] = N_het[j] + 2*N_hom[j];

				nbar += n[i];
				pbar[j] += p[i][j];
				hbar[j] += N_het[j];
			}
			for (unsigned int j=0; j<N_alleles; j++)
				p[i][j] /= (2.0*n[i]);	// diploid

			sum_nsqr += (n[i] * n[i]);
		}
		n_sum = accumulate(n.begin(),n.end(),0);
		nbar = n_sum / N_pops;

		for (unsigned int j=0; j<N_alleles; j++)
		{
			pbar[j] /= (n_sum*2.0); //diploid
			hbar[j] /= n_sum;
		}

		for (unsigned int j=0; j<N_alleles; j++)
		{
			for (unsigned int i=0; i<N_pops; i++)
				ssqr[j] += (n[i]*(p[i][j] - pbar[j])*(p[i][j] - pbar[j]));
			ssqr[j] /= ((N_pops-1.0)*nbar);
		}
		double nc = (n_sum - (sum_nsqr / n_sum)) / (N_pops - 1.0);

		vector<double> snp_Fst(N_alleles, 0.0);
		vector<double> a(N_alleles, 0.0);
		vector<double> b(N_alleles, 0.0);
		vector<double> c(N_alleles, 0.0);
		double r = double(N_pops);
		double sum_a = 0.0;
		double sum_all = 0.0;

		for(unsigned int j=0; j<N_alleles; j++)
		{
			a[j] = (ssqr[j] - ( pbar[j]*(1.0-pbar[j]) - (((r-1.0)*ssqr[j])/r) - (hbar[j]/4.0) )/(nbar-1.0))*nbar/nc;
			b[j] = (pbar[j]*(1.0-pbar[j]) - (ssqr[j]*(r-1.0)/r) - hbar[j]*( ((2.0*nbar)-1.0) / (4.0*nbar) ))*nbar / (nbar-1.0) ;
			c[j] = hbar[j] / 2.0;
			snp_Fst[j] = a[j]/(a[j]+b[j]+c[j]);

			if ((!std::isnan(a[j])) && (!std::isnan(b[j])) && (!std::isnan(c[j])))
			{
				sum_a += a[j];
				sum_all += (a[j]+b[j]+c[j]);
			}
		}
		double fst = sum_a/sum_all;
		if (!std::isnan(fst))
		{
			int pos = (int)e->get_POS();
			CHROM = e->get_CHROM();
			if (CHROM != last_chr)
			{
				chrs.push_back(CHROM);
				last_chr = CHROM;
			}

			int first = (int) ceil((pos - fst_window_size)/double(fst_window_step));
			if (first < 0)
				first = 0;
			int last = (int) ceil(pos/double(fst_window_step));
			for(int idx = first; idx < last; idx++)
			{
				if (idx >= (int)bins[CHROM].size())
					bins[CHROM].resize(idx+1, empty_vector);

				bins[CHROM][idx][0] += sum_a;
				bins[CHROM][idx][1] += sum_all;
				bins[CHROM][idx][2] += fst;
				bins[CHROM][idx][3]++;
			}

			sum1 += sum_a;
			sum2 += sum_all;
			sum3 += fst;
			count++;
		}
	}

	double weighted_Fst = sum1 / sum2;
	double mean_Fst = sum3 / count;

	LOG.printLOG("Weir and Cockerham mean Fst estimate: " + output_log::dbl2str(mean_Fst, 5) + "\n");
	LOG.printLOG("Weir and Cockerham weighted Fst estimate: " + output_log::dbl2str(weighted_Fst, 5) + "\n");

	string output_file = params.output_prefix + ".windowed.weir.fst";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Fst Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST" << endl;
	for (unsigned int ui=0; ui<chrs.size(); ui++)
	{
		CHROM = chrs[ui];
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			if ((bins[CHROM][s][1] != 0) && (!std::isnan(bins[CHROM][s][0])) && (!std::isnan(bins[CHROM][s][1])) && (bins[CHROM][s][3] > 0))
			{
				double weighted_Fst = bins[CHROM][s][0] / bins[CHROM][s][1];
				double mean_Fst = bins[CHROM][s][2] / bins[CHROM][s][3];

				out << CHROM << "\t"
				<< s*fst_window_step + 1 << "\t"
				<< (s*fst_window_step + fst_window_size) << "\t"
				<< bins[CHROM][s][3] << "\t"
				<< weighted_Fst << "\t" << mean_Fst << endl;
			}
		}
	}
	delete e;
}

void variant_file::output_per_site_nucleotide_diversity(const parameters &params)
{
	// Output nucleotide diversity, calculated on a per-site basis.
	// Pi = average number of pairwise differences
	// Assumes a constant distance of 1 between all possible mutations
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Nucleotide Diversity Statistics.");

	LOG.printLOG("Outputting Per-Site Nucleotide Diversity Statistics...\n");
	string output_file = params.output_prefix + ".sites.pi";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Nucleotide Diversity Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS\tPI" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object();
	vector<int> allele_counts;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tsitePi: Only using fully diploid sites.");
			continue;
		}

		unsigned int N_non_missing_chr;
		e->get_allele_counts(allele_counts, N_non_missing_chr);
		unsigned int total_alleles = std::accumulate(allele_counts.begin(), allele_counts.end(), 0);

		unsigned int N_alleles = e->get_N_alleles();
		int mismatches = 0;
		for(unsigned int allele = 0; allele < N_alleles; allele++)
		{
			int other_alleles_count = (total_alleles - allele_counts[allele]);
			mismatches += (allele_counts[allele] * other_alleles_count);
		}

		int pairs = (total_alleles * (total_alleles - 1));
		double pi = (mismatches/static_cast<double>(pairs));

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << pi << endl;
	}
	delete e;
}

//Output Tajima's D
//Carlson et al. Genome Res (2005)
void variant_file::output_Tajima_D(const parameters &params)
{
	int window_size = params.output_Tajima_D_bin_size;

	if (window_size <= 0)
		return;

	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Tajima's D Statistic.");

	LOG.printLOG("Outputting Tajima's D Statistic...\n");
	string output_file = params.output_prefix + ".Tajima.D";

	double a1=0.0, a2=0.0, b1, b2, c1, c2, e1, e2;
	unsigned int n = N_kept_individuals()*2;
	if (n < 2)
		LOG.error("Require at least two chromosomes!");

	for (unsigned int ui=1; ui<n; ui++)
	{
		a1 += 1.0 / double(ui);
		a2 += 1.0 / double(ui * ui);
	}
	b1 = double(n+1) / 3.0 / double(n-1);
	b2 = 2.0 * double(n*n + n + 3) / 9.0 / double(n) / double(n-1);
	c1 = b1 - (1.0 / a1);
	c2 = b2 - (double(n+2)/double(a1*n)) + (a2/a1/a1);
	e1 = c1 / a1;
	e2 = c2 / ((a1*a1) + a2);

	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object();
	map<string, vector< pair<int, double> > > bins;

	unsigned int idx;
	double C = 1.0 / double(window_size);
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	string prev_chr = "";
	vector<string> chrs;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tTajimaD: Only using bialleleic sites.");
			continue;
		}

		CHROM = e->get_CHROM();
		idx = (unsigned int)(e->get_POS() * C);
		e->parse_genotype_entries(true);

		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tTajimaD: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		double p = double(allele_counts[0]) / N_non_missing_chr;

		if(idx>=bins[CHROM].size())
			bins[CHROM].resize(idx+1, make_pair(0,0));
		if(CHROM != prev_chr)
		{
			chrs.push_back(CHROM);
			prev_chr = CHROM;
		}

		if ((p > 0.0) && (p < 1.0))
		{
			bins[CHROM][idx].first++;
			bins[CHROM][idx].second += p * (1.0-p);
		}
	}
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Tajima D Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tBIN_START\tN_SNPS\tTajimaD" << endl;

	for (unsigned int ui=0; ui<chrs.size(); ui++)
	{
		CHROM = chrs[ui];
		bool output = false;
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			int S = bins[CHROM][s].first;
			double D = numeric_limits<double>::quiet_NaN();
			if (S > 0)
			{
				double pi = 2.0*bins[CHROM][s].second*n/double(n-1);
				double tw = double(S) / a1;
				double var = (e1*S) + e2*S*(S-1);
				D = (pi - tw) / sqrt(var);
				output = true;
			}
			if (output == true)
				out << CHROM << "\t" << s*window_size << "\t" << bins[CHROM][s].first << "\t" << D << endl;
		}
	}
	delete e;
}

void variant_file::output_windowed_nucleotide_diversity(const parameters &params)
{
	// Output nucleotide diversity, as calculated in windows.
	// Average number of pairwise differences in windows.
	int window_size = params.pi_window_size;
	int window_step = params.pi_window_step;

	if (window_size <= 0)
		return;

	if ((window_step <= 0) || (window_step > window_size))
		window_step = window_size;

	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Nucleotide Diversity Statistics.");

	LOG.printLOG("Outputting Windowed Nucleotide Diversity Statistics...\n");
	string output_file = params.output_prefix + ".windowed.pi";

	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object();

	// Calculate number of bins for each chromosome and allocate memory for them.
	// Each bin is a vector with four entries:
	// N_variant_sites: Number of sites in a window that have VCF entries
	// N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
	// N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
	// N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference allele
	const unsigned int N_variant_sites = 0;
	const unsigned int N_variant_site_pairs = 1;
	const unsigned int N_mismatches = 2;
	const unsigned int N_polymorphic_sites = 3;
	const vector< unsigned long > empty_vector(4, 0);
	map<string, vector< vector< unsigned long> > > bins;
	vector<string> chrs;
	string prev_chr;

	// Count polymorphic sites and pairwise mismatches
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned long N_comparisons;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);
		CHROM = e->get_CHROM();

		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\twindowPi: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);

		unsigned int N_site_mismatches = 0;
		for (vector<int>::iterator ac = allele_counts.begin(); ac != allele_counts.end(); ++ac)
		{
			N_site_mismatches += (*ac * (N_non_missing_chr - *ac));
		}

		if (N_site_mismatches == 0)
				continue;	// Site is actually fixed.

		// Place the counts into bins
		int pos = (int)e->get_POS();
		int first = (int) ceil((pos - window_size)/double(window_step));
		if (first < 0)
			first = 0;
		int last = (int) ceil(pos/double(window_step));
		N_comparisons = N_non_missing_chr * (N_non_missing_chr - 1);

		if(CHROM != prev_chr)
		{
			chrs.push_back(CHROM);
			prev_chr = CHROM;
			bins[CHROM].resize(1,empty_vector);
		}
		if(last>= (int)bins[CHROM].size())
			bins[CHROM].resize(last+1,empty_vector);

		for(int idx = first; idx < last; idx++)
		{
			bins[CHROM][idx][N_variant_sites]++;
			bins[CHROM][idx][N_variant_site_pairs] += N_comparisons;
			bins[CHROM][idx][N_mismatches] += N_site_mismatches;

			if(allele_counts[0] < (signed)N_non_missing_chr)
				bins[CHROM][idx][N_polymorphic_sites]++;
		}
	}

	// Calculate and print nucleotide diversity statistics
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Windowed Nucleotide Diversity Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" << endl;

	unsigned long N_monomorphic_sites = 0;
	int N_kept_chr = 2*N_kept_individuals();
	N_comparisons = (N_kept_chr * (N_kept_chr - 1)); 	// Number of pairwise comparisons at a monomorphic site
	unsigned long N_pairs = 0; 								// Number of pairwise comparisons within a window
	double pi = 0;

	for (unsigned int ui=0; ui<chrs.size(); ui++)
	{
		CHROM = chrs[ui];
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			if( (bins[CHROM][s][N_polymorphic_sites] > 0) || (bins[CHROM][s][N_mismatches] > 0) )
			{
				// This number can be slightly off for the last bin since the
				// window size can go off the end of the chromosome.
				N_monomorphic_sites = window_size - bins[CHROM][s][N_variant_sites];

				// The total number of possible pairwise comparisons is the sum of
				// pairwise comparisons at polymorphic sites and pairwise
				// comparisons at monomorphic sites.
				N_pairs = bins[CHROM][s][N_variant_site_pairs] + (N_monomorphic_sites * N_comparisons);

				pi = bins[CHROM][s][N_mismatches] / double(N_pairs);
				out << CHROM << "\t"
				    << s*window_step + 1 << "\t"
				    << (s*window_step + window_size) << "\t"
				    << bins[CHROM][s][N_polymorphic_sites] << "\t"
				    << pi << endl;
			}
		}
	}
	delete e;
}

void variant_file::output_kept_sites(const parameters &params)
{
	// Output lists of sites that have been filtered (or not).
	LOG.printLOG("Outputting Kept Sites...\n");
	string output_file = params.output_prefix + ".kept.sites";

	string CHROM;
	vector<char> variant_line;
	int POS;
	entry *e = get_entry_object();

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Kept Site Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS" << endl;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();
		POS = e->get_POS();
		CHROM = e->get_CHROM();
		out << CHROM << "\t" << POS << endl;
	}
	delete e;
}

void variant_file::output_removed_sites(const parameters &params)
{
	// Output lists of sites that have been filtered (or not).
	LOG.printLOG("Outputting Removed Sites...\n");
	string output_file = params.output_prefix + ".removed.sites";

	string CHROM;
	vector<char> variant_line;
	int POS;
	entry *e = get_entry_object();

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Removed Site Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS" << endl;

	while (!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry();
		POS = e->get_POS();
		CHROM = e->get_CHROM();

		if(!eof())
			out << CHROM << "\t" << POS << endl;
	}
	delete e;
}

void variant_file::output_LROH(const parameters &params)
{
	// Detect and output Long Runs of Homozygosity, following the method
	// developed by Adam Boyko, and described in Auton et al., Genome Research, 2009
	// (Although using Forward-backwards algorithm in place of Viterbi).
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LROH.");

	LOG.printLOG("Outputting Long Runs of Homozygosity (Experimental)... \n");
	string output_file = params.output_prefix + ".LROH";

	unsigned int nGen=4;				// Number of generations since common ancestry
	double genotype_error_rate = 0.01;	// Assumed genotype error rate
	double p_auto_prior = 0.05;			// Prior probability of being in autozygous state
	double p_auto_threshold = 0.99;		// Threshold for reporting autozygous region
	int min_SNPs=0;						// Threshold for reporting autozygous region

	string CHROM;
	vector<char> variant_line;
	int POS;
	entry *e = get_entry_object();
	pair<int, int> alleles;
	vector< vector< int > > s_vector;
	vector< vector<pair<double, double> > > p_emission;
	vector< vector< vector<double> > > p_trans;
	vector<int> last_POS;
	vector<vector<bool> > is_het;

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open LROH Output file: " + output_file, 12);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV" << endl;

	s_vector.resize(meta_data.N_indv);
	p_emission.resize(meta_data.N_indv);
	p_trans.resize(meta_data.N_indv);
	last_POS.resize(meta_data.N_indv,-1);
	is_het.resize(meta_data.N_indv);

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		//if (e->get_N_alleles() != 2)
		//{
		//	LOG.one_off_warning("\tLROH: Only using bialleleic sites.");
		//	continue;	// TODO: Probably could do without this...
		//}

		CHROM = e->get_CHROM();
		POS = e->get_POS();
		double r = 0;

		unsigned int N_genotypes = 0;
		unsigned int N_hets = 0;

		vector<int> indv_alleles(meta_data.N_indv, -1);

		bool has_non_ref = false;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if ((include_indv[ui] == false) || (e->include_genotype[ui] == false))
				continue;

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);

			if (e->get_indv_ploidy(ui) != 2)
			{
				LOG.one_off_warning("\tLROH: Only using diploid sites.");
				continue;
			}

			if ((alleles.first < 0) || (alleles.second < 0))
				continue;	// Skip missing genotypes

			if ((alleles.first > 0) || (alleles.second > 0))
				has_non_ref = true;

			N_genotypes++;
			bool is_het = (alleles.first != alleles.second);
			if (is_het == true)
				N_hets++;

			indv_alleles[ui] = (int)is_het;
		}

		if (has_non_ref == false)
			continue;

		double h = N_hets / double(N_genotypes);	// Heterozygosity
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if ((include_indv[ui] == false) || (e->include_genotype[ui] == false))
				continue;

			double p_emission_given_nonauto;
			double p_emission_given_auto;

			if (indv_alleles[ui] < 0)
				continue;
			else if (indv_alleles[ui] == 1)
			{	// Heterozygote
				p_emission_given_nonauto = h;
				p_emission_given_auto = genotype_error_rate;
				p_emission[ui].push_back(make_pair(p_emission_given_auto, p_emission_given_nonauto));
				is_het[ui].push_back(true);
			}
			else
			{	// Homozygote
				p_emission_given_nonauto = 1.0-h;
				p_emission_given_auto = 1.0-genotype_error_rate;
				p_emission[ui].push_back(make_pair(p_emission_given_auto, p_emission_given_nonauto));
				is_het[ui].push_back(false);
			}

			if (last_POS[ui] > 0)
			{	// Assume 1cM/Mb.
				r = (POS - last_POS[ui]) / 1000000.0 / 100.0; // Morgans
			}

			double e = (1.0 - exp(-2.0*nGen*r));
			double p_trans_auto_to_nonauto = (1.0 - p_auto_prior) * e;	//A[1]
			double p_trans_nonauto_to_auto = p_auto_prior * e;	//A[2]
			double p_trans_auto_to_auto = 1.0 - p_trans_nonauto_to_auto;	//A[0]
			double p_trans_nonauto_to_nonauto = 1.0 - p_trans_auto_to_nonauto; // A[3]
			vector<double> A(4);
			A[0] = p_trans_auto_to_auto;
			A[1] = p_trans_auto_to_nonauto;
			A[2] = p_trans_nonauto_to_auto;
			A[3] = p_trans_nonauto_to_nonauto;

			p_trans[ui].push_back(A);
			s_vector[ui].push_back(POS);
			last_POS[ui] = POS;
		}
	}
	delete e;

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		LOG.printLOG("\t" + meta_data.indv[ui] + "\n");

		// Forward-backward algorithm
		int N_obs = (int)p_emission[ui].size();
		if (N_obs == 0)
			continue;

		vector<vector<double> > alpha(N_obs, vector<double>(2,0));
		vector<vector<double> > beta(N_obs, vector<double>(2,0));

		alpha[0][0] = p_emission[ui][0].first;
		alpha[0][1] = p_emission[ui][0].second;

		for (int i=1; i<N_obs; i++)
		{
			alpha[i][0] = alpha[i-1][0] * p_trans[ui][i-1][0] * p_emission[ui][i].first;
			alpha[i][0] += alpha[i-1][1] * p_trans[ui][i-1][2] * p_emission[ui][i].first;

			alpha[i][1] = alpha[i-1][1] * p_trans[ui][i-1][3] * p_emission[ui][i].second;
			alpha[i][1] += alpha[i-1][0] * p_trans[ui][i-1][1] * p_emission[ui][i].second;

			while (alpha[i][0] + alpha[i][1] < 1e-20)
			{	// Renormalise to prevent underflow
				alpha[i][0] *= 1e20;
				alpha[i][1] *= 1e20;
			}
		}
		beta[N_obs-1][0] = 1.0;
		beta[N_obs-1][1] = 1.0;
		for (int i=N_obs-2; i>=0; i--)
		{
			beta[i][0] = beta[i+1][0] * p_trans[ui][i][0] * p_emission[ui][i].first;
			beta[i][0] += beta[i+1][1] * p_trans[ui][i][2] * p_emission[ui][i].first;

			beta[i][1] = beta[i+1][1] * p_trans[ui][i][3] * p_emission[ui][i].second;
			beta[i][1] += beta[i+1][0] * p_trans[ui][i][1] * p_emission[ui][i].second;

			while (beta[i][0] + beta[i][1] < 1e-20)
			{	// Renormalise to prevent underflow
				beta[i][0] *= 1e20;
				beta[i][1] *= 1e20;
			}
		}
		// Calculate probability of each site being autozygous
		vector<double> p_auto(N_obs);
		for (int i=0; i<N_obs; i++)
			p_auto[i] = alpha[i][0] * beta[i][0] / (alpha[i][0] * beta[i][0] + alpha[i][1] * beta[i][1]);

		// Generate output
		bool in_auto=false;
		int start_pos=0, end_pos=0;
		int N_SNPs = 0;
		int N_SNPs_between_hets = 0;
		int N_hets_in_region = 0;
		int last_het_pos = s_vector[ui][0];
		int next_het_pos = -1;
		for (int i=0; i<N_obs; i++)
		{
			if (p_auto[i] > p_auto_threshold)
			{
				if (in_auto == false)
				{	// Start of autozygous region
					start_pos = s_vector[ui][i];
				}
				N_SNPs++;
				N_SNPs_between_hets++;
				if (is_het[ui][i] == true)
					N_hets_in_region++;
				in_auto = true;
			}
			else
			{
				if (in_auto == true)
				{	// end of autozygous region
					// Find next_het position
					next_het_pos = s_vector[ui][N_obs-1];
					for (int j=i; j<N_obs; j++)
					{
						if (is_het[ui][j] == true)
						{
							next_het_pos = s_vector[ui][j];
							break;
						}
						N_SNPs_between_hets++;
					}

					end_pos = s_vector[ui][i-1];
					if (N_SNPs >= min_SNPs)
					{
						out << CHROM << "\t" << start_pos << "\t" << end_pos << "\t" << (last_het_pos+1) << "\t" << (next_het_pos-1) << "\t" << N_SNPs_between_hets << "\t" << N_hets_in_region << "\t" << meta_data.indv[ui] << endl;
					}
				}
				in_auto = false;
				N_SNPs = 0;
				N_hets_in_region = 0;
				if (is_het[ui][i] == true)
				{
					last_het_pos = s_vector[ui][i];
					N_SNPs_between_hets = 0;
				}
			}
		}
		if (in_auto == true)
		{	// Report final region if needed
			end_pos = s_vector[ui][N_obs-1];
			next_het_pos = s_vector[ui][N_obs-1];
			if (N_SNPs >= min_SNPs)
				out << CHROM << "\t" << start_pos << "\t" << end_pos << "\t" << (last_het_pos+1) << "\t" << next_het_pos << "\t" << N_SNPs_between_hets << "\t" << N_hets_in_region << "\t" << meta_data.indv[ui] << endl;
		}
	}
}

void variant_file::output_indv_relatedness_Manichaikul(const parameters &params)
{
	// Calculate and output a relatedness statistic based on the method of
	// Manichaikul et al., BIOINFORMATICS 2010
	// doi:10.1093/bioinformatics/btq559
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Individual Relatedness.");

	LOG.printLOG("Outputting Individual Relatedness\n");
	string output_file = params.output_prefix + ".relatedness2";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Individual Relatedness Output file: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object();
	vector<int> allele_counts;
	unsigned int N_alleles;
	pair<int, int> geno_id;
	pair<int, int> geno_id2;
	vector<vector<double> > phi(meta_data.N_indv, vector<double>(meta_data.N_indv, 0.0));
	vector<vector<double> > N_AaAa(meta_data.N_indv, vector<double>(meta_data.N_indv, 0.0));
	vector<vector<double> > N_AAaa(meta_data.N_indv, vector<double>(meta_data.N_indv, 0.0));
	vector<double > N_Aa(meta_data.N_indv, 0.0);

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tRelatedness: Only using biallelic sites.");
			continue;	// Only use biallelic loci
		}

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tRelatedness: Only using fully diploid sites.");
			continue;
		}

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			if ((geno_id.first != geno_id.second) && (geno_id.first >= 0) && (geno_id.second >= 0))
			{
				N_Aa[ui]++;
			}

			for (unsigned int uj=0; uj<meta_data.N_indv; uj++)
			{
				if (include_indv[uj] == false)
					continue;

				e->get_indv_GENOTYPE_ids(uj, geno_id2);
				if ((geno_id.first != geno_id.second) && (geno_id.first >= 0) && (geno_id.second >= 0))
				{
					if ((geno_id2.first != geno_id2.second) && (geno_id2.first >= 0) && (geno_id2.second >= 0))
					{
						N_AaAa[ui][uj]++;
					}
				}
				if ((geno_id.first == geno_id.second) && (geno_id.first >= 0) && (geno_id.second >= 0))
				{
					if ((geno_id2.first == geno_id2.second) && (geno_id2.first >= 0) && (geno_id2.second >= 0))
					{
						if (geno_id.first != geno_id2.first)
						{
							N_AAaa[ui][uj]++;
						}
					}
				}
			}
		}
	}

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		for (unsigned int uj=0; uj<meta_data.N_indv; uj++)
		{
			if (include_indv[uj] == false)
				continue;
			phi[ui][uj] = (N_AaAa[ui][uj] - 2.0*N_AAaa[ui][uj]) / (N_Aa[ui] + N_Aa[uj]);
			out << meta_data.indv[ui] << "\t" << meta_data.indv[uj];
			out << "\t" << N_AaAa[ui][uj] << "\t" << N_AAaa[ui][uj] << "\t" << N_Aa[ui] << "\t" << N_Aa[uj];
			out << "\t" << phi[ui][uj] << endl;
		}
	}
	delete e;
}

void variant_file::output_indv_relatedness_Yang(const parameters &params)
{
	// Calculate and output a relatedness statistic based on the method of
	// Yang et al, 2010 (doi:10.1038/ng.608). Specifically, calculate the
	// unadjusted Ajk statistic (equation 6 of paper).
	// Expectation of Ajk is zero for individuals within a populations, and
	// one for an individual with themselves.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Individual Relatedness.");

	LOG.printLOG("Outputting Individual Relatedness\n");
	string output_file = params.output_prefix + ".relatedness";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Individual Relatedness Output file: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "INDV1\tINDV2\tRELATEDNESS_AJK" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object();
	vector<int> allele_counts;
	unsigned int N_alleles, N_non_missing_chr;
	double freq;
	pair<int, int> geno_id;
	vector<vector<double> > Ajk(meta_data.N_indv, vector<double>(meta_data.N_indv, 0.0));
	vector<vector<double> > N_sites(meta_data.N_indv, vector<double>(meta_data.N_indv, 0.0));

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tRelatedness: Only using biallelic sites.");
			continue;	// Only use biallelic loci
		}

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tRelatedness: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

		if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
			continue;

		vector<double> x(meta_data.N_indv, -1.0);
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			x[ui] = geno_id.first + geno_id.second;
		}

		double div = 1.0/(2.0*freq*(1.0-freq));
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if ((include_indv[ui] == false) || (e->include_genotype[ui] == false) || (x[ui] < 0))
				continue;
			Ajk[ui][ui] += (x[ui]*x[ui] - (1 + 2.0*freq)*x[ui] + 2.0*freq*freq) * div;
			N_sites[ui][ui]++;
			for (unsigned int uj=(ui+1); uj<meta_data.N_indv; uj++)
			{
				if ((include_indv[uj] == false) || (e->include_genotype[uj] == false) || (x[uj] < 0))
					continue;
				Ajk[ui][uj] += (x[ui] - 2.0*freq) * (x[uj] - 2.0*freq) * div;
				N_sites[ui][uj]++;
			}
		}
	}

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		Ajk[ui][ui] = 1.0 + (Ajk[ui][ui] / N_sites[ui][ui]);
		out << meta_data.indv[ui] << "\t" << meta_data.indv[ui] << "\t" << Ajk[ui][ui] << endl;
		for (unsigned int uj=(ui+1); uj<meta_data.N_indv; uj++)
		{
			if (include_indv[uj] == false)
				continue;
			Ajk[ui][uj] /= N_sites[ui][uj];
			out << meta_data.indv[ui] << "\t" << meta_data.indv[uj] << "\t" << Ajk[ui][uj] << endl;
		}
	}
	delete e;
}

void variant_file::output_PCA(const parameters &params)
{
#ifndef VCFTOOLS_PCA
	string out = "Cannot run PCA analysis. Vcftools has been compiled without PCA enabled (requires LAPACK).";
	LOG.error(out);
#else
    bool use_normalisation = !params.PCA_no_normalisation;
	// Output PCA, following method of Patterson, Price and Reich 2006.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to perform PCA.");

	if (use_normalisation)
		LOG.printLOG("Outputting Principal Component Analysis (with normalisation)\n");
	else
		LOG.printLOG("Outputting Principal Component Analysis (without normalisation)\n");
	string output_file = params.output_prefix + ".pca";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Principal Component Analysis Output file: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);

	unsigned int N_indvs = N_kept_individuals();
	vector<char> variant_line;
	entry *e = get_entry_object();
	pair<int, int> geno_id;
	double x, freq;
	vector<int> allele_counts;
	unsigned int N_alleles, N_non_missing_chr;

	// Store list of included individuals
	vector<string> included_indvs(N_indvs);
	unsigned int ui_prime = 0;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		included_indvs[ui_prime] = meta_data.indv[ui];
		ui_prime++;
	}

	// Potentially uses a lot of memory. Should issue a warning about this.
	vector< vector<double> > M(N_indvs);

	// Populate M
	unsigned int s_prime = 0;
	unsigned int N_sites = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;

		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();
		if (N_alleles != 2)
			LOG.error("PCA only works for biallelic sites.");

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
			LOG.error("PCA only works for fully diploid sites. Non-diploid site at " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()));

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

		if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
			continue;

		double mu = freq*2.0;
		double div = 1.0 / sqrt(freq * (1.0-freq));

		ui_prime = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			x = geno_id.first + geno_id.second;
			if (x > -1)
			{
				if (use_normalisation == true)
					M[ui_prime].push_back((x - mu) * div);
				else
					M[ui_prime].push_back((x - mu));
			}
			ui_prime++;
		}
		s_prime++;
		N_sites++;
	}

	if (N_indvs >= N_sites)
		LOG.error("PCA computation requires that there are more sites than individuals.");

	// Now construct X = (1/n)MM'.
	double **X = new double *[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		X[ui] = new double[N_indvs];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] = 0;

	// Only populate one half of matrix
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=ui; uj<N_indvs; uj++)
			for (unsigned int s=0; s<N_sites; s++)
				X[ui][uj] += M[ui][s] * M[uj][s];

	// Populate other half
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<ui; uj++)
			X[ui][uj] = X[uj][ui];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] /= N_sites;

	double *Er = new double[N_indvs];
	double *Ei = new double[N_indvs];
	double **Evecs = new double*[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		Evecs[ui] = new double[N_indvs];

	// Call LAPACK routine to calculate eigenvectors and eigenvalues
	dgeev(X, N_indvs, Er, Ei, Evecs);

	// Check there are no complex eigenvalues.
	for (unsigned int ui=0; ui<N_indvs; ui++)
		if (Ei[ui] != 0)
			LOG.error("Complex eigenvalue.");

	// Output results
	out << "INDV";
	for (unsigned int ui=0; ui<N_indvs; ui++)
		out << "\tEIG_" << ui;
	out << endl;

	out << "EIGENVALUE";
	for (unsigned int ui=0; ui<N_indvs; ui++)
		out << "\t" << Er[ui];
	out << endl;

	// Output eigenvectors (as columns)
	for (unsigned int ui=0; ui<N_indvs; ui++)
	{
		out << included_indvs[ui];
		for (unsigned int uj=0; uj<N_indvs; uj++)
			out << "\t" << Evecs[ui][uj];
		out << endl;
	}

	delete e;
	delete [] Er;
	delete [] Ei;
	delete [] Evecs;
	delete [] X;
#endif
}

void variant_file::output_PCA_SNP_loadings(const parameters &params)
{
    // TODO: This function duplicates a lot of what is in the output PCA function. Would be better to combine in a more
    // sensible fashion.
#ifndef VCFTOOLS_PCA
	string out = "Cannot run PCA analysis. Vcftools has been compiled without PCA enabled (requires LAPACK).";
	LOG.error(out);
#else
    int SNP_loadings_N_PCs = params.output_N_PCA_SNP_loadings;
    bool use_normalisation = !params.PCA_no_normalisation;

	// Output PCA, following method of Patterson, Price and Reich 2006.
	if ((meta_data.has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to perform PCA.");

	if (use_normalisation)
		LOG.printLOG("Outputting " + header::int2str(SNP_loadings_N_PCs) + " SNP loadings (with normalisation)\n");
	else
		LOG.printLOG("Outputting " + header::int2str(SNP_loadings_N_PCs) + " SNP loadings (without normalisation)\n");
	string output_file = params.output_prefix + ".pca.loadings";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Principal Component SNP Loading Output file: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS";
	for (unsigned int ui=0; ui<(unsigned int)SNP_loadings_N_PCs; ui++)
		out << "\tGAMMA_" << ui;
	out << endl;

	unsigned int N_indvs = N_kept_individuals();
	vector<string> included_indvs(N_indvs);
	unsigned int ui_prime = 0;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		included_indvs[ui_prime] = meta_data.indv[ui];
		ui_prime++;
	}

	vector< vector<double> > M(N_indvs);
	vector< vector<double> > ids(N_indvs);
    
    map<int, string> idx_to_chrom;
    map<string, int> chrom_to_idx;
    string chr;
    int pos;
    int chrom_idx = 0;
    vector<int> CHROMidx_list;
    vector<int> pos_list;

	vector<char> variant_line;
	entry *e = get_entry_object();
	pair<int, int> geno_id;
	double x, freq;
	vector<int> allele_counts;
	unsigned int N_alleles, N_non_missing_chr;

	ui_prime = 0;
	unsigned int s_prime = 0;
	unsigned int N_sites = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();
		if (N_alleles != 2)
			LOG.error("PCA only works for biallelic sites.");

		e->parse_genotype_entries(true);
		if (e->is_diploid() == false)
			LOG.error("PCA only works for fully diploid sites. Non-diploid site at " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()));

		e->get_allele_counts(allele_counts, N_non_missing_chr);
		freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

		if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
			continue;

		double mu = freq*2.0;
		double div = 1.0 / sqrt(freq * (1.0-freq));
        
        chr = e->get_CHROM();
        pos = e->get_POS();
        
        if (chrom_to_idx.find(chr) == chrom_to_idx.end())
        {
            chrom_to_idx[chr] = chrom_idx;
            idx_to_chrom[chrom_idx] = chr;
            chrom_idx++;
        }
        
        CHROMidx_list.push_back(chrom_to_idx[chr]);
        pos_list.push_back(pos);

		ui_prime = 0;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			x = geno_id.first + geno_id.second;
			if (x > -1)
			{
				if (use_normalisation == true)
					M[ui_prime].push_back((x - mu) * div);
				else
					M[ui_prime].push_back((x - mu));
			}
			ids[ui_prime].push_back(x);
			ui_prime++;
		}
		s_prime++;
		N_sites++;
	}

	if (N_indvs >= N_sites)
		LOG.error("PCA computation requires that there are more sites than individuals.");

	// Now construct X = (1/n)MM'.
	double **X = new double *[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		X[ui] = new double[N_indvs];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] = 0;

	// Only populate one half of matrix
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=ui; uj<N_indvs; uj++)
			for (unsigned int s=0; s<N_sites; s++)
				X[ui][uj] += M[ui][s] * M[uj][s];

	// Populate other half
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<ui; uj++)
			X[ui][uj] = X[uj][ui];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] /= N_sites;

	double *Er = new double[N_indvs];
	double *Ei = new double[N_indvs];
	double **Evecs = new double*[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		Evecs[ui] = new double[N_indvs];

	// Call LAPACK routine to calculate eigenvectors and eigenvalues
	dgeev(X, N_indvs, Er, Ei, Evecs);

	// Check there are no complex eigenvalues.
	for (unsigned int ui=0; ui<N_indvs; ui++)
		if (Ei[ui] != 0)
			LOG.error("Complex eigenvalue.");

	for (unsigned int ui=0; ui<N_sites; ui++)
	{
		vector<double> gamma(SNP_loadings_N_PCs, 0.0);
		vector<double> a_sum(SNP_loadings_N_PCs, 0.0);
        
        chr = idx_to_chrom[CHROMidx_list[ui]];
        pos = pos_list[ui];
        
		out << chr << "\t" << pos;

		for (unsigned int uj_prime=0; uj_prime<N_indvs; uj_prime++)
		{
			x = ids[uj_prime][ui];
			if (x > -1)
			{
				for (unsigned int uk=0; uk<(unsigned int)SNP_loadings_N_PCs; uk++)
				{
					gamma[uk] += (x * Evecs[uj_prime][uk]);
					a_sum[uk] += (Evecs[uj_prime][uk]*Evecs[uj_prime][uk]);
				}
			}
		}
		for (unsigned int uj=0; uj<(unsigned int)SNP_loadings_N_PCs; uj++)
			out << "\t" << gamma[uj] / a_sum[uj];
		out << endl;
	}
	delete e;
    delete [] Er;
    delete [] Ei;
    delete [] Evecs;
    delete [] X;
#endif
}

void variant_file::output_indel_hist(const parameters &params)
{
	vector<char> variant_line;
	entry *e = get_entry_object();
	string allele;
	unsigned int ref_len, N_alleles;
	int indel_len, smallest_len, largest_len, snp_count;
	vector<int> s_vector;

	string output_file = params.output_prefix + ".indel.hist";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Indel Histogram Output file: " + output_file, 2);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	LOG.printLOG("Outputting Indel Histogram\n");
	out << "LENGTH\tCOUNT\tPRCT" << endl;
	largest_len = 0;
	smallest_len = 0;
	snp_count = 0;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		allele = e->get_REF();
		ref_len = allele.size();
		N_alleles = e->get_N_alleles();
		if (e->is_SNP() )
			snp_count++;

		for (unsigned int ui=1; ui<N_alleles; ui++)
		{
			e->get_allele(ui, allele);
			if (allele.size() != ref_len)
			{
				if (allele.find_first_not_of("acgtACGT") == string::npos)
				{	// Check all bases are ATCGatcg
					indel_len = allele.size() - ref_len;
					s_vector.push_back (indel_len);
					if (indel_len > largest_len)
						largest_len = indel_len;
					else if (indel_len < smallest_len)
						smallest_len = indel_len;
				}
			}
		}
	}

	double total = s_vector.size() + snp_count;
	double pct;
	for (int i=smallest_len; i<=largest_len; i++)
	{
		int icount = (int) count (s_vector.begin(), s_vector.end(), i);
		if (icount > 0)
		{
			pct = 100.0*icount/total;
			out << i << "\t" << icount << "\t" << pct << endl;

		}
		else if ((i == 0) and (snp_count>0))
		{
			pct = 100.0*snp_count/total;
			out << i << "\t" << snp_count << "\t" << pct << endl;
		}
	}
}

void variant_file::output_mendel_inconsistencies(const parameters &params)
{
	LOG.printLOG("Outputting Mendel Errors.\n");

	ifstream PED(params.mendel_ped_file.c_str());
	if (!PED.is_open())
		LOG.error("Could not open PED file: " + params.mendel_ped_file);

	string line; stringstream ss;
	string family, child, mother, father;
	vector<int> child_idx, mother_idx, father_idx;
	vector<string> family_ids;
	PED.ignore(numeric_limits<streamsize>::max(), '\n');

	while (!PED.eof())
	{
		getline(PED, line);
		if ((line[0] == '#') || (line.size() == 0))
			continue;

		ss.clear();
		ss.str(line);
		ss >> family >> child >> father >> mother;

		if ((child == "0") || (father == "0") || (mother == "0"))
			continue;

		int idx1 = -1, idx2 = -1, idx3 = -1;
		vector<string>::iterator it = find(meta_data.indv.begin(), meta_data.indv.end(), child);
		if (it != meta_data.indv.end())
			idx1 = distance(meta_data.indv.begin(), it);
		it = find(meta_data.indv.begin(), meta_data.indv.end(), mother);
		if (it != meta_data.indv.end())
			idx2 = distance(meta_data.indv.begin(), it);
		it = find(meta_data.indv.begin(), meta_data.indv.end(), father);
		if (it != meta_data.indv.end())
			idx3 = distance(meta_data.indv.begin(), it);

		if ((idx1 != -1) && (idx2 != -1) && (idx3 != -1))
		{	// Trio is in the VCF
			if (include_indv[idx1] == false)
				continue;
			if (include_indv[idx2] == false)
				continue;
			if (include_indv[idx3] == false)
				continue;

			child_idx.push_back(idx1); mother_idx.push_back(idx2); father_idx.push_back(idx3);
			family_ids.push_back(child + "_" + father + "_" + mother);
		}
	}
	PED.close();

	LOG.printLOG("Found " + LOG.int2str(child_idx.size()) + " trios in the VCF file.\n");

	if (child_idx.size() == 0)
		LOG.error("No PED individuals found in VCF.\n", 5);

	string output_file = params.output_prefix + ".mendel";

	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Mendel Error output file: " + output_file, 4);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);


	pair<int, int> child_alleles, child_alleles2;
	pair<int, int> mother_alleles;
	pair<int, int> father_alleles;
	string CHROM; int POS;
	string REF, ALT;
	out << "CHR\tPOS\tREF\tALT\tFAMILY\tCHILD\tFATHER\tMOTHER" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object();
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);

		CHROM = e->get_CHROM();
		POS = e->get_POS();

		for (unsigned int trio=0; trio<child_idx.size(); trio++)
		{
			int idx1 = child_idx[trio], idx2 = mother_idx[trio], idx3 = father_idx[trio];
			if ((e->include_genotype[idx1] == false) || (e->include_genotype[idx2] == false) || (e->include_genotype[idx3] == false))
				continue;

			e->get_indv_GENOTYPE_ids(idx1, child_alleles);

			e->get_indv_GENOTYPE_ids(idx2, mother_alleles);

			e->get_indv_GENOTYPE_ids(idx3, father_alleles);

			if ((child_alleles.first == -1) || (child_alleles.second == -1) ||
					(mother_alleles.first == -1) || (mother_alleles.second == -1) ||
					(father_alleles.first == -1) || (father_alleles.second == -1))
				continue;

		//	cout << CHROM << "\t" << POS << "\t" << REF << "\t" << ALT << "\t" << family_ids[trio] << "\t" << child_alleles.first << "/" << child_alleles.second;
		//	cout << "\t" << father_alleles.first << "/" << father_alleles.second << "\t" << mother_alleles.first << "/" << mother_alleles.second << endl;

			set<pair<int, int> > possible_child_genotypes;
			possible_child_genotypes.insert(make_pair(mother_alleles.first, father_alleles.first));
			possible_child_genotypes.insert(make_pair(mother_alleles.first, father_alleles.second));
			possible_child_genotypes.insert(make_pair(mother_alleles.second, father_alleles.first));
			possible_child_genotypes.insert(make_pair(mother_alleles.second, father_alleles.second));

			child_alleles2 = make_pair(child_alleles.second, child_alleles.first);

			if ((possible_child_genotypes.find(child_alleles) == possible_child_genotypes.end()) && (possible_child_genotypes.find(child_alleles2) == possible_child_genotypes.end()))
			{	// Mendel error!
				CHROM = e->get_CHROM();
				POS = e->get_POS();
				REF = e->get_REF();
				ALT = e->get_ALT();
				out << CHROM << "\t" << POS << "\t" << REF << "\t" << ALT << "\t" << family_ids[trio] << "\t" << child_alleles.first << "/" << child_alleles.second;
				out << "\t" << father_alleles.first << "/" << father_alleles.second << "\t" << mother_alleles.first << "/" << mother_alleles.second << endl;
			}
		}
	}
	delete e;
}

void variant_file::write_stats(const parameters &params)
{
	vector<char> variant_line;
	entry *e = get_entry_object();

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);
	}
	delete e;
}
