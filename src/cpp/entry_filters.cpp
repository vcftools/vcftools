/*
 * entry_filters.cpp
 *
 *  Created on: Aug 9, 2013
 *      Author: amarcketta
 */

#include "entry.h"

set<string> entry::local_snps_to_keep;
set<string> entry::snps_to_exclude;
vector< set<int > > entry::keep_positions;
vector< set<int > > entry::exclude_positions;
map<string,int> entry::chr_to_idx;
vector< deque<pair<int,int> > > entry::lims;
ifstream entry::mask;
string entry::mask_chr;
string entry::mask_line;
int entry::mask_pos;
string entry::thin_chrom;
int entry::thin_pos;

int entry::apply_filters(const parameters &params)
{
	if (line.empty())
	{
		passed_filters = false;
		return 0;
	}

// 	Apply all filters in turn.
	filter_sites_by_allele_type(params.keep_only_indels, params.remove_indels);
	filter_sites(params.snps_to_keep, params.snps_to_keep_file, params.snps_to_exclude_file);
	filter_sites_by_filter_status(params.site_filter_flags_to_exclude, params.site_filter_flags_to_keep, params.remove_all_filtered_sites);
	string chr_to_keep = "";
	if (params.chrs_to_keep.size() == 1)
		chr_to_keep = *(params.chrs_to_keep.begin()); // Get first chromosome in list (there should only be one).
	filter_sites_by_position(chr_to_keep, params.start_pos, params.end_pos);
	filter_sites_by_positions(params.positions_file, params.exclude_positions_file);
	filter_sites_by_overlap_positions(params.positions_overlap_file, params.exclude_positions_overlap_file);
	filter_sites_by_chromosome(params.chrs_to_keep, params.chrs_to_exclude);
	filter_sites_by_BED_file(params.BED_file, params.BED_exclude);
	filter_sites_by_number_of_alleles(params.min_alleles, params.max_alleles);
	filter_sites_by_INFO(params.site_INFO_flags_to_remove, params.site_INFO_flags_to_keep);
	filter_sites_by_quality(params.min_quality);
	filter_sites_by_mean_depth(params.min_mean_depth, params.max_mean_depth);
	filter_sites_by_mask(params.mask_file, params.invert_mask, params.min_kept_mask_value);
	if (params.phased_only == true)
		filter_sites_by_phase();
	filter_genotypes_by_quality_value(params.min_genotype_quality);
	filter_genotypes_by_depth_range(params.min_genotype_depth, params.max_genotype_depth);
	filter_genotypes_by_filter_flag(params.geno_filter_flags_to_exclude, params.remove_all_filtered_genotypes);
	filter_sites_by_frequency_and_call_rate(params.min_maf, params.max_maf, params.min_non_ref_af, params.max_non_ref_af, params.min_non_ref_af_any, params.max_non_ref_af_any, params.min_site_call_rate);
	filter_sites_by_allele_count(params.min_mac, params.max_mac, params.min_non_ref_ac, params.max_non_ref_ac, params.min_non_ref_ac_any, params.max_non_ref_ac_any, params.max_missing_call_count);
	filter_sites_by_HWE_pvalue(params.min_HWE_pvalue);
	filter_sites_by_thinning(params.min_interSNP_distance);

	return 1;
}

void entry::filter_genotypes_by_quality_value(double min_genotype_quality)
{
	// Filter genotypes by quality
	if (passed_filters == false)
		return;

	if (min_genotype_quality <= 0)
		return;

	parse_genotype_entries(false, true);

	if (entry_header.has_genotypes == false)
		LOG.error("Require Genotypes in variant file in order to filter genotypes by Quality.");

	filter_genotypes_by_quality(min_genotype_quality);
}

void entry::filter_genotypes_by_depth_range(int min_depth, int max_depth)
{
	// Filter genotypes by depth
	if (passed_filters == false)
		return;

	if ((min_depth <= 0) && (max_depth == numeric_limits<int>::max()))
		return;
	if (entry_header.has_genotypes == false)
		LOG.error("Require Genotypes in variant file in order to filter genotypes by Depth.");

	parse_genotype_entries(false, false, true);
	filter_genotypes_by_depth(min_depth, max_depth);

}

void entry::filter_genotypes_by_filter_flag(const set<string> &filter_flags_to_remove, bool remove_all)
{
	// Filter genotypes by Filter Flags
	if (passed_filters == false)
		return;

	if ((remove_all == false) && (filter_flags_to_remove.empty()))
		return;

	parse_genotype_entries(false, false, false, true);
	if (entry_header.has_genotypes == false)
		LOG.error("Require Genotypes in variant file in order to filter genotypes by Filter Flag.");

	filter_genotypes_by_filter_status(filter_flags_to_remove, remove_all);
}

void entry::filter_sites(const set<string> &snps_to_keep, const string &snps_to_keep_file, const string &snps_to_exclude_file, bool keep_then_exclude)
{
	// Filter sites by user provided lists
	if (keep_then_exclude)
	{
		filter_sites_to_keep(snps_to_keep, snps_to_keep_file);
		filter_sites_to_exclude(snps_to_exclude_file);
	}
	else
	{
		filter_sites_to_exclude(snps_to_exclude_file);
		filter_sites_to_keep(snps_to_keep, snps_to_keep_file);
	}
}

void entry::filter_sites_to_keep(const set<string> &snps_to_keep, const string &snps_to_keep_file)
{
	// Filter sites by user provided list
	if(passed_filters == false)
		return;

	if ((snps_to_keep.empty()) && (snps_to_keep_file == ""))
		return;

	if (snps_to_keep_file != "" && local_snps_to_keep.empty())
	{
		ifstream in(snps_to_keep_file.c_str());
		string tmp;
		local_snps_to_keep = snps_to_keep;

		if (!in.is_open())
		{
			LOG.error("Could not open SNPs to Keep file" + snps_to_keep_file, 0);
		}
		while (!in.eof())
		{
			in >> tmp;
			local_snps_to_keep.insert(tmp);
			in.ignore(numeric_limits<streamsize>::max(), '\n');
		}

		in.close();
	}
	parse_basic_entry();

	if ( (local_snps_to_keep.find(ID) == local_snps_to_keep.end()) && (snps_to_keep.find(ID) == snps_to_keep.end()) )
		passed_filters = false;
}

void entry::filter_sites_to_exclude(const string &snps_to_exclude_file)
{
	// Filter sites by user provided list
	if(passed_filters == false)
		return;

	if (snps_to_exclude_file == "")
		return;

	if (snps_to_exclude_file != "" && snps_to_exclude.empty())
	{
		ifstream in(snps_to_exclude_file.c_str());
		string tmp;
		if (!in.is_open())
		{
			LOG.error("Could not open SNPs to Exclude file" + snps_to_exclude_file, 0);
		}
		while (!in.eof())
		{
			in >> tmp;
			snps_to_exclude.insert(tmp);
			in.ignore(numeric_limits<streamsize>::max(), '\n');
		}
		in.close();
	}

	parse_basic_entry();

	if (snps_to_exclude.find(ID) != snps_to_exclude.end())
		passed_filters = false;
}

void entry::filter_sites_by_quality(double min_quality)
{
	// Filter sites by quality
	if (passed_filters == false)
		return;

	if (min_quality < 0)
		return;

	parse_basic_entry(true);

	string alt_allele = get_ALT_allele(0);
	// The QUAL field has different definitions depending on the state of the
	// alternative allele. Here I treat them separately, although in this case
	// it is unnecessary.
	if ((alt_allele == ".") || (alt_allele == ""))
	{	// The case that the alternative allele is unknown
		// QUAL is -10log_10 p(variant)
			if (QUAL < min_quality)
				passed_filters = false;
	}
	else
	{	// The normal case
		// QUAL is -10log_10 p(no variant)
		if (QUAL < min_quality)
			passed_filters = false;
	}
}

void entry::filter_sites_by_mean_depth(double min_mean_depth, double max_mean_depth)
{
	// Filter sites by mean depth
	if (passed_filters == false)
		return;

	if ((min_mean_depth <= 0) && (max_mean_depth == numeric_limits<double>::max()))
		return;

	int depth;

	unsigned int N_indv_included = 0;
	double depth_sum = 0.0;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		if (include_genotype[ui] == true)
		{
			parse_genotype_entry(ui, false, false, true);

			if (entry_header.has_genotypes == false)
				LOG.error("Require Genotypes in VCF file in order to filter sites by mean depth");

			depth = get_indv_DEPTH(ui);
			if (depth >= 0)
			{
				depth_sum += depth;
			}
			N_indv_included++;
		}
	}
	double mean_depth = depth_sum / N_indv_included;

	if ((mean_depth < min_mean_depth) || (mean_depth > max_mean_depth))
		passed_filters = false;
}

void entry::filter_sites_by_position(const string &chr, int start_pos, int end_pos)
{
	// Filter sites by user provided position range
	if (passed_filters == false)
		return;

	if ((chr == "") || ((start_pos == -1) && (end_pos==numeric_limits<int>::max())))
		return;

	parse_basic_entry();

	if (CHROM == chr)
	{
		if ((POS < start_pos) || (POS > end_pos))
			passed_filters = false;
	}
	else
		passed_filters = false;
}

void entry::filter_sites_by_positions(const string &positions_file, const string &exclude_positions_file)
{
	// Filter sites by a user defined file containing a list of positions
	if (passed_filters == false)
		return;

	if ((positions_file == "") && (exclude_positions_file == ""))
		return;

	int idx;
	if (keep_positions.empty() && positions_file != "")
	{
		string chr;
		int pos1;
		unsigned int N_chr=chr_to_idx.size();
		stringstream ss;
		string line;
		unsigned int gzMAX_LINE_LEN = 1024*1024;
		char *gz_readbuffer = new char[gzMAX_LINE_LEN];

		gzFile gz_in = gzopen(positions_file.c_str(), "rb");
		if (gz_in == NULL)
			LOG.error("Could not open Positions file: " + positions_file);

		keep_positions.resize(N_chr);
		while (!gzeof(gz_in))
		{
			line = "";
			bool again = true;
			while (again == true)
			{
				gzgets(gz_in, gz_readbuffer, gzMAX_LINE_LEN);
				line.append(gz_readbuffer);
				if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
					again = false;
			}
			if (line[0] == '#')
				continue;
			line.erase( line.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)

			ss.clear();
			ss.str(line);
			ss >> chr >> pos1;

			if (chr_to_idx.find(chr) == chr_to_idx.end())
			{
				N_chr++;
				chr_to_idx[chr] = (N_chr-1);
				keep_positions.resize(N_chr);
			}

			idx = chr_to_idx[chr];
			keep_positions[idx].insert(pos1);
		}
		gzclose(gz_in);
		delete [] gz_readbuffer;
	}
	if (exclude_positions.empty() && exclude_positions_file != "")
	{
		string chr;
		int pos1;
		unsigned int N_chr=chr_to_idx.size();
		stringstream ss;
		string line;
		unsigned int gzMAX_LINE_LEN = 1024*1024;
		char *gz_readbuffer = new char[gzMAX_LINE_LEN];

		gzFile gz_in = gzopen(exclude_positions_file.c_str(), "rb");
		if (gz_in == NULL)
			LOG.error("Could not open Positions file: " + exclude_positions_file);

		exclude_positions.resize(N_chr);
		while (!gzeof(gz_in))
		{
			line = "";
			bool again = true;
			while (again == true)
			{
				gzgets(gz_in, gz_readbuffer, gzMAX_LINE_LEN);
				line.append(gz_readbuffer);
				if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
					again = false;
			}
			if (line[0] == '#')
				continue;
			line.erase( line.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)

			ss.clear();
			ss.str(line);
			ss >> chr >> pos1;

			if (chr_to_idx.find(chr) == chr_to_idx.end())
			{
				N_chr++;
				chr_to_idx[chr] = (N_chr-1);
				exclude_positions.resize(N_chr);
			}

			idx = chr_to_idx[chr];
			exclude_positions[idx].insert(pos1);
		}
		gzclose(gz_in);
		delete [] gz_readbuffer;
	}

	parse_basic_entry();

	if (!keep_positions.empty())
	{	// Check to see if position is in keep list
			if (chr_to_idx.find(CHROM) == chr_to_idx.end())
				passed_filters = false;
			else
			{
				idx = chr_to_idx[CHROM];
				if (keep_positions[idx].find(POS) == keep_positions[idx].end())
					passed_filters = false;
			}
	}
	if (!exclude_positions.empty())
	{	// Check to see if position is in exclude list
		if (chr_to_idx.find(CHROM) != chr_to_idx.end())
		{
			idx = chr_to_idx[CHROM];
			if (exclude_positions[idx].find(POS) != exclude_positions[idx].end())
                passed_filters = false;
		}
	}
}

void entry::filter_sites_by_overlap_positions(const string &positions_overlap_file, const string &exclude_positions_overlap_file)
{
	// Filter sites by overlapping with a user defined file containing a list of positions
	if (passed_filters == false)
		return;

	if ((positions_overlap_file == "") && (exclude_positions_overlap_file == ""))
		return;

	int idx;
	if (keep_positions.empty() && positions_overlap_file != "")
	{
		string chr;
		int pos1;
		unsigned int N_chr=chr_to_idx.size();
		stringstream ss;
		string line;
		unsigned int gzMAX_LINE_LEN = 1024*1024;
		char *gz_readbuffer = new char[gzMAX_LINE_LEN];

		gzFile gz_in = gzopen(positions_overlap_file.c_str(), "rb");
		if (gz_in == NULL)
			LOG.error("Could not open Positions file: " + positions_overlap_file);

		while (!gzeof(gz_in))
		{
			line = "";
			bool again = true;
			while (again == true)
			{
				gzgets(gz_in, gz_readbuffer, gzMAX_LINE_LEN);
				line.append(gz_readbuffer);
				if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
					again = false;
			}
			if (line[0] == '#')
				continue;
			line.erase( line.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)

			ss.clear();
			ss.str(line);
			ss >> chr >> pos1;

			if (chr_to_idx.find(chr) == chr_to_idx.end())
			{
				N_chr++;
				chr_to_idx[chr] = (N_chr-1);
				keep_positions.resize(N_chr);
			}

			idx = chr_to_idx[chr];
			keep_positions[idx].insert(pos1);
		}
		gzclose(gz_in);
		delete [] gz_readbuffer;
	}
	if (exclude_positions.empty() && exclude_positions_overlap_file != "")
	{
		string chr;
		int pos1;
		unsigned int N_chr=0;
		stringstream ss;
		string line;
		unsigned int gzMAX_LINE_LEN = 1024*1024;
		char *gz_readbuffer = new char[gzMAX_LINE_LEN];

		gzFile gz_in = gzopen(exclude_positions_overlap_file.c_str(), "rb");
		if (gz_in == NULL)
			LOG.error("Could not open Positions file: " + exclude_positions_overlap_file);

		while (!gzeof(gz_in))
		{
			line = "";
			bool again = true;
			while (again == true)
			{
				gzgets(gz_in, gz_readbuffer, gzMAX_LINE_LEN);
				line.append(gz_readbuffer);
				if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
					again = false;
			}
			if (line[0] == '#')
				continue;
			line.erase( line.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)

			ss.clear();
			ss.str(line);
			ss >> chr >> pos1;

			if (chr_to_idx.find(chr) == chr_to_idx.end())
			{
				N_chr++;
				chr_to_idx[chr] = (N_chr-1);
				exclude_positions.resize(N_chr);
			}

			idx = chr_to_idx[chr];
			exclude_positions[idx].insert(pos1);
		}
		gzclose(gz_in);
		delete [] gz_readbuffer;
	}

	parse_basic_entry();

	if (!keep_positions.empty())
	{	// Check to see if position is in keep list
			if (chr_to_idx.find(CHROM) == chr_to_idx.end())
				passed_filters = false;
			else
			{
				idx = chr_to_idx[CHROM];
				bool found=false;

				for (unsigned int ui=POS; ui<(POS+REF.size()); ui++)
					if (keep_positions[idx].find(ui) != keep_positions[idx].end())
					{
						found = true;
						break;
					}

				if (found == false)
					passed_filters = false;
			}
	}
	if (!exclude_positions.empty())
	{	// Check to see if position is in exclude list
		if (chr_to_idx.find(CHROM) != chr_to_idx.end())
		{
			idx = chr_to_idx[CHROM];
			bool found=false;

			for (unsigned int ui=POS; ui<(POS+REF.size()); ui++)
				if (exclude_positions[idx].find(ui) != exclude_positions[idx].end())
					found = true;

				if (found == true)
					passed_filters = false;
		}
	}
}

void entry::filter_sites_by_chromosome(const set<string> &chrs_to_keep, const set<string> &chrs_to_exclude)
{
	if (passed_filters == false)
		return;

	if (chrs_to_keep.empty() && chrs_to_exclude.empty())
		return;

	parse_basic_entry();

	if (!chrs_to_keep.empty())
	{
		if (chrs_to_keep.find(CHROM) == chrs_to_keep.end())
			passed_filters = false;
	}
	else
	{
		if (chrs_to_exclude.find(CHROM) != chrs_to_exclude.end())
			passed_filters = false;
	}
}

void entry::filter_sites_by_BED_file(const string &bed_file, bool BED_exclude)
{
	// Filter sites depending on positions in a BED file.
	if (passed_filters == false)
		return;

	if (bed_file == "")
		return;

	int pos1, pos2, idx;
	if (lims.empty())
	{
		ifstream BED(bed_file.c_str());
		if (!BED.is_open())
			LOG.error("Could not open BED file: " + bed_file);

		string chr;
		unsigned int N_chr=chr_to_idx.size();
		BED.ignore(numeric_limits<streamsize>::max(), '\n');	// Ignore header
		unsigned int N_BED_entries=0;
		while (!BED.eof())
		{
			BED >> chr >> pos1 >> pos2;
			BED.ignore(numeric_limits<streamsize>::max(), '\n');

			if (chr_to_idx.find(chr) == chr_to_idx.end())
			{
				N_chr++;
				chr_to_idx[chr] = (N_chr-1);
				lims.resize(N_chr);
			}

			idx = chr_to_idx[chr];
			lims[idx].push_back(make_pair(pos1,pos2));
			N_BED_entries++;
		}
		BED.close();

		LOG.printLOG("\tRead " + output_log::int2str(N_BED_entries) + " BED file entries.\n");

		for (unsigned int ui=0; ui<lims.size(); ui++)
			sort(lims[ui].begin(), lims[ui].end());
	}
	vector<unsigned int> min_ui(lims.size(), 0);
	parse_basic_entry(true);

	pos1 = POS;
	pos2 = pos1;
	unsigned int N_alleles = get_N_alleles();
	for (int i=0; i<(int)N_alleles; i++)
		pos2 = max(pos2, (int)(pos1 + get_allele(i).length() - 1));

	if (BED_exclude == false)
	{	// Exclude sites not in BED file
		if (chr_to_idx.find(CHROM) == chr_to_idx.end())
			passed_filters = false;
		else
		{
			idx = chr_to_idx[CHROM];
			bool found=false;
			unsigned int max_ui = lims[idx].size();
			for (unsigned int ui=min_ui[idx]; ui<max_ui; ui++)
			{	// No need to start this loop at zero every time...
				if (((pos1 > lims[idx][ui].first) && (pos1 <= lims[idx][ui].second)) ||	// Start pos inside bin
					((pos2 > lims[idx][ui].first) && (pos2 <= lims[idx][ui].second)) ||	// End pos inside bin
					((pos1 <= lims[idx][ui].first) && (pos2 >= lims[idx][ui].second)))	// Variant spans bin
				{
					found=true;
					break;
				}
				else if (pos1 > lims[idx][ui].second)
					min_ui[idx] = ui+1;
			}
			if (found == false)
				passed_filters = false;
		}
	}
	else
	{	// Exclude sites in BED file
		if (chr_to_idx.find(CHROM) != chr_to_idx.end())
		{
			idx = chr_to_idx[CHROM];
			bool found=false;
			unsigned int max_ui = lims[idx].size();
			for (unsigned int ui=min_ui[idx]; ui<max_ui; ui++)
			{	// No need to start this loop at zero every time...
				if (((pos1 > lims[idx][ui].first) && (pos1 <= lims[idx][ui].second)) ||	// Start pos inside bin
					((pos2 > lims[idx][ui].first) && (pos2 <= lims[idx][ui].second)) ||	// End pos inside bin
					((pos1 <= lims[idx][ui].first) && (pos2 >= lims[idx][ui].second)))	// Variant spans bin
				{
					found=true;
					break;
				}
				else if (pos1 > lims[idx][ui].second)
					min_ui[idx] = ui+1;
			}
			if (found == true)
				passed_filters = false;
		}
	}
}

void entry::filter_sites_by_mask(const string &mask_file, bool invert_mask, int min_kept_mask_value)
{
	// Filter sites on the basis of a fasta-like mask file.
	if (passed_filters == false || mask_file == "")
		return;

	if (!mask.is_open())
	{
		mask.open(mask_file.c_str());
		mask_chr = "";
		mask_line = "";
		mask_pos = 1;

		if (!mask.is_open())
			LOG.error("Could not open mask file: " + mask_file);
	}

	string line;
	string next_chr="";
	unsigned int next_pos = 0;

	parse_basic_entry();
	next_chr = CHROM;

	while (mask_chr != next_chr && !mask.eof())
	{
		getline(mask, line);
		line.erase( line.find_last_not_of(" \t") + 1);

		if (line[0] == '>')
		{
			mask_chr = line.substr(1, line.find_first_of(" \t")-1);
			mask_pos = 1;
			getline(mask, line);
			mask_line = line;
		}
	}

	if (next_chr == mask_chr)
		next_pos = (unsigned)POS;
	else
	{
		passed_filters = false;
		return;
	}

	while (next_pos > (mask_pos + mask_line.size()) && !mask.eof())
	{
		getline(mask, line);
		line.erase( line.find_last_not_of(" \t") + 1);

		if (line[0] == '>')
		{
			mask_chr = line.substr(1, line.find_first_of(" \t")-1);
			mask_pos = 1;
			passed_filters = false;
			return;
		}
		else
		{
			mask_pos += mask_line.size();
			mask_line = line;
		}
	}

	if (next_chr == mask_chr && next_pos <= (mask_pos+mask_line.size()))
	{
		char mask_base = mask_line[next_pos-mask_pos]-48;
		bool keep = (mask_base <= min_kept_mask_value);

		if (invert_mask == true)
			keep = !keep;

		if (keep == false)
			passed_filters = false;
	}
	else
		passed_filters = false;
}

void entry::filter_sites_by_number_of_alleles(int min_alleles, int max_alleles)
{
	// Filter sites by the number of alleles (e.g. 2 for bi-allelic)
	if (passed_filters == false)
		return;

	if ((min_alleles <= 0) && (max_alleles == numeric_limits<int>::max()))
		return;

	int N_alleles;
	parse_basic_entry(true);
	N_alleles = get_N_alleles();
	if ((N_alleles < min_alleles) || (N_alleles > max_alleles))
		passed_filters = false;
}

void entry::filter_sites_by_frequency_and_call_rate(double min_maf, double max_maf,
		double min_non_ref_af, double max_non_ref_af,
		double min_non_ref_af_any, double max_non_ref_af_any,
		double min_site_call_rate)
{
	// Filter sites so that all allele frequencies are between limits
	if (passed_filters == false)
		return;

	if ((min_maf <= 0.0) && (max_maf >= 1.0) &&
			(min_site_call_rate <= 0) &&
			(min_non_ref_af <= 0.0) && (max_non_ref_af >= 1.0) &&
			(min_non_ref_af_any <= 0.0) && (max_non_ref_af_any >= 1.0))
		return;

	unsigned int N_alleles;
	unsigned int N_non_missing_chr;
	parse_basic_entry(true);
	parse_genotype_entries(true);

	if (GT_idx == -1)
		LOG.error("Require Genotypes in variant file to filter by frequency and/or call rate");

	N_alleles = get_N_alleles();

	vector<int> allele_counts;
	get_allele_counts(allele_counts, N_non_missing_chr);

	double freq, folded_freq;
	double maf=numeric_limits<double>::max();
	int N_failed = 0;
	for (unsigned int ui=0; ui<N_alleles; ui++)
	{
		freq = allele_counts[ui] / (double)N_non_missing_chr;
		folded_freq = min(freq, 1.0 - freq);

		maf = min(maf, folded_freq);
		if ((ui > 0) && ((freq < min_non_ref_af) || (freq > max_non_ref_af)))
			passed_filters = false;

		if ((ui > 0) && ((freq < min_non_ref_af_any) || (freq > max_non_ref_af_any)))
			N_failed++;
	}

	if (((min_non_ref_af > 0.0) || (max_non_ref_af < 1.0)) && (N_failed == (N_alleles-1)))
		passed_filters = false;

	if ((maf < min_maf) || (maf > max_maf))
		passed_filters = false;

	double call_rate = N_non_missing_chr / double(get_N_chr());

	if (call_rate < min_site_call_rate)
		passed_filters = false;
}

void entry::filter_sites_by_allele_type(bool keep_only_indels, bool remove_indels)
{
	if (passed_filters == false)
		return;

	if ((keep_only_indels == false) && (remove_indels == false))
		return;
	if ((keep_only_indels == true) && (remove_indels == true))
		LOG.error("Can't both keep and remove all indels!");

	string allele;
	unsigned int ref_len, N_alleles;
	bool is_indel;

	parse_basic_entry(true);
	is_indel = false;
	allele = REF;
	ref_len = allele.size();
	if (ref_len != 1)
		is_indel = true;
	N_alleles = get_N_alleles();
	for (unsigned int ui=1; ui<N_alleles; ui++)
	{
		get_allele(ui, allele);
		if (allele.size() != ref_len)
		{
			is_indel = true;
			break;
		}
	}

	if (keep_only_indels == true)
	{
		if (is_indel == false)
			passed_filters = false;
	}
	else if (remove_indels == true)
	{
		if (is_indel == true)
			passed_filters = false;
	}
}

void entry::filter_sites_by_allele_count(double min_mac, double max_mac,
		double min_non_ref_ac, double max_non_ref_ac,
		double min_non_ref_ac_any, double max_non_ref_ac_any,
		double max_missing_call_count)
{
	// Filter sites so that all allele counts are between limits
	if (passed_filters == false)
		return;

	if ((min_mac <= 0) && (max_mac == numeric_limits<int>::max()) &&
			(min_non_ref_ac <= 0) && (max_non_ref_ac == numeric_limits<int>::max()) &&
			(min_non_ref_ac_any <= 0) && (max_non_ref_ac_any == numeric_limits<int>::max()) &&
			(max_missing_call_count == numeric_limits<int>::max()))
		return;

	unsigned int N_alleles, N_chr, N_non_missing_chr;
	parse_basic_entry(true);
	parse_genotype_entries(true);

	if (entry_header.has_genotypes == false)
		LOG.error("Require Genotypes in variant file to filter by allele counts and/or missing data");

	N_alleles = get_N_alleles();

	if (N_alleles <= 1 && min_mac > 0)
		passed_filters = false;

	vector<int> allele_counts;
	get_allele_counts(allele_counts, N_non_missing_chr);
	N_chr = get_N_chr();

	int mac = numeric_limits<int>::max();
	int N_failed = 0;
	for (unsigned int ui=0; ui<N_alleles; ui++)
	{
		mac = min(allele_counts[ui], mac);
		if ((ui > 0) && ((allele_counts[ui] < min_non_ref_ac) || (allele_counts[ui] > max_non_ref_ac)))
			passed_filters = false;

		if ((ui > 0) && ((allele_counts[ui] < min_non_ref_ac_any) || (allele_counts[ui] > max_non_ref_ac_any)))
			N_failed++;
	}

	if (((min_non_ref_ac_any > 0) || (max_non_ref_ac_any < numeric_limits<int>::max())) && (N_failed == (N_alleles-1)))
		passed_filters = false;

	if ((mac < min_mac) || (mac > max_mac))
		passed_filters = false;

	if ((N_chr-N_non_missing_chr) > max_missing_call_count)
		passed_filters = false;
}

void entry::filter_sites_by_HWE_pvalue(double min_HWE_pvalue)
{
	// Filter sites by HWE p-value
	// Note this assumes Biallelic SNPs.
	if(passed_filters == false)
		return;

	if (min_HWE_pvalue <= 0)
		return;

	unsigned int b11, b12, b22;
	double p_hwe, p_lo, p_hi;

	parse_basic_entry(true);
	parse_genotype_entries(true);

	if (entry_header.has_genotypes == false)
		LOG.error("Require Genotypes in variant file to filter sites by HWE.");

	get_genotype_counts(b11, b12, b22);
	entry::SNPHWE(b12, b11, b22, p_hwe, p_lo, p_hi);

	if (p_hwe < min_HWE_pvalue)
		passed_filters = false;
}

void entry::filter_sites_by_filter_status(const set<string> &filter_flags_to_remove, const set<string> &filter_flags_to_keep, bool remove_all)
{
	// Filter sites by entries in the FILTER field.
	if (passed_filters == false)
		return;

	if ((remove_all == false) && (filter_flags_to_remove.empty()) && (filter_flags_to_keep.empty()))
		return;

	vector<string> FILTERs;
	unsigned int N_to_remove = filter_flags_to_remove.size();
	unsigned int N_to_keep = filter_flags_to_keep.size();
	parse_basic_entry(false, true);
	get_FILTER_vector(FILTERs);

	if (N_to_keep > 0)
	{
		bool keep = false;
		for (unsigned int ui=0; ui<FILTERs.size(); ui++)
			if (filter_flags_to_keep.find(FILTERs[ui]) != filter_flags_to_keep.end())
			{
				keep = true; break;
			}

		passed_filters = keep;
	}
	if (passed_filters==false)
		return;

	if ( (FILTERs.size() >= 1) && (FILTERs[0] == "PASS") )
		return;
	else if ((remove_all == true) && (!FILTERs.empty()))
		passed_filters = false;
	else if (N_to_remove > 0)
	{
		for (unsigned int ui=0; ui<FILTERs.size(); ui++)
			if (filter_flags_to_remove.find(FILTERs[ui]) != filter_flags_to_remove.end())
				passed_filters = false;
	}
}

void entry::filter_sites_by_phase()
{
	// Filter out sites with unphased entries
	// TODO: Alter this to allow for a max/min level of unphased-ness.
	if (passed_filters == false)
		return;

	unsigned int count_unphased = 0;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		parse_genotype_entry(ui, true);

		if (get_indv_PHASE(ui) != '|')
			count_unphased++;
	}

	if (count_unphased > 0)
		passed_filters = false;
}

void entry::filter_sites_by_thinning(int min_SNP_distance)
{
	// Filter sites so that no two SNPs are within some minimum distance
	if (passed_filters == false)
		return;

	if (min_SNP_distance < 1)
		return;

	parse_basic_entry();
	if (CHROM == thin_chrom)
	{
		int distance_from_last_SNP = POS - thin_pos;
		if (distance_from_last_SNP < min_SNP_distance)
			passed_filters = false;
	}
	if (passed_filters == true)
		thin_pos = POS;
	thin_chrom = CHROM;
}

void entry::filter_sites_by_INFO(const set<string> &flags_to_remove, const set<string> &flags_to_keep)
{
	// Filter sites by entries in the INFO field.
	if (passed_filters == false)
		return;

	if ((flags_to_remove.empty()) && (flags_to_keep.empty()))
		return;

	string value;
	unsigned int N_to_remove = flags_to_remove.size();
	unsigned int N_to_keep = flags_to_keep.size();

	parse_basic_entry(false, false, true);

	if (N_to_keep > 0)
	{
		bool keep = false;
		for (set<string>::iterator it=flags_to_keep.begin(); it != flags_to_keep.end(); ++it)
		{
			if (entry_header.INFO_map[ entry_header.INFO_reverse_map[*it] ].Type != Flag)
				LOG.error("Using INFO flag filtering on non flag type "+*it+" will not work correctly.");
			else
			{
				value = get_INFO_value(*it);
				if (value == "1")
					keep = true;
			}
		}
		passed_filters = keep;
	}

	if (passed_filters==false)
		return;

	if (N_to_remove > 0)
	{
		for (set<string>::iterator it=flags_to_remove.begin(); it != flags_to_remove.end(); ++it)
		{
			if (entry_header.INFO_map[ entry_header.INFO_reverse_map[*it] ].Type != Flag)
				LOG.error("Using INFO flag filtering on non flag type "+*it+" will not work correctly.");
			else
			{
				value = get_INFO_value(*it);
				if (value == "1")
				{
					passed_filters = false;
					continue;
				}
			}
		}
	}
}
