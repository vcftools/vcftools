/*
 * variant_file_filters.cpp
 *
 *      Author: amarcketta
 */

#include "variant_file.h"

void variant_file::apply_filters(const parameters &params)
{
	filter_individuals(params.indv_to_keep, params.indv_to_exclude, params.indv_keep_files, params.indv_exclude_files);
	filter_individuals_randomly(params.max_N_indv);
}

void variant_file::filter_individuals(const set<string> &indv_to_keep, const set<string> &indv_to_exclude, const vector<string> &indv_to_keep_filenames, const vector<string> &indv_to_exclude_filenames, bool keep_then_exclude)
{
	// Filter individuals by user provided lists
	if (keep_then_exclude)
	{
		filter_individuals_by_keep_list(indv_to_keep, indv_to_keep_filenames);
		filter_individuals_by_exclude_list(indv_to_exclude, indv_to_exclude_filenames);
	}
	else
	{
		filter_individuals_by_exclude_list(indv_to_exclude, indv_to_exclude_filenames);
		filter_individuals_by_keep_list(indv_to_keep, indv_to_keep_filenames);
	}
}

void variant_file::filter_individuals_by_keep_list(const set<string> &indv_to_keep, const vector<string> &indv_to_keep_filenames)
{
	// Filter individuals by user provided list
	if ((indv_to_keep_filenames.size() == 0) && (indv_to_keep.size() == 0))
		return;

	LOG.printLOG("Keeping individuals in 'keep' list\n");
	set<string> indv_to_keep_copy = indv_to_keep;
	if (indv_to_keep_filenames.size() != 0)
	{
		for (unsigned int ui=0; ui<indv_to_keep_filenames.size(); ui++)
		{
			ifstream infile(indv_to_keep_filenames[ui].c_str());
			if (!infile.is_open())
				LOG.error("Could not open Individual file:" + indv_to_keep_filenames[ui], 1);
			string line;
			string tmp_indv;
			stringstream ss;
			while (!infile.eof())
			{
				getline(infile, line);
				ss.str(line);
				ss >> tmp_indv;
				indv_to_keep_copy.insert(tmp_indv);
				ss.clear();
			}
			infile.close();
		}
	}
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (indv_to_keep_copy.find(meta_data.indv[ui]) == indv_to_keep_copy.end())
			include_indv[ui] = false;
	}
}

void variant_file::filter_individuals_by_exclude_list(const set<string> &indv_to_exclude, const vector<string> &indv_to_exclude_filenames)
{
	// Filter individuals by user provided list
	if ((indv_to_exclude_filenames.size() == 0) && (indv_to_exclude.size() == 0))
		return;
	LOG.printLOG("Excluding individuals in 'exclude' list\n");
	set<string> indv_to_exclude_copy = indv_to_exclude;
	if (indv_to_exclude_filenames.size() != 0)
	{
		for (unsigned int ui=0; ui<indv_to_exclude_filenames.size(); ui++)
		{
			ifstream infile(indv_to_exclude_filenames[ui].c_str());
			if (!infile.is_open())
				LOG.error("Could not open Individual file:" + indv_to_exclude_filenames[ui], 1);
			string line;
			string tmp_indv;
			stringstream ss;
			while (!infile.eof())
			{
				getline(infile, line);
				ss.str(line);
				ss >> tmp_indv;
				indv_to_exclude_copy.insert(tmp_indv);
				ss.clear();
			}
			infile.close();
		}
	}
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (indv_to_exclude_copy.find(meta_data.indv[ui]) != indv_to_exclude_copy.end())
			include_indv[ui] = false;
	}
}

void variant_file::filter_individuals_randomly(int max_N_indv)
{
	// Filter individuals randomly until have a random subset
	if (max_N_indv < 0)
		return;
	LOG.printLOG("Filtering Individuals Randomly\n");

	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in variant file filter individuals.");

	unsigned int N_kept_indv = N_kept_individuals();

	srand ( time(NULL) );
	vector<unsigned int> keep_index(N_kept_indv);
	int count = 0;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == true)
		{
			keep_index[count] = ui;
			count++;
		}
	}

	random_shuffle(keep_index.begin(), keep_index.end());			// Get a random order
	keep_index.resize(min(max_N_indv, (signed)keep_index.size()));	// Only keep a subset

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		bool found = false;
		for (unsigned int uj=0; uj<keep_index.size(); uj++)
		{
			if (keep_index[uj] == ui)
			{
				found = true;
			}
		}
		if (found == false)
			include_indv[ui] = false;
	}
}
