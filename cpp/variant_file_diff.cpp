/*
 * variant_file_diff.cpp
 *
 *  Created on: Oct 30, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "variant_file.h"

void variant_file::return_indv_union(variant_file &file2, map<string, pair< int, int> > &combined_individuals, const string &indv_ID_map_file)
{
	map<string, string> indv_map;
	bool use_map = false;
	if (indv_ID_map_file != "")
	{
		LOG.printLOG("Reading individual mapping file. ");
		ifstream map(indv_ID_map_file.c_str());
		if (!map.is_open())
			LOG.error("Could not open map file: " + indv_ID_map_file);
		while (!map.eof())
		{
			string indv1, indv2;
			map >> indv1 >> indv2;
			map.ignore(numeric_limits<streamsize>::max(), '\n');
			if ((indv1 != "") && (indv1.substr(0,1) != "#"))
			{
				indv_map[indv1] = indv2;
			}
		}
		map.close();
		use_map = true;
		LOG.printLOG("Read " + LOG.int2str(indv_map.size()) + " entries.\n");
	}

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		if (include_indv[ui] == true)
		{
			combined_individuals[meta_data.indv[ui]] = make_pair<int,int>((int)ui, -1);
		}

	for (unsigned int ui=0; ui<file2.meta_data.N_indv; ui++)
		if (file2.include_indv[ui] == true)
		{
			string indv_id = file2.meta_data.indv[ui];
			if (use_map == true)
			{
				if (indv_map.find(indv_id) != indv_map.end())
					indv_id = indv_map[indv_id];
			}
			if (combined_individuals.find(indv_id) != combined_individuals.end())
			{
				combined_individuals[indv_id].second = ui;
			}
			else
				combined_individuals[indv_id] = make_pair<int,int>(-1, (int)ui);
		}
}

void variant_file::output_sites_in_files(const parameters &params, variant_file &diff_variant_file)
{
	string CHROM;
	vector<char> variant_line;
	entry *e1 = get_entry_object();
	entry *e2 = diff_variant_file.get_entry_object();

	bool new_e1 = true;
	bool new_e2 = true;
	string CHROM1 = "";
	string CHROM2 = "";
	string curr_CHROM = "";
	vector<string> all_CHROM;
	int POS1 = -1;
	int POS2 = -1;
	string REF1 = "";
	string REF2 = "";
	string ALT1 = "";
	string ALT2 = "";
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0, N_overlap_SNPs = 0;
	string output_file = params.output_prefix + ".diff.sites_in_files";

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

	ostream sites_in_files(buf);

	sites_in_files << "CHROM\tPOS1\tPOS2\tIN_FILE\tREF1\tREF2\tALT1\tALT2" << endl;

	LOG.printLOG("Comparing sites in VCF files...\n");
	while(true)
	{
		if(new_e1)
		{
			while(!eof())
			{
				get_entry(variant_line);
				e1->reset(variant_line);
				N_entries += e1->apply_filters(params);

				if(!e1->passed_filters)
					continue;
				N_kept_entries++;
				e1->parse_basic_entry(true);
				CHROM1 = e1->get_CHROM();
				POS1 = e1->get_POS();
				REF1 = e1->get_REF();
				ALT1 = e1->get_ALT();
				break;
			}
			new_e1 = false;
		}
		if(new_e2)
		{
			while(!diff_variant_file.eof())
			{
				diff_variant_file.get_entry(variant_line);
				e2->reset(variant_line);
				diff_variant_file.N_entries += e2->apply_filters(params);

				if(!e2->passed_filters)
					continue;
				diff_variant_file.N_kept_entries++;
				e2->parse_basic_entry(true);
				CHROM2 = e2->get_CHROM();
				POS2 = e2->get_POS();
				REF2 = e2->get_REF();
				ALT2 = e2->get_ALT();
				break;
			}
			new_e2 = false;
		}

		if(eof() && diff_variant_file.eof())
			break;
		else if(diff_variant_file.eof())
		{
			if(CHROM1 == curr_CHROM)
			{
				sites_in_files << CHROM1 << "\t" << POS1 << "\t.\t1\t" << REF1 << "\t.\t" << ALT1 << "\t." << endl;
				N_SNPs_file1_only++;
				new_e1 = true;

			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM1;
					all_CHROM.push_back(CHROM1);
					sites_in_files << CHROM1 << "\t" << POS1 << "\t.\t1\t" << REF1 << "\t.\t" << ALT1 << "\t." << endl;
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
		}
		else if(eof())
		{
			if(CHROM2 == curr_CHROM)
			{
				sites_in_files << CHROM2 << "\t.\t" << POS2 << "\t2\t.\t" << REF2 << "\t.\t" << ALT2 << endl;
				N_SNPs_file2_only++;
				new_e2 = true;

			}
			else
			{
				if (find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM2;
					all_CHROM.push_back(CHROM2);
					sites_in_files << CHROM2 << "\t.\t" << POS2 << "\t2\t.\t" << REF2 << "\t.\t" << ALT2 << endl;
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else if(CHROM1 == CHROM2)
		{
			if (CHROM1 != curr_CHROM)
			{
				curr_CHROM = CHROM1;
				all_CHROM.push_back(curr_CHROM);
			}

			if(POS1 == POS2)
			{
				if ((REF1 == "N") || (REF1 == ".") || (REF1 == "") )
					REF1 = REF2;
				if ((REF2 == "N") || (REF2 == ".") || (REF2 == "") )
					REF2 = REF1;

				new_e1 = true;
				new_e2 = true;
				if ((REF1 != REF2) && (REF2 != "N") && (REF1 != "N") && (REF1 != ".") && (REF2 != ".") && (REF1 != "") && (REF2 != ""))
				{
					sites_in_files << CHROM1 << "\t" << POS1 << "\t" << POS2 << "\tO\t" << REF1 << "\t" << REF2 << "\t" << ALT1 << "\t" << ALT2 << endl;
					N_overlap_SNPs++;
				}
				else
				{
					sites_in_files << CHROM1 << "\t" << POS1 << "\t" << POS2 << "\tB\t" << REF1 << "\t" << REF2 << "\t" << ALT1 << "\t" << ALT2 << endl;
					N_common_SNPs++;
				}
			}
			else if(POS1 < POS2)
			{
				if (POS2 < (POS1+REF1.size()))
				{
					sites_in_files << CHROM1 << "\t" << POS1 << "\t" << POS2 << "\tO\t" << REF1 << "\t" << REF2 <<"\t" << ALT1 << "\t" << ALT2 << endl;
					N_overlap_SNPs++;
					new_e1 = true;
					new_e2 = true;
				}
				else
				{
					sites_in_files << CHROM1 << "\t" << POS1 << "\t.\t1\t" << REF1 << "\t.\t" << ALT1 << "\t." << endl;
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
			else
			{
				if (POS1 < (POS2+REF2.size()))
				{
					sites_in_files << CHROM1 << "\t" << POS1 << "\t" << POS2 << "\tO\t" << REF1 << "\t" << REF2 <<"\t" << ALT1 << "\t" << ALT2 << endl;
					N_overlap_SNPs++;
					new_e1 = true;
					new_e2 = true;
				}
				else
				{
					sites_in_files << CHROM2 << "\t.\t" << POS2 << "\t2\t.\t" << REF2 << "\t.\t" << ALT2 << endl;
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else
		{
			if (CHROM1 == curr_CHROM)
			{
				sites_in_files << CHROM1 << "\t" << POS1 << "\t.\t1\t" << REF1 << "\t.\t" << ALT1 << "\t." << endl;
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else if (CHROM2 == curr_CHROM)
			{
				sites_in_files << CHROM2 << "\t.\t" << POS2 << "\t2\t.\t" << REF2 << "\t.\t" << ALT2 << endl;
				N_SNPs_file2_only++;
				new_e2 = true;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");

				LOG.error("Cannot determine chromosomal ordering of files, both files must contain the same chromosomes to use the diff functions.\nFound "+CHROM1+" in file 1 and "+CHROM2+" in file 2.\nUse option --not-chr to filter out chromosomes only found in one file.");
			}
		}
	}

	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " sites common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " sites only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " sites only in second file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_overlap_SNPs) + " non-matching overlapping sites.\n");
	delete e1;
	delete e2;
}

void variant_file::output_indv_in_files(const parameters &params, variant_file &diff_variant_file)
{
	LOG.printLOG("Comparing individuals in VCF files...\n");

	string output_file = params.output_prefix + ".diff.indv_in_files";

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
	out << "INDV\tFILES" << endl;

	// Build a list of individuals contained in each file
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, params.diff_indv_map_file);

	unsigned int N_combined_indv = combined_individuals.size();
	unsigned int N[3]={0,0,0};
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		if ((combined_individuals_it->second.first != -1) && (combined_individuals_it->second.second != -1))
		{
			N[0]++;
			out << combined_individuals_it->first << "\tB" << endl;
		}
		else if (combined_individuals_it->second.first != -1)
		{
			N[1]++;
			out << combined_individuals_it->first << "\t1" << endl;
		}
		else if (combined_individuals_it->second.second != -1)
		{
			N[2]++;
			out << combined_individuals_it->first << "\t2" << endl;
		}
		else
			LOG.error("Unhandled case");
	}

	LOG.printLOG("N_combined_individuals:\t" + output_log::int2str(N_combined_indv) + "\n");
	LOG.printLOG("N_individuals_common_to_both_files:\t" + output_log::int2str(N[0]) + "\n");
	LOG.printLOG("N_individuals_unique_to_file1:\t" + output_log::int2str(N[1]) + "\n");
	LOG.printLOG("N_individuals_unique_to_file2:\t" + output_log::int2str(N[2]) + "\n");
}

void variant_file::output_discordance_by_indv(const parameters &params, variant_file &diff_variant_file)
{
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, params.diff_indv_map_file);

	LOG.printLOG("Outputting Discordance By Individual...\n");
	map<string, pair<int, int> > indv_sums;

	vector<char> variant_line;
	int indv1, indv2;
	entry * e1 = get_entry_object();
	entry * e2 = diff_variant_file.get_entry_object();
	string CHROM;
	bool new_e1 = true;
	bool new_e2 = true;
	string CHROM1 = "";
	string CHROM2 = "";
	string curr_CHROM = "";
	vector<string> all_CHROM;
	int POS1 = -1;
	int POS2 = -1;
	string REF1 = "";
	string REF2 = "";
	string ALT1 = "";
	string ALT2 = "";
	bool alleles_match = false;
	pair<string, string> genotype1, genotype2;
	pair<int,int> geno_ids1, geno_ids2;
	pair<string, string> missing_genotype(".",".");
	pair<int, int> missing_id(-1,-1);
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;

	while(true)
	{
		if(new_e1)
		{
			while(!eof())
			{
				get_entry(variant_line);
				e1->reset(variant_line);
				N_entries += e1->apply_filters(params);

				if(!e1->passed_filters)
					continue;
				N_kept_entries++;
				e1->parse_basic_entry(true);
				CHROM1 = e1->get_CHROM();
				POS1 = e1->get_POS();
				REF1 = e1->get_REF();
				ALT1 = e1->get_ALT();
				break;
			}
			new_e1 = false;
		}
		if(new_e2)
		{
			while(!diff_variant_file.eof())
			{
				diff_variant_file.get_entry(variant_line);
				e2->reset(variant_line);
				diff_variant_file.N_entries += e2->apply_filters(params);

				if(!e2->passed_filters)
					continue;
				diff_variant_file.N_kept_entries++;
				e2->parse_basic_entry(true);
				CHROM2 = e2->get_CHROM();
				POS2 = e2->get_POS();
				REF2 = e2->get_REF();
				ALT2 = e2->get_ALT();
				break;
			}
			new_e2 = false;
		}

		if(eof() && diff_variant_file.eof())
			break;
		else if(diff_variant_file.eof())
		{
			if(CHROM1 == curr_CHROM)
			{
				N_SNPs_file1_only++;
				new_e1 = true;

			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM1;
					all_CHROM.push_back(CHROM1);
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
		}
		else if(eof())
		{
			if(CHROM2 == curr_CHROM)
			{
				N_SNPs_file2_only++;
				new_e2 = true;

			}
			else
			{
				if (find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM2;
					all_CHROM.push_back(CHROM2);
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else if(CHROM1 == CHROM2)
		{
			if (CHROM1 != curr_CHROM)
			{
				curr_CHROM = CHROM1;
				all_CHROM.push_back(curr_CHROM);
			}

			if(POS1 == POS2)
			{
				new_e1 = true;
				new_e2 = true;
				N_common_SNPs++;
			}
			else if(POS1 < POS2)
			{
				new_e1 = true;
				N_SNPs_file1_only++;
			}
			else
			{
				new_e2 = true;
				N_SNPs_file2_only++;
			}
		}
		else
		{
			if (CHROM1 == curr_CHROM)
			{
				new_e1 = true;
				N_SNPs_file1_only++;
			}
			else if (CHROM2 == curr_CHROM)
			{
				new_e2 = true;
				N_SNPs_file2_only++;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");

				LOG.error("Cannot determine chromosomal ordering of files, both files must contain the same chromosomes to use the diff functions.\nFound "+CHROM1+" in file 1 and "+CHROM2+" in file 2.\nUse option --not-chr to filter out chromosomes only found in one file.");
			}
		}

		if(new_e1 && new_e2)
		{
			if (REF1 == "N")
				REF1 = REF2;
			if (REF2 == "N")
				REF2 = REF1;

			if ((REF1.size() != REF2.size()) || ((REF1 != REF2) && (REF2 != "N") && (REF1 != "N")))
			{
				LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
				continue;
			}
			alleles_match = (ALT1 == ALT2) && (REF1 == REF2);
			e1->parse_full_entry(true);
			e1->parse_genotype_entries(true);
			e2->parse_full_entry(true);
			e2->parse_genotype_entries(true);

			for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
			{
				indv1 = combined_individuals_it->second.first;
				indv2 = combined_individuals_it->second.second;

				if ((indv1 == -1) || (indv2 == -1))
					continue;	// Individual not found in one of the files

				if (alleles_match)
				{	// Alleles match, so can compare ids instead of strings
					e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
					e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

					if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
					{
						indv_sums[combined_individuals_it->first].first++;
						if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
								((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
						{	// Match
							// Don't do anything
						}
						else
						{	// Mismatch
							indv_sums[combined_individuals_it->first].second++;
						}
					}
					else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
					{	// Both missing
						// Don't do anything.
					}
					else if (geno_ids1 != missing_id)
					{	// Genotype 1 is not missing, genotype 2 is.
						// Don't do anything.
					}
					else if (geno_ids2 != missing_id)
					{	// Genotype 2 is not missing, genotype 1 is.
						// Don't do anything.
					}
					else
						LOG.error("Unknown condition");
				}
				else
				{	// Alleles don't match, so need to be more careful and compare strings
					e1->get_indv_GENOTYPE_strings(indv1, genotype1);
					e2->get_indv_GENOTYPE_strings(indv2, genotype2);

					if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
					{	// No missing data
						indv_sums[combined_individuals_it->first].first++;
						if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
								((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
						{	// Match
							// Don't do anything
						}
						else
						{	// Mismatch
							indv_sums[combined_individuals_it->first].second++;
						}
					}
					else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
					{	// Both missing
						// Don't do anything
					}
					else if (genotype1 != missing_genotype)
					{	// Genotype 1 is not missing, genotype 2 is.
						// Don't do anything
					}
					else if (genotype2 != missing_genotype)
					{	// Genotype 2 is not missing, genotype 1 is.
						// Don't do anything
					}
					else
						LOG.error("Unknown condition");
				}
			}
		}
	}

	string output_file = params.output_prefix + ".diff.indv";
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
	out << "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	int N, N_discord;
	double discordance;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		out << combined_individuals_it->first;
		N = indv_sums[combined_individuals_it->first].first;
		N_discord = indv_sums[combined_individuals_it->first].second;
		discordance = N_discord / double(N);
		out << "\t" << N << "\t" << N_discord << "\t" << discordance << endl;
	}
	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " sites common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " sites only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " sites only in second file.\n");

	delete e1;
	delete e2;
}

void variant_file::output_discordance_by_site(const parameters &params, variant_file &diff_variant_file)
{
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, params.diff_indv_map_file);

	LOG.printLOG("Outputting Discordance By Site...\n");

	string CHROM;
	vector<char> variant_line;
	int indv1, indv2;
	entry *e1 = get_entry_object();
	entry *e2 = diff_variant_file.get_entry_object();

	bool new_e1 = true;
	bool new_e2 = true;
	string CHROM1 = "";
	string CHROM2 = "";
	string curr_CHROM = "";
	vector<string> all_CHROM;
	int POS1 = -1;
	int POS2 = -1;
	string REF1 = "";
	string REF2 = "";
	string ALT1 = "";
	string ALT2 = "";
	bool alleles_match = false;
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;

	string output_file = params.output_prefix + ".diff.sites";
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

	ostream diffsites(buf);
	diffsites << "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	while(true)
	{
		if(new_e1)
		{
			while(!eof())
			{
				get_entry(variant_line);
				e1->reset(variant_line);
				N_entries += e1->apply_filters(params);

				if(!e1->passed_filters)
					continue;
				N_kept_entries++;
				e1->parse_basic_entry(true);
				CHROM1 = e1->get_CHROM();
				POS1 = e1->get_POS();
				REF1 = e1->get_REF();
				ALT1 = e1->get_ALT();
				break;
			}
			new_e1 = false;
		}
		if(new_e2)
		{
			while(!diff_variant_file.eof())
			{
				diff_variant_file.get_entry(variant_line);
				e2->reset(variant_line);
				diff_variant_file.N_entries += e2->apply_filters(params);

				if(!e2->passed_filters)
					continue;
				diff_variant_file.N_kept_entries++;
				e2->parse_basic_entry(true);
				CHROM2 = e2->get_CHROM();
				POS2 = e2->get_POS();
				REF2 = e2->get_REF();
				ALT2 = e2->get_ALT();
				break;
			}
			new_e2 = false;
		}

		if(eof() && diff_variant_file.eof())
			break;
		else if(diff_variant_file.eof())
		{
			if(CHROM1 == curr_CHROM)
			{
				diffsites << CHROM1 << "\t" << POS1 << "\t1\t";
				N_SNPs_file1_only++;
				new_e1 = true;

			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM1;
					all_CHROM.push_back(CHROM1);
					diffsites << CHROM1 << "\t" << POS1 << "\t1\t";
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
		}
		else if(eof())
		{
			if(CHROM2 == curr_CHROM)
			{
				diffsites << CHROM2 << "\t" << POS2 << "\t2\t";
				N_SNPs_file2_only++;
				new_e2 = true;

			}
			else
			{
				if (find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM2;
					all_CHROM.push_back(CHROM2);
					diffsites << CHROM2 << "\t" << POS2 << "\t2\t";
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else if(CHROM1 == CHROM2)
		{
			if (CHROM1 != curr_CHROM)
			{
				curr_CHROM = CHROM1;
				all_CHROM.push_back(curr_CHROM);
			}

			if(POS1 == POS2)
			{
				if ((REF1 == "N") || (REF1 == ".") || (REF1 == "") )
					REF1 = REF2;
				if ((REF2 == "N") || (REF2 == ".") || (REF2 == "") )
					REF2 = REF1;

				new_e1 = true;
				new_e2 = true;
				if ((REF1 != REF2) && (REF2 != "N") && (REF1 != "N") && (REF1 != ".") && (REF2 != ".") && (REF1 != "") && (REF2 != ""))
				{
					LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
					continue;
				}
				diffsites << CHROM1 << "\t" << POS1 << "\tB\t";
				N_common_SNPs++;
			}
			else if(POS1 < POS2)
			{
				diffsites << CHROM1 << "\t" << POS1 << "\t1\t";
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else
			{
				diffsites << CHROM2 << "\t" << POS2 << "\t2\t";
				N_SNPs_file2_only++;
				new_e2 = true;
			}
		}
		else
		{
			if (CHROM1 == curr_CHROM)
			{
				diffsites << CHROM1 << "\t" << POS1 << "\t1\t";
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else if (CHROM2 == curr_CHROM)
			{
				diffsites << CHROM2 << "\t" << POS2 << "\t2\t";
				N_SNPs_file2_only++;
				new_e2 = true;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");

				LOG.error("Cannot determine chromosomal ordering of files, both files must contain the same chromosomes to use the diff functions.\nFound "+CHROM1+" in file 1 and "+CHROM2+" in file 2.\nUse option --not-chr to filter out chromosomes only found in one file.");
			}
		}

		pair<string, string> genotype1, genotype2;
		pair<int,int> geno_ids1, geno_ids2;
		pair<string, string> missing_genotype(".",".");
		pair<int, int> missing_id(-1,-1);

		unsigned int N_common_called=0;	// Number of genotypes called in both files
		unsigned int N_missing_1=0, N_missing_2=0;
		unsigned int N_discord=0;
		unsigned int N_concord_non_missing=0;

		if(new_e1 && new_e2)
		{
			alleles_match = (ALT1 == ALT2) && (REF1 == REF2);
			diffsites << alleles_match;

			e1->parse_full_entry(true);
			e1->parse_genotype_entries(true);
			e2->parse_full_entry(true);
			e2->parse_genotype_entries(true);

			for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
			{
				indv1 = combined_individuals_it->second.first;
				indv2 = combined_individuals_it->second.second;

				if ((indv1 == -1) || (indv2 == -1))
					continue;	// Individual not found in one of the files

				if (alleles_match)
				{	// Alleles match, so can compare ids instead of strings
					e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
					e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

					if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
					{
						N_common_called++;
						if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
								((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
						{	// Match
							N_concord_non_missing++;
						}
						else
						{	// Mismatch
							N_discord++;
						}
					}
					else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
					{	// Both missing
						N_missing_1++; N_missing_2++;
					}
					else if (geno_ids1 != missing_id)
					{	// Genotype 1 is not missing, genotype 2 is.
						N_missing_2++;
					}
					else if (geno_ids2 != missing_id)
					{	// Genotype 2 is not missing, genotype 1 is.
						N_missing_1++;
					}
					else
						LOG.error("Unknown condition");
				}
				else
				{	// Alleles don't match, so need to be more careful and compare strings
					e1->get_indv_GENOTYPE_strings(indv1, genotype1);
					e2->get_indv_GENOTYPE_strings(indv2, genotype2);

					if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
					{	// No missing data
						N_common_called++;
						if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
								((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
						{	// Match
							N_concord_non_missing++;
						}
						else
						{	// Mismatch
							N_discord++;
						}
					}
					else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
					{	// Both missing
						N_missing_1++; N_missing_2++;
					}
					else if (genotype1 != missing_genotype)
					{	// Genotype 1 is not missing, genotype 2 is.
						N_missing_2++;
					}
					else if (genotype2 != missing_genotype)
					{	// Genotype 2 is not missing, genotype 1 is.
						N_missing_1++;
					}
					else
						LOG.error("Unknown condition");
				}
			}
		}
		else
			diffsites << "0";

		double discordance = N_discord / double(N_common_called);
		diffsites << "\t" << N_common_called << "\t" << N_discord << "\t" << discordance;
		diffsites << endl;
	}

	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " sites common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " sites only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " sites only in second file.\n");
	delete e1;
	delete e2;
}

void variant_file::output_discordance_matrix(const parameters &params, variant_file &diff_variant_file)
{
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, params.diff_indv_map_file);

	LOG.printLOG("Outputting Discordance Matrix\n\tFor bi-allelic loci, called in both files, with matching alleles only...\n");
	string CHROM;
	vector<char> variant_line;
	int indv1, indv2;
	entry *e1 = get_entry_object();
	entry *e2 = diff_variant_file.get_entry_object();

	bool new_e1 = true;
	bool new_e2 = true;
	string CHROM1 = "";
	string CHROM2 = "";
	string curr_CHROM = "";
	vector<string> all_CHROM;
	int POS1 = -1;
	int POS2 = -1;
	string REF1 = "";
	string REF2 = "";
	string ALT1 = "";
	string ALT2 = "";
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;
	vector<vector<int> > discordance_matrix(4, vector<int>(4, 0));

	if (combined_individuals.size() <= 0)
		LOG.error("No overlapping individuals can be found.");

	while(true)
	{
		if(new_e1)
		{
			while(!eof())
			{
				get_entry(variant_line);
				e1->reset(variant_line);
				N_entries += e1->apply_filters(params);

				if(!e1->passed_filters)
					continue;
				N_kept_entries++;
				e1->parse_basic_entry(true);
				CHROM1 = e1->get_CHROM();
				POS1 = e1->get_POS();
				REF1 = e1->get_REF();
				ALT1 = e1->get_ALT();
				break;
			}
			new_e1 = false;
		}
		if(new_e2)
		{
			while(!diff_variant_file.eof())
			{
				diff_variant_file.get_entry(variant_line);
				e2->reset(variant_line);
				diff_variant_file.N_entries += e2->apply_filters(params);

				if(!e2->passed_filters)
					continue;
				diff_variant_file.N_kept_entries++;
				e2->parse_basic_entry(true);
				CHROM2 = e2->get_CHROM();
				POS2 = e2->get_POS();
				REF2 = e2->get_REF();
				ALT2 = e2->get_ALT();
				break;
			}
			new_e2 = false;
		}

		if(eof() && diff_variant_file.eof())
			break;
		else if(diff_variant_file.eof())
		{
			if(CHROM1 == curr_CHROM)
			{
				N_SNPs_file1_only++;
				new_e1 = true;

			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM1;
					all_CHROM.push_back(CHROM1);
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
		}
		else if(eof())
		{
			if(CHROM2 == curr_CHROM)
			{
				N_SNPs_file2_only++;
				new_e2 = true;

			}
			else
			{
				if (find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM2;
					all_CHROM.push_back(CHROM2);
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else if(CHROM1 == CHROM2)
		{
			if (CHROM1 != curr_CHROM)
			{
				curr_CHROM = CHROM1;
				all_CHROM.push_back(curr_CHROM);
			}

			if(POS1 == POS2)
			{
				if ((REF1 == "N") || (REF1 == ".") || (REF1 == "") )
					REF1 = REF2;
				if ((REF2 == "N") || (REF2 == ".") || (REF2 == "") )
					REF2 = REF1;

				new_e1 = true;
				new_e2 = true;
				if ((REF1 != REF2) && (REF2 != "N") && (REF1 != "N") && (REF1 != ".") && (REF2 != ".") && (REF1 != "") && (REF2 != ""))
				{
					LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
					continue;
				}
				N_common_SNPs++;
			}
			else if(POS1 < POS2)
			{
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else
			{
				N_SNPs_file2_only++;
				new_e2 = true;
			}
		}
		else
		{
			if (CHROM1 == curr_CHROM)
			{
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else if (CHROM2 == curr_CHROM)
			{
				N_SNPs_file2_only++;
				new_e2 = true;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");

				LOG.error("Cannot determine chromosomal ordering of files, both files must contain the same chromosomes to use the diff functions.\nFound "+CHROM1+" in file 1 and "+CHROM2+" in file 2.\nUse option --not-chr to filter out chromosomes only found in one file.");
			}
		}

		if(new_e1 && new_e2)
		{
			if (e1->get_N_alleles() != 2 || e2->get_N_alleles() != 2)
				continue;

			if (ALT1 != ALT2)
			{
				LOG.one_off_warning("Non-matching ALT. Skipping all such sites.");
				continue;
			}
			e1->parse_full_entry(true);
			e1->parse_genotype_entries(true);

			e2->parse_full_entry(true);
			e2->parse_genotype_entries(true);

			pair<int,int> geno_ids1, geno_ids2;
			int N1, N2;

			for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
			{
				indv1 = combined_individuals_it->second.first;
				indv2 = combined_individuals_it->second.second;

				if ((indv1 == -1) || (indv2 == -1))
				{
					LOG.one_off_warning("Non-matching individual found. Skipping all such combinations.");
					continue;	// Individual not found in one of the files
				}

				// Alleles match, so can compare ids instead of strings
				e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
				e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

				if (((geno_ids1.first != -1) && (geno_ids1.second == -1)) ||
						((geno_ids2.first != -1) && (geno_ids2.second == -1)))
				{	// Haploid
					LOG.one_off_warning("***Warning: Haploid chromosomes not counted!***");
					continue;
				}

				N1 = geno_ids1.first + geno_ids1.second;
				N2 = geno_ids2.first + geno_ids2.second;

				if ((N1 == -1) || (N1 < -2) || (N1 > 2))
					LOG.error("Unhandled case");
				if ((N2 == -1) || (N2 < -2) || (N2 > 2))
					LOG.error("Unhandled case");

				if (N1 == -2)
					N1 = 3;

				if (N2 == -2)
					N2 = 3;

				discordance_matrix[N1][N2]++;
			}
		}
	}

	string output_file = params.output_prefix + ".diff.discordance_matrix";
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

	out << "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1" << endl;
	out << "N_0/0_file2\t" << discordance_matrix[0][0] << "\t" << discordance_matrix[1][0] << "\t" << discordance_matrix[2][0] << "\t" << discordance_matrix[3][0] << endl;
	out << "N_0/1_file2\t" << discordance_matrix[0][1] << "\t" << discordance_matrix[1][1] << "\t" << discordance_matrix[2][1] << "\t" << discordance_matrix[3][1] << endl;
	out << "N_1/1_file2\t" << discordance_matrix[0][2] << "\t" << discordance_matrix[1][2] << "\t" << discordance_matrix[2][2] << "\t" << discordance_matrix[3][2] << endl;
	out << "N_./._file2\t" << discordance_matrix[0][3] << "\t" << discordance_matrix[1][3] << "\t" << discordance_matrix[2][3] << "\t" << discordance_matrix[3][3] << endl;

	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " sites common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " sites only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " sites only in second file.\n");
	delete e1;
	delete e2;
}

void variant_file::output_switch_error(const parameters &params, variant_file &diff_variant_file)
{
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, params.diff_indv_map_file);

	LOG.printLOG("Outputting Phase Switch Errors...\n");

	vector<char> variant_line;
	int indv1, indv2;
	entry *e1 = get_entry_object();
	entry *e2 = diff_variant_file.get_entry_object();

	bool new_e1 = true;
	bool new_e2 = true;
	string CHROM1 = "";
	string CHROM2 = "";
	string curr_CHROM = "";
	vector<string> all_CHROM;
	int POS1 = -1;
	int POS2 = -1;
	string REF1 = "";
	string REF2 = "";
	string ALT1 = "";
	string ALT2 = "";
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;

	string output_file = params.output_prefix + ".diff.switch";
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

	ostream switcherror(buf);

	switcherror << "CHROM\tPOS_START\tPOS_END\tINDV" << endl;

	unsigned int N_combined_indv = combined_individuals.size();
	vector<int> N_phased_het_sites(N_combined_indv, 0);
	vector<int> N_switch_errors(N_combined_indv, 0);

	pair<string, string> missing_genotype(".",".");
	pair<string, int> missing_loc(".",-1);
	vector<pair<string, string> > prev_geno_file1(N_combined_indv, missing_genotype);
	vector<pair<string, string> > prev_geno_file2(N_combined_indv, missing_genotype);
	vector<pair<string, int> > prev_pos_file1(N_combined_indv, missing_loc);
	vector<pair<string, int> > prev_pos_file2(N_combined_indv, missing_loc);
	pair<string, string> file1_hap1, file1_hap2, file2_hap1;

	if (N_combined_indv <= 0)
		LOG.error("No overlapping individuals can be found.");

	while(true)
	{
		if(new_e1)
		{
			while(!eof())
			{
				get_entry(variant_line);
				e1->reset(variant_line);
				N_entries += e1->apply_filters(params);

				if(!e1->passed_filters)
					continue;
				N_kept_entries++;
				e1->parse_basic_entry(true);
				CHROM1 = e1->get_CHROM();
				POS1 = e1->get_POS();
				REF1 = e1->get_REF();
				ALT1 = e1->get_ALT();
				break;
			}
			new_e1 = false;
		}
		if(new_e2)
		{
			while(!diff_variant_file.eof())
			{
				diff_variant_file.get_entry(variant_line);
				e2->reset(variant_line);
				diff_variant_file.N_entries += e2->apply_filters(params);

				if(!e2->passed_filters)
					continue;
				diff_variant_file.N_kept_entries++;
				e2->parse_basic_entry(true);
				CHROM2 = e2->get_CHROM();
				POS2 = e2->get_POS();
				REF2 = e2->get_REF();
				ALT2 = e2->get_ALT();
				break;
			}
			new_e2 = false;
		}

		if(eof() && diff_variant_file.eof())
			break;
		else if(diff_variant_file.eof())
		{
			if(CHROM1 == curr_CHROM)
			{
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM1;
					all_CHROM.push_back(CHROM1);
					N_SNPs_file1_only++;
					new_e1 = true;
				}
			}
		}
		else if(eof())
		{
			if(CHROM2 == curr_CHROM)
			{
				N_SNPs_file2_only++;
				new_e2 = true;
			}
			else
			{
				if (find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");
				else
				{
					curr_CHROM = CHROM2;
					all_CHROM.push_back(CHROM2);
					N_SNPs_file2_only++;
					new_e2 = true;
				}
			}
		}
		else if(CHROM1 == CHROM2)
		{
			if (CHROM1 != curr_CHROM)
			{
				curr_CHROM = CHROM1;
				all_CHROM.push_back(curr_CHROM);
			}

			if(POS1 == POS2)
			{
				N_common_SNPs++;
				new_e1 = true;
				new_e2 = true;
			}
			else if(POS1 < POS2)
			{
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else
			{
				N_SNPs_file2_only++;
				new_e2 = true;
			}
		}
		else
		{
			if (CHROM1 == curr_CHROM)
			{
				N_SNPs_file1_only++;
				new_e1 = true;
			}
			else if (CHROM2 == curr_CHROM)
			{
				N_SNPs_file2_only++;
				new_e2 = true;
			}
			else
			{
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM1) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM1+" in file 1 appears to be out of order.");
				if(find(all_CHROM.begin(), all_CHROM.end(), CHROM2) != all_CHROM.end())
					LOG.error("Both files must be sorted in the same chromosomal order.\n"+CHROM2+" in file 2 appears to be out of order.");

				LOG.error("Cannot determine chromosomal ordering of files, both files must contain the same chromosomes to use the diff functions.\nFound "+CHROM1+" in file 1 and "+CHROM2+" in file 2.\nUse option --not-chr to filter out chromosomes only found in one file.");			}
		}
		if(new_e1 && new_e2)
		{
			e1->parse_full_entry(true);
			e1->parse_genotype_entries(true);

			e2->parse_full_entry(true);
			e2->parse_genotype_entries(true);

			pair<string, string> genotype1, genotype2;
			pair<string, string> missing_genotype(".",".");

			unsigned int N_common_called=0;	// Number of genotypes called in both files
			unsigned int indv_count=0;

			for (combined_individuals_it=combined_individuals.begin();
					combined_individuals_it!=combined_individuals.end();
					++combined_individuals_it, indv_count++)
			{
				indv1 = combined_individuals_it->second.first;
				indv2 = combined_individuals_it->second.second;

				if ((indv1 == -1) || (indv2 == -1))
				{
					LOG.one_off_warning("Non-matching individual found. Skipping all such combinations.");
					continue;	// Individual not found in one of the files
				}

				e1->get_indv_GENOTYPE_strings(indv1, genotype1);
				e2->get_indv_GENOTYPE_strings(indv2, genotype2);

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{	// No missing data
					N_common_called++;
					if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
							((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
					{	// Have a matching genotypes in files 1 and 2
						if (genotype1.first != genotype1.second)
						{	// It's a heterozgote
							char phase1, phase2;
							phase1 = e1->get_indv_PHASE(indv1);
							phase2 = e2->get_indv_PHASE(indv2);
							if ((phase1 == '|') && (phase2 == '|'))
							{	// Calculate Phasing error (switch error)
								N_phased_het_sites[indv_count]++;
								file1_hap1 = make_pair<string,string>((string)prev_geno_file1[indv_count].first, (string)genotype1.first);
								file1_hap2 = make_pair<string,string>((string)prev_geno_file1[indv_count].second, (string)genotype1.second);
								file2_hap1 = make_pair<string,string>((string)prev_geno_file2[indv_count].first, (string)genotype2.first);

								if ((file2_hap1 != file1_hap1) && (file2_hap1 != file1_hap2))
								{	// Must be a switch error
									string indv_id;
									N_switch_errors[indv_count]++;
									if (indv1 != -1)
										indv_id = meta_data.indv[indv1];
									else
										indv_id = diff_variant_file.meta_data.indv[indv2];

									if (prev_pos_file1[indv_count].first == prev_pos_file2[indv_count].first)
									{
										if (prev_pos_file1[indv_count].second <= prev_pos_file2[indv_count].second)
											switcherror << prev_pos_file1[indv_count].first << "\t" << prev_pos_file1[indv_count].second << "\t" << POS1 << "\t" << indv_id << endl;
										else
											switcherror << prev_pos_file1[indv_count].first << "\t" << prev_pos_file2[indv_count].second << "\t" << POS1 << "\t" << indv_id << endl;
									}
								}
								prev_geno_file1[indv_count] = genotype1;
								prev_geno_file2[indv_count] = genotype2;
								prev_pos_file1[indv_count] = std::pair<string,int>(CHROM1,POS1);
								prev_pos_file2[indv_count] = std::pair<string,int>(CHROM2,POS2);
							}
						}
					}
				}
			}
		}
	}
	delete e1;
	delete e2;

	output_file = params.output_prefix + ".diff.indv.switch";
	ofstream idiscord(output_file.c_str());
	if (!idiscord.is_open())
		LOG.error("Could not open Individual Discordance File: " + output_file, 3);

	idiscord << "INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH" << endl;
	unsigned int indv_count=0;
	double switch_error;
	string indv_id;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		indv1 = combined_individuals_it->second.first;
		indv2 = combined_individuals_it->second.second;

		if (indv1 != -1)
			indv_id = meta_data.indv[indv1];
		else
			indv_id = diff_variant_file.meta_data.indv[indv2];

		if (N_phased_het_sites[indv_count] > 0)
			switch_error = double(N_switch_errors[indv_count]) / N_phased_het_sites[indv_count];
		else
			switch_error = 0;
		idiscord << indv_id << "\t" << N_phased_het_sites[indv_count] << "\t" << N_switch_errors[indv_count] << "\t" << switch_error << endl;

		indv_count++;
	}
	idiscord.close();

	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " sites common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " sites only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " sites only in second file.\n");
}
