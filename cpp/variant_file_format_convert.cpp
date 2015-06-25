/*
 * variant_file_format_convert.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */
#include "variant_file.h"

void variant_file::output_as_plink(const parameters &params)
{
	// Output as PLINK formatted PED/MAP files.
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output as PLINK.");

	LOG.printLOG("Writing PLINK PED and MAP files ... \n");
	string ped_file = params.output_prefix + ".ped";
	string map_file = params.output_prefix + ".map";
	int fd = -1;

	vector<ofstream *> tmp_files(meta_data.N_indv);
	vector<string> tmp_filenames(meta_data.N_indv);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
		char tmpname[new_tmp.size()];
		strcpy(tmpname, new_tmp.c_str());

		fd = mkstemp(tmpname);
		if (fd == -1)
			LOG.error(" Could not open temporary file.\n", 12);
		::close(fd);

		ofstream *tmp_file = new ofstream(tmpname);
		if (!tmp_file->good())
			LOG.error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n"
					"Alternatively, try the --plink-tped command.", 12);
		(*tmp_file) << meta_data.indv[ui] << "\t" << meta_data.indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = tmpname;
	}

	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) LOG.error("Could not open output file: " + map_file, 12);
	int POS;
	string ID, CHROM, CHROM2;
	map<string, string> CHROM_to_PLINK;

	if (params.chrom_map_file != "")
	{
		ifstream chrom_map(params.chrom_map_file.c_str());
		if (!chrom_map.is_open())
			LOG.error("Could not open chromosome mapping file: " + params.chrom_map_file);
		string chr, plink;
		unsigned int N_chrom_entries=0;
		while (!chrom_map.eof())
		{
			chrom_map >> chr >> plink;
			CHROM_to_PLINK[chr] = plink;
			N_chrom_entries++;
		}
		chrom_map.close();
		LOG.printLOG("\tRead " + output_log::int2str(N_chrom_entries) + " chromosome mapping file entries.\n");
	}
	else
	{
		for (int i=1; i<23; i++)
		{
			ostringstream convert;
			convert << i;
			CHROM_to_PLINK["chr" + convert.str()] = convert.str();
			CHROM_to_PLINK[convert.str()] = convert.str();
		}
		CHROM_to_PLINK["chrX"] = "X";
		CHROM_to_PLINK["chrY"] = "Y";
		CHROM_to_PLINK["chrXY"] = "XY";
		CHROM_to_PLINK["chrMT"] = "MT";
		CHROM_to_PLINK["chrM"] = "M";
		CHROM_to_PLINK["X"] = "X";
		CHROM_to_PLINK["Y"] = "Y";
		CHROM_to_PLINK["XY"] = "XY";
		CHROM_to_PLINK["MT"] = "MT";
		CHROM_to_PLINK["M"] = "M";
	}

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
	vector<char> variant_line;
	entry *e = get_entry_object();
	ofstream *tmp_file;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() > 2)
		{
			LOG.one_off_warning("\tPLINK: Only outputting biallelic loci.");
			continue;
		}

		POS = e->get_POS();
		ID = e->get_ID();
		CHROM = e->get_CHROM();

		if (CHROM_to_PLINK.find(CHROM) == CHROM_to_PLINK.end())
		{
			string tmp = "";
			if (CHROM.compare(0,3,"chr") == 0)
				tmp = CHROM.substr(3, string::npos);
			else
				tmp = CHROM;

			bool isNumber = true;
			for(unsigned int ui=0; ui<tmp.size(); ui++)
			    isNumber = isNumber && isdigit(tmp[ui]);

			if (isNumber)
				CHROM_to_PLINK[CHROM] = tmp;
			else
			{
				LOG.one_off_warning("\nUnrecognized values used for CHROM: " + CHROM + " - Replacing with 0.");
				CHROM_to_PLINK[CHROM] = "0";
			}
		}

		CHROM2 = CHROM_to_PLINK[CHROM];
		if (ID == ".")
			MAP << CHROM2 << "\t" << CHROM << ":" << POS << "\t0\t" << POS << endl;
		else
			MAP << CHROM2 << "\t" << ID << "\t0\t" << POS << endl;

		e->get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			genotype = make_pair(-1,-1);
			phase = '/';
			if (e->include_genotype[ui] == true)
			{
				e->parse_genotype_entry(ui, true);
				e->get_indv_GENOTYPE_ids(ui, genotype);
				phase = e->get_indv_PHASE(ui);
			}

			if (genotype.first < 0)
				(*tmp_file) << "\t0";
			else if (genotype.first > 1)
				LOG.error("File contains entries with nonexistent genotypes at " + CHROM + ":" + output_log::int2str(e->get_POS()) );
			else
				(*tmp_file) << "\t" << alleles[genotype.first];

			if (genotype.second < 0)
			{
				if (phase == '/')
					(*tmp_file) << "\t0";
				else if (genotype.first > -1)
					(*tmp_file) << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
				else
					(*tmp_file) << "\t0";
			}
			else if (genotype.second > 1)
				LOG.error("File contains entries with nonexistent genotypes at " + CHROM + ":" + output_log::int2str(e->get_POS()) );
			else
				(*tmp_file) << "\t" << alleles[genotype.second];
		}
	}
	MAP.close();

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) LOG.error("Could not open output file: " + ped_file, 12);
	string tmp_line;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();
		delete tmp_file;

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			LOG.error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n"
					"Alternatively, try the --plink-tped command.", 12);
		getline(read_file, tmp_line);
		PED << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}
	PED.close();

	delete e;
	LOG.printLOG("Done.\n");
}

// Output as Plink Transposed file
void variant_file::output_as_plink_tped(const parameters &params)
{
	// Output as PLINK formatted PED/MAP files.
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output as PLINK TPED.");

	LOG.printLOG("Writing PLINK TPED file ... ");
	string tped_file = params.output_prefix + ".tped";
	string tfam_file = params.output_prefix + ".tfam";

	ofstream TPED(tped_file.c_str());
	if (!TPED.is_open()) LOG.error("Could not open output file: " + tped_file, 12);

	string CHROM, CHROM2;
	map<string, string> CHROM_to_PLINK;
	if (params.chrom_map_file != "")
	{
		ifstream chrom_map(params.chrom_map_file.c_str());
		if (!chrom_map.is_open())
			LOG.error("Could not open chromosome mapping file: " + params.chrom_map_file);
		string chr, plink;
		unsigned int N_chrom_entries=0;
		while (!chrom_map.eof())
		{
			chrom_map >> chr >> plink;
			CHROM_to_PLINK[chr] = plink;
			N_chrom_entries++;
		}
		chrom_map.close();
		LOG.printLOG("\n\tRead " + output_log::int2str(N_chrom_entries) + " chromosome mapping file entries.\n");
	}
	else
	{
		for (int i=1; i<23; i++)
		{
			ostringstream convert;
			convert << i;
			CHROM_to_PLINK["chr" + convert.str()] = convert.str();
			CHROM_to_PLINK[convert.str()] = convert.str();
		}
		CHROM_to_PLINK["chrX"] = "X";
		CHROM_to_PLINK["chrY"] = "Y";
		CHROM_to_PLINK["chrXY"] = "XY";
		CHROM_to_PLINK["chrMT"] = "MT";
		CHROM_to_PLINK["chrM"] = "M";
		CHROM_to_PLINK["X"] = "X";
		CHROM_to_PLINK["Y"] = "Y";
		CHROM_to_PLINK["XY"] = "XY";
		CHROM_to_PLINK["MT"] = "MT";
		CHROM_to_PLINK["M"] = "M";
	}

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
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

		if (e->get_N_alleles() > 2)	// Only output sites with at most one alternative allele
		{
			LOG.one_off_warning("\tPLINK-TPED: Only outputting biallelic loci.");
			continue;
		}

		CHROM = e->get_CHROM();
		if (CHROM_to_PLINK.find(CHROM) == CHROM_to_PLINK.end())
		{
			string tmp = "";
			if (CHROM.compare(0,3,"chr") == 0)
				tmp = CHROM.substr(3, string::npos);
			else
				tmp = CHROM;

			bool isNumber = true;
			for(unsigned int ui=0; ui<tmp.size(); ui++)
			    isNumber = isNumber && isdigit(tmp[ui]);

			if (isNumber)
				CHROM_to_PLINK[CHROM] = tmp;
			else
			{
				LOG.one_off_warning("\nUnrecognized values used for CHROM: " + CHROM + " - Replacing with 0.");
				CHROM_to_PLINK[CHROM] = "0";
			}
		}
		CHROM2 = CHROM_to_PLINK[CHROM];

		if (e->get_ID() == ".")
			TPED << CHROM2 << "\t" << e->get_CHROM() << ":" << e->get_POS() << "\t0\t" << e->get_POS();
		else
			TPED << CHROM2 << "\t" << e->get_ID() << "\t0\t" << e->get_POS();

		e->get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			genotype = make_pair(-1,-1);
			phase = '/';
			if (e->include_genotype[ui] == true)
			{
				e->parse_genotype_entry(ui, true);
				e->get_indv_GENOTYPE_ids(ui, genotype);
				phase = e->get_indv_PHASE(ui);
			}

			if (genotype.first < 0)
				TPED << "\t0";
			else if (genotype.first > 1)
				LOG.error("File contains entries with nonexistent genotypes at " + CHROM + ":" + output_log::int2str(e->get_POS()) );
			else
				TPED << "\t" << alleles[genotype.first];

			if (genotype.second < 0)
			{
				if (phase == '/')
					TPED << "\t0";
				else if (genotype.first > -1)
					TPED << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
				else
					TPED << "\t0";
			}
			else if (genotype.second > 1)
				LOG.error("File contains entries with nonexistent genotypes at " + CHROM + ":" + output_log::int2str(e->get_POS()) );
			else
				TPED << "\t" << alleles[genotype.second];
		}
		TPED << endl;
	}
	TPED.close();

	LOG.printLOG("Writing PLINK TFAM file ... ");
	ofstream TFAM(tfam_file.c_str());
	if (!TFAM.is_open()) LOG.error("Could not open output file: " + tfam_file, 12);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		TFAM << meta_data.indv[ui] << "\t" << meta_data.indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	}

	delete e;
	TFAM.close();
	LOG.printLOG("Done.\n");
}

void variant_file::output_as_012_matrix(const parameters &params)
{
	// Output as PLINK formatted PED/MAP files.
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output as 0/1/2 matrix.");

	LOG.printLOG("Writing 012 matrix files ... ");
	string ped_file = params.output_prefix + ".012";
	string map_file = params.output_prefix + ".012.pos";
	string fam_file = params.output_prefix + ".012.indv";
	int fd = -1;

	ofstream FAM(fam_file.c_str());
	if (!FAM.is_open()) LOG.error("Could not open output file: " + fam_file, 12);
	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) LOG.error("Could not open output file: " + map_file, 12);

	vector<ofstream *> tmp_files(meta_data.N_indv);
	vector<string> tmp_filenames(meta_data.N_indv);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		FAM << meta_data.indv[ui] << endl;

		string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
		char tmpname[new_tmp.size()];
		strcpy(tmpname, new_tmp.c_str());

		fd = mkstemp(tmpname);
		if (fd == -1)
			LOG.error(" Could not open temporary file.\n", 12);
		::close(fd);

		ofstream *tmp_file = new ofstream(tmpname);
		if (!tmp_file->good())
			LOG.error("\n\nCould not open temporary file.\n\n"
			"Most likely this is because the system is not allowing me to open enough temporary files.\n"
			"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		(*tmp_file) << ui;
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = tmpname;
	}
	FAM.close();

	vector<string> alleles;
	pair<int, int> genotype;
	vector<char> variant_line;
	entry *e = get_entry_object();
	ofstream *tmp_file;
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		if (e->get_N_alleles() > 2)
		{
			LOG.one_off_warning("\t012: Only outputting biallelic loci.");
			continue;
		}
		MAP << e->get_CHROM() << "\t" << e->get_POS() << endl;

		e->get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			genotype = make_pair(-1,-1);
			if (e->include_genotype[ui] == true)
			{
				e->parse_genotype_entry(ui, true);
				e->get_indv_GENOTYPE_ids(ui, genotype);
			}

			if ((genotype.first < 0) && (genotype.second < 0))
				(*tmp_file) << "\t-1";	// Missing data
			else if ((genotype.first == 0) && (genotype.second == 0))
				(*tmp_file) << "\t0";	// No copies of the alternative allele
			else
			{
				if ((genotype.first == 1) && (genotype.second == 1))
					(*tmp_file) << "\t2";	// Two copies of the alternative allele
				else
					(*tmp_file) << "\t1";	// Must be one copy of the alternative allele.
			}
		}
	}

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) LOG.error("Could not open output file: " + ped_file, 12);
	string tmp_line;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();
		delete tmp_file;

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			LOG.error("\n\nCould not open temporary file.\n\n"
			"Most likely this is because the system is not allowing me to open enough temporary files.\n"
			"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		getline(read_file, tmp_line);
		PED << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}
	delete e;
	MAP.close();
	PED.close();
	LOG.printLOG("Done.\n");
}

// Output as IMPUTE format
void variant_file::output_as_IMPUTE(const parameters &params)
{
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output IMPUTE format.");

	LOG.printLOG("Outputting in IMPUTE format (bi-allelic, completely phased SNPs only)\n");
	unsigned int ui;
	string legend_file = params.output_prefix + ".impute.legend";
	string haplotype_file = params.output_prefix + ".impute.hap";
	string indv_file = params.output_prefix + ".impute.hap.indv";
	ofstream legend(legend_file.c_str());
	if (!legend.is_open())
		LOG.error("Could not open IMPUTE Legend Output File: " + legend_file, 2);
	legend << "ID pos allele0 allele1" << endl;

	ofstream hap(haplotype_file.c_str());
	if (!hap.is_open())
		LOG.error("Could not open IMPUTE Haplotype Output File: " + haplotype_file, 2);

	ofstream indv_out(indv_file.c_str());
	if (!indv_out.is_open())
		LOG.error("Could not open IMPUTE Individual Output File: " + indv_file, 2);

	for (ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		indv_out << meta_data.indv[ui] << endl;
	}
	indv_out.close();

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
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tIMPUTE: Only outputting biallelic loci.");
			continue;
		}

		// Exclude entries with missing data and/or unphased
		bool missing = false;
		for (ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == false)
			{
				missing = true;
				break;
			}

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);
			if ((alleles.first < 0) || (alleles.second < 0))
			{
				missing = true;
				break;
			}

			if (e->get_indv_PHASE(ui) != '|')
			{
				missing = true;
				break;
			}
		}
		if (missing == true)
			continue;

		if (e->get_ID() == ".")
		{
			legend << e->get_CHROM() << "-" << e->get_POS() << " " << e->get_POS() << " " << e->get_REF() << " " << e->get_ALT_allele(0) << endl;
		}
		else
			legend << e->get_ID() << " " << e->get_POS() << " " << e->get_REF() << " " << e->get_ALT_allele(0) << endl;

		bool first = true;
		for (ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);
			if (first == true)
			{
				hap << alleles.first << " " << alleles.second;
				first = false;
			}
			else
				hap << " " << alleles.first << " " << alleles.second;
		}
		hap << endl;
	}
	delete e;
	hap.close();
	legend.close();
}

void variant_file::output_as_LDhat_phased(const parameters &params)
{
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output LDhat format.");

	LOG.printLOG("Outputting in phased LDhat format\n");

	unsigned int n_sites = 0;
	int max_pos = -1;
	char* fd;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	fd = mktemp(tmpname);
	string locs_tmp_filename(tmpname);
	ofstream *locs_tmp_file = new ofstream(tmpname);
	if (!locs_tmp_file->good())
		LOG.error("\nCould not open temporary file.\n", 12);

	string sites_file = params.output_prefix + ".ldhat.sites";
	string locs_file = params.output_prefix + ".ldhat.locs";

	ofstream sites(sites_file.c_str());
	if (!sites.is_open())
		LOG.error("Could not open LDhat sites Output File: " + sites_file, 2);
	ofstream locs(locs_file.c_str());
	if (!locs.is_open())
		LOG.error("Could not open LDhat locs Output File: " + locs_file, 2);

	unsigned int n_indv = N_kept_individuals();
	pair<int, int> alleles;

	vector<ofstream *> tmp_files(2*meta_data.N_indv);
	vector<string> tmp_filenames(2*meta_data.N_indv);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		char tmpname[new_tmp.size()];
		strcpy(tmpname, new_tmp.c_str());
		fd = mktemp(tmpname);

		ofstream *tmp_file = new ofstream(tmpname);
		if (!tmp_file->good())
		{	// Clean up temp files.
			tmp_file->close(); remove(tmpname);
			locs_tmp_file->close(); remove(locs_tmp_filename.c_str());
			for (unsigned int uj=0; uj<ui; uj++)
			{
				(tmp_files[2*uj])->close(); remove(tmp_filenames[2*uj].c_str());
				(tmp_files[2*uj+1])->close(); remove(tmp_filenames[2*uj+1].c_str());
			}
			LOG.error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		}
		tmp_files[2*ui] = tmp_file;
		tmp_filenames[2*ui] = tmpname;

		char tmpname2[new_tmp.size()];
		strcpy(tmpname2, new_tmp.c_str());
		fd = mktemp(tmpname2);

		ofstream *tmp_file2 = new ofstream(tmpname2);
		if (!tmp_file2->good())
		{	// Clean up temp files.
			tmp_file2->close(); remove(tmpname2);
			locs_tmp_file->close(); remove(locs_tmp_filename.c_str());
			for (unsigned int uj=0; uj<ui; uj++)
			{
				(tmp_files[2*uj])->close(); remove(tmp_filenames[2*uj].c_str());
				(tmp_files[2*uj+1])->close(); remove(tmp_filenames[2*uj+1].c_str());
			}
			(tmp_files[2*ui])->close(); remove(tmp_filenames[2*ui].c_str());
			LOG.error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		}
		tmp_files[2*ui+1] = tmp_file2;
		tmp_filenames[2*ui+1] = tmpname2;
	}

	vector<char> variant_line;
	entry *e = get_entry_object();
	ofstream *tmp_file;
	int POS;

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
			LOG.one_off_warning("\tLDhat: Only outputting biallelic loci.");
			continue;
		}
		POS = e->get_POS();
		max_pos = max(POS, max_pos);
		(*locs_tmp_file) << POS << endl;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);

			for (unsigned int k=0; k<2; k++)
			{
				tmp_file = tmp_files[(2*ui)+k];

				int geno;
				if (k == 0)
					geno = alleles.first;
				else
					geno = alleles.second;

				if ((geno >= 0) && (e->include_genotype[ui]==true))
					(*tmp_file) << geno;
				else
					(*tmp_file) << "?";
			}
		}
		n_sites++;
	}

	locs << n_sites;
	locs.setf(ios::fixed,ios::floatfield);
	locs.precision(4);
	locs << "\t" << max_pos / 1000.0 << "\tL" << endl;
	ifstream locs_read_file(locs_tmp_filename.c_str());
	string tmp_line;
	for (unsigned int ui=0; ui<n_sites; ui++)
	{
		getline(locs_read_file,tmp_line);
		POS = header::str2int(tmp_line);
		locs << POS / 1000.0 << endl;
	}
	locs.close();

	sites << n_indv*2 << "\t" << n_sites << "\t1" << endl;	// Note - this is incorrect for the X-chr.
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		for (unsigned int k=0; k<2; k++)
		{
			ofstream *tmp_file = tmp_files[2*ui+k];
			(*tmp_file) << endl;
			tmp_file->close();
			delete tmp_file;

			ifstream read_file(tmp_filenames[2*ui+k].c_str());
			if (!read_file.good())
				LOG.error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
			getline(read_file, tmp_line);
			sites << ">" << meta_data.indv[ui] << "-" << k << endl;
			sites << tmp_line << endl;
			read_file.close();
			remove(tmp_filenames[2*ui+k].c_str());
		}
	}
	delete e;
	delete locs_tmp_file;
	remove(locs_tmp_filename.c_str());
	sites.close();
}

void variant_file::output_as_LDhat_unphased(const parameters &params)
{
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output LDhat format.");

	LOG.printLOG("Outputting in unphased LDhat format\n");

	unsigned int n_sites = 0;
	int max_pos = -1;
	char* fd;

	string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
	char tmpname[new_tmp.size()];
	strcpy(tmpname, new_tmp.c_str());
	fd = mktemp(tmpname);
	string locs_tmp_filename(tmpname);
	ofstream *locs_tmp_file = new ofstream(tmpname);
	if (!locs_tmp_file->good())
		LOG.error("\nCould not open temporary file.\n", 12);

	string sites_file = params.output_prefix + ".ldhat.sites";
	string locs_file = params.output_prefix + ".ldhat.locs";

	ofstream sites(sites_file.c_str());
	if (!sites.is_open())
		LOG.error("Could not open LDhat sites Output File: " + sites_file, 2);
	ofstream locs(locs_file.c_str());
	if (!locs.is_open())
		LOG.error("Could not open LDhat locs Output File: " + locs_file, 2);

	unsigned int n_indv = N_kept_individuals();
	pair<int, int> alleles;

	vector<ofstream *> tmp_files(meta_data.N_indv);
	vector<string> tmp_filenames(meta_data.N_indv);
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		char tmpname[new_tmp.size()];
		strcpy(tmpname, new_tmp.c_str());
		fd = mktemp(tmpname);
		string filename(tmpname);
		ofstream *tmp_file = new ofstream(tmpname);
		if (!tmp_file->good())
		{	// Clean up temp files
			tmp_file->close(); remove(tmpname);
			locs_tmp_file->close(); remove(locs_tmp_filename.c_str());
			for (unsigned int uj=0; uj<ui; uj++)
			{
				(tmp_files[uj])->close(); remove(tmp_filenames[uj].c_str());
			}
			LOG.error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		}
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = filename;
	}

	vector<char> variant_line;
	entry *e = get_entry_object();
	ofstream *tmp_file;
	int POS;

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
			LOG.one_off_warning("\tLDhat: Only outputting biallelic loci.");
			continue;
		}
		POS = e->get_POS();
		max_pos = max(POS, max_pos);
		(*locs_tmp_file) << POS << endl;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			if (e->include_genotype[ui] == false)
				(*tmp_file) << "?";
			else
			{
				e->parse_genotype_entry(ui, true);
				e->get_indv_GENOTYPE_ids(ui, alleles);

				switch (alleles.first)
				{
				case -2:
					(*tmp_file) << "?"; break;
				case -1:
					(*tmp_file) << "?"; break;
				case 0:
					if (alleles.second == 0)
						(*tmp_file) << 0;
					else if (alleles.second == 1)
						(*tmp_file) << 2;
					else if ((alleles.second == -1) && (e->get_indv_PHASE(ui) == '|'))
						(*tmp_file) << 0;	// Haploid case
					else if (alleles.second == -2)
						(*tmp_file) << 0;	// Haploid case
					else
						(*tmp_file) << '?';
					break;
				case 1:
					if (alleles.second == 0)
						(*tmp_file) << 2;
					else if (alleles.second == 1)
						(*tmp_file) << 1;
					else if ((alleles.second == -1) && (e->get_indv_PHASE(ui) == '|'))
						(*tmp_file) << 1;	// Haploid case
					else if (alleles.second == -2)
						(*tmp_file) << 1;	// Haploid case
					else
						(*tmp_file) << '?';
					break;
				default:
					(*tmp_file) << '?';
					break;
				}
			}
		}
		n_sites++;
	}
	locs << n_sites;
	locs.setf(ios::fixed,ios::floatfield);
	locs.precision(4);
	locs << "\t" << max_pos / 1000.0 << "\tL" << endl;
	ifstream locs_read_file(locs_tmp_filename.c_str());
	string tmp_line;
	for (unsigned int ui=0; ui<n_sites; ui++)
	{
		getline(locs_read_file,tmp_line);
		POS = header::str2int(tmp_line);
		locs << POS / 1000.0 << endl;
	}
	locs.close();

	sites << n_indv << "\t" << n_sites << "\t2" << endl;
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();
		delete tmp_file;

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			LOG.error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		getline(read_file, tmp_line);
		sites << ">" << meta_data.indv[ui] << endl;
		sites << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}
	delete e;
	delete locs_tmp_file;
	remove(locs_tmp_filename.c_str());
	sites.close();
}

// Output INFO fields in tab-delimited format
void variant_file::output_INFO_for_each_site(const parameters &params)
{
	LOG.printLOG("Outputting INFO for each site\n");
	string output_file = params.output_prefix + ".INFO";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open INFO output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS\tREF\tALT";
	for (unsigned int ui=0; ui<params.INFO_to_extract.size(); ui++)
		out << "\t" << params.INFO_to_extract[ui];
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

		e->parse_basic_entry(true, false, true);
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e->get_REF() << "\t" << e->get_ALT();

		for (unsigned int ui=0; ui<params.INFO_to_extract.size(); ui++)
			out << "\t" << e->get_INFO_value(params.INFO_to_extract[ui]);
		out << endl;
	}
	delete e;
}

// Output FORMAT information in tab-delimited format.
void variant_file::output_FORMAT_information(const parameters &params)
{
	string FORMAT_id = params.FORMAT_id_to_extract;
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output FORMAT information.");

	LOG.printLOG("Outputting FORMAT information for " + FORMAT_id + "\n");
	string output_file = params.output_prefix + "." + FORMAT_id + ".FORMAT";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open FORMAT Output file: " + output_file, 7);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "CHROM\tPOS";
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == true)
			out << "\t" << meta_data.indv[ui];
	}
	out << endl;

	string FORMAT_out;
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
		e->parse_full_entry(true);

		if (e->FORMAT_id_exists(FORMAT_id) == false)
			continue;

		out << e->get_CHROM() << "\t" << e->get_POS();

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->read_indv_generic_entry(ui, FORMAT_id, FORMAT_out);
			out << "\t" << FORMAT_out;
		}
		out << endl;
	}
	delete e;
}

// Output genotype likelihoods from GL or PL FORMAT tag, ready for input into BEAGLE
// using the Genotype likelihoods file format.
void variant_file::output_BEAGLE_genotype_likelihoods(const parameters &params, int GL_or_PL)
{
	if (meta_data.has_genotypes == false)
		LOG.error("Require Genotypes in VCF file in order to output BEAGLE genotype likelihoods.");

	if (GL_or_PL == 0)
		LOG.printLOG("Outputting GLs in BEAGLE Genotype Likelihood format (bi-allelic SNPs with GL tags only)\n");
	else if (GL_or_PL == 1)
		LOG.printLOG("Outputting PLs in BEAGLE Genotype Likelihood format (bi-allelic SNPs with PL tags only)\n");
	else
		LOG.error("Unknown GL or PL option.");

	string output_file = params.output_prefix + ".BEAGLE.GL";
	if (GL_or_PL == 1)
		output_file = params.output_prefix + ".BEAGLE.PL";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open Beagle GL/PL Output file: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	out << "marker\talleleA\talleleB";
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
	{
		if (include_indv[ui] == true)
			out << "\t" << meta_data.indv[ui] << "\t" << meta_data.indv[ui] << "\t" << meta_data.indv[ui];
	}
	out << endl;

	string GL_entry, tmp_string;
	vector<char> variant_line;
	entry *e = get_entry_object();
	double lk1, lk2, lk3;
	bool found_GL=false;
	istringstream ss;

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
			LOG.one_off_warning("\tBEAGLE: Only outputting biallelic loci.");
			continue;
		}

		e->parse_full_entry(true);

		if (GL_or_PL == 0)
			if (e->FORMAT_id_exists("GL") == false)
				continue;
		if (GL_or_PL == 1)
			if (e->FORMAT_id_exists("PL") == false)
				continue;

		found_GL = true;

		out << e->get_CHROM() << ":" << e->get_POS() << "\t" << e->get_REF() << "\t" << e->get_ALT();

		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				if (GL_or_PL == 0)
					e->read_indv_generic_entry(ui, "GL", GL_entry);
				else
					e->read_indv_generic_entry(ui, "PL", GL_entry);
				ss.clear();
				ss.str(GL_entry);
				getline(ss, tmp_string, ',');
				lk1 = atof(tmp_string.c_str());
				getline(ss, tmp_string, ',');
				lk2 = atof(tmp_string.c_str());
				getline(ss, tmp_string);
				lk3 = atof(tmp_string.c_str());
				if (GL_or_PL == 0)
					out << "\t" << pow(10,lk1) << "\t" << pow(10,lk2) << "\t" << pow(10,lk3);
				else
					out << "\t" << pow(10,-lk1*0.1) << "\t" << pow(10,-lk2*0.1) << "\t" << pow(10,-lk3*0.1);
			}
			else
			{
				out << "\t1\t1\t1";	// Mark as unknown
			}
		}
		out << endl;
	}
	delete e;
	if (found_GL == false)
		LOG.error("Require GL or PL FORMAT tags in VCF file to output BEAGLE input.");
}
