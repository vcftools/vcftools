/*
 * vcf_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "vcf_file.h"

vcf_file::vcf_file(const parameters &p, bool diff)
{
	if (!diff)
	{
		filename = p.vcf_filename;
		compressed = p.vcf_compressed;
		stream = p.stream_in;
	}
	else
	{
		filename = p.diff_file;
		compressed = p.diff_file_compressed;
		stream = false;
	}

	gzMAX_LINE_LEN = 0;
	N_entries = 0; N_kept_entries = 0;
	meta_data = header();

	if (stream && compressed)
		open_gz();
	else if (stream)
	{
		char first = cin.peek();
		if (first == 0x1f)
			LOG.error("File starts with gzip magic string. Shouldn't you be using --gzvcf?\n");

		file_in = &std::cin;
	}
	else
		open();

	read_header();
	include_indv = vector<bool>(meta_data.N_indv,true);
}

vcf_file::~vcf_file()
{
	close();
}

void vcf_file::read_header()
{
	string line;
	unsigned int line_index = 0;
	line_index += meta_data.add_FILTER_descriptor("ID=PASS,Description=PASS", line_index);

	while (!eof())
	{
		read_line(line);
		if (line[0] == '#')
			if (line[1] == '#')
				meta_data.parse_meta(line, line_index);
			else
			{
				meta_data.parse_header(line);
				return;
			}
		else
			return;
	}
}

void vcf_file::print(const parameters &params)
{
	LOG.printLOG("Outputting VCF file...\n");

	string output_file = params.output_prefix + ".recode.vcf";
	streambuf * buf;
	ofstream temp_out;
	if (!params.stream_out)
	{
		temp_out.open(output_file.c_str(), ios::out);
		if (!temp_out.is_open()) LOG.error("Could not open VCF Output file: " + output_file, 3);
		buf = temp_out.rdbuf();
	}
	else
		buf = cout.rdbuf();

	ostream out(buf);
	for (unsigned int ui=0; ui<meta_data.lines.size(); ui++)
		out << meta_data.lines[ui] << endl;

	out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (meta_data.N_indv > 0)
		out << "\tFORMAT";
	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		if (include_indv[ui])
			out << "\t" << meta_data.indv[ui];
	out << endl;

	vector<char> variant_line;
	entry * e = new vcf_entry(meta_data, include_indv);
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);
		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true, true, true);
		e->parse_full_entry(true);
		e->parse_genotype_entries(true,true,true,true);
		e->print(out, params.recode_INFO_to_keep, params.recode_all_INFO);
	}
	delete e;
}

void vcf_file::print_bcf(const parameters &params)
{
	LOG.printLOG("Outputting BCF file...\n");
	BGZF * out;
	if(!params.stream_out)
	{
		string output_file = params.output_prefix + ".recode.bcf";
		out = bgzf_open(output_file.c_str(), "w");
	}
	else
		out = bgzf_dopen(1, "w");

	string header_str;
	uint32_t len_text = 0;
	vector<char> header;

	char magic[5] = {'B','C','F','\2', '\2'};
	bgzf_write(out, magic, 5);

	if (meta_data.has_idx)
	{
		LOG.warning("VCF file contains IDX values in header. These are being removed for conversion to BCF.");
		meta_data.reprint();
		meta_data.reparse();
	}
	for (unsigned int ui=0; ui<meta_data.lines.size(); ui++)
	{
		for (unsigned int uj=0; uj<meta_data.lines[ui].length(); uj++)
			header.push_back( meta_data.lines[ui][uj] );
		header.push_back('\n');
	}

	if (meta_data.has_contigs == false)
	{
		vector<string> contig_vector;
		get_contigs(params.contigs_file, contig_vector);

		for(unsigned int ui=0; ui<contig_vector.size(); ui++)
		{
			meta_data.add_CONTIG_descriptor(contig_vector[ui].substr(10, contig_vector[ui].size()-8),int(ui));
			for(unsigned int uj=0; uj<contig_vector[ui].size(); uj++)
				header.push_back(contig_vector[ui][uj]);
			header.push_back('\n');
		}
	}

	header_str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (meta_data.N_indv > 0)
		header_str += "\tFORMAT";

	for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		if (include_indv[ui])
		{
			header_str += "\t";
			header_str += meta_data.indv[ui];
		}
	header_str += "\n";

	for (unsigned int ui=0; ui<header_str.length(); ui++)
		header.push_back( header_str[ui] );

	header.push_back( '\0' );
	len_text = header.size();

	bgzf_write(out, (char *)&len_text, sizeof(len_text) );
	bgzf_write(out, (char *)&header[0], len_text );

	vector<char> variant_line;
	entry * e = new vcf_entry(meta_data, include_indv);
	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);
		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true, true, true);
		e->parse_full_entry(true);
		e->parse_genotype_entries(true,true,true,true);
		e->print_bcf(out, params.recode_INFO_to_keep, params.recode_all_INFO);
	}
	delete e;
	bgzf_close(out);
}

void vcf_file::open()
{
	struct stat buf;

	int i = stat(filename.c_str(), &buf);
	if (i != 0)
	{
		perror("stat error");
		LOG.error("Can't determine file type of " + filename, 0);
	}
	if (!S_ISREG(buf.st_mode))
		LOG.error("Does not appear to be a regular file: " + filename, 0);

	if (filename.substr(filename.size()-4) == ".bcf")
		LOG.error("Filename ends in '.bcf'. Shouldn't you be using --bcf?\n");

	if (!compressed)
	{
		if (filename.substr(filename.size()-3) == ".gz")
			LOG.error("Filename ends in '.gz'. Shouldn't you be using --gzvcf or --gzdiff?\n");
		file_tmp.open(filename.c_str(), ios::in);
		if (!file_tmp.is_open())
			LOG.error("Could not open VCF file: " + filename, 0);

		file_in = &file_tmp;
	}
	else
		open_gz();
}

void vcf_file::open_gz()
{
	gzMAX_LINE_LEN = 1024*1024;
	gz_readbuffer = new char[gzMAX_LINE_LEN];

	if (stream)
		gzfile_in = gzdopen(fileno(stdin), "r");
	else
		gzfile_in = gzopen(filename.c_str(), "rb");

	if (gzfile_in == NULL)
		LOG.error("Could not open GZVCF file: " + filename, 0);
	#ifdef ZLIB_VERNUM
		string tmp(ZLIB_VERSION);
		LOG.printLOG("Using zlib version: " + tmp + "\n");
		#if (ZLIB_VERNUM >= 0x1240)
			gzbuffer(gzfile_in, gzMAX_LINE_LEN); // Included in zlib v1.2.4 and makes things MUCH faster
		#else
			LOG.printLOG("Versions of zlib >= 1.2.4 will be *much* faster when reading zipped VCF files.\n");
		#endif
	#endif
}

void vcf_file::close()
{
	if (compressed)
	{
		gzclose(gzfile_in);
		delete [] gz_readbuffer;
	}
}

bool vcf_file::eof()
{
	bool out;
	if (!compressed)
		out = file_in->eof();
	else
		out = gzeof(gzfile_in);	// Returns 1 when EOF has previously been detected reading the given input stream, otherwise zero.
	return out;
}

void vcf_file::get_entry(vector<char> &out)
{
	out.resize(0);
	read_line(out);
}

entry* vcf_file::get_entry_object()
{
	return new vcf_entry(meta_data, include_indv);
}

void vcf_file::read_line(string &out)
{
	char * tmp;
	out = "";
	if (!compressed)
	{
		getline(*file_in, out);
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
	else
	{
		bool again = true;
		while (again == true)
		{
			tmp = gzgets(gzfile_in, gz_readbuffer, gzMAX_LINE_LEN);
			if (tmp == NULL)
				return;

			out.append(gz_readbuffer);
			if ((strlen(gz_readbuffer) != gzMAX_LINE_LEN-1) || (gz_readbuffer[gzMAX_LINE_LEN-2] == '\n'))
				again = false;
		}
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)
	}
}

void vcf_file::read_line(vector<char> &out)
{
	static string tmp;
	tmp="";
	out.resize(0);
	read_line(tmp);
	vector<char> tmp_char(tmp.begin(),tmp.end());
	out = tmp_char;
}
