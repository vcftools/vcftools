/*
 * vcf_file.h
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#ifndef VCF_FILE_H_
#define VCF_FILE_H_

#include "output_log.h"
#include "vcf_entry.h"
#include "parameters.h"
#include "variant_file.h"

extern output_log LOG;

using namespace std;

class vcf_file : public variant_file
{
public:
	vcf_file(const parameters &params, bool diff=false);

	void get_entry(vector<char> &out);
	entry* get_entry_object();

	void print(const parameters &params);
	void print_bcf(const parameters &params);

protected:
	~vcf_file();

private:
	char *gz_readbuffer;
	void open();
	void open_gz();
	void close();
	bool eof();
	inline void read_line(string &out);
	inline void read_line(vector<char> &out);
	void read_header();

	bool stream;
};

#endif /* VCF_FILE_H_ */
