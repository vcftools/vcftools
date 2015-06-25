/*
 * vcf_file.h
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#ifndef BCF_FILE_H_
#define BCF_FILE_H_

#include "output_log.h"
#include "parameters.h"
#include "variant_file.h"
#include "bgzf.h"

extern output_log LOG;
using namespace std;

class bcf_file : public variant_file
{
public:
	bcf_file(const parameters &params, bool diff=false);

	void get_entry(vector<char> &out);
	entry* get_entry_object();

	void print(const parameters &params);
	void print_bcf(const parameters &params);

protected:
	~bcf_file();

private:
	bool is_BGZF;
	bool big_endian;
	bool stream;

	int read(void *buffer, unsigned int len, size_t size);
	void read_header();
	void read_file();
	void open();
	void open_gz();
	void close();
	bool eof();
	void check_bcf();
};

#endif /* BCF_FILE_H_ */
