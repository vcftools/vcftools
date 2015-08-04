/*
 * entry.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: amarcketta
 */

#include "entry.h"

/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/
void entry::SNPHWE(int obs_hets, int obs_hom1, int obs_hom2, double &p_hwe, double &p_lo, double &p_hi)
{
	p_hwe = 1.0; p_lo = 1.0; p_hi = 1.0; //p_hi_lo = 1.0;
	if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return;

	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
		LOG.error("Internal error: negative count in HWE test", 91);

	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;

	double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL)
		LOG.error("Internal error: SNP-HWE: Unable to allocate array", 90);

	for (int i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
	  het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
	  sum += het_probs[curr_hets - 2];

	  /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
	  curr_homr++;
	  curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
	  het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc	/((curr_hets + 2.0) * (curr_hets + 1.0));
	  sum += het_probs[curr_hets + 2];

	  /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
	  curr_homr--;
	  curr_homc--;
	}

	for (int i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	// alternate p-value calculation for p_hi/p_lo
	p_hi = het_probs[obs_hets];
	for (int i = obs_hets + 1; i <= rare_copies; i++)
	 p_hi += het_probs[i];

	p_lo = het_probs[obs_hets];
	for (int i = obs_hets - 1; i >= 0; i--)
	  p_lo += het_probs[i];

	//p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;

	p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (int i = 0; i <= rare_copies; i++)
	{
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	free(het_probs);
}

void entry::make_typed_string(vector<char> &out, const string &in, bool typed)
{
	vector<char> tmp_vector;
	out.resize(0);

	if (in == "." or in == " " or in == "")
	{
		if (typed == false)
			return;

		int8_t tmp = (int8_t)0;
		tmp = tmp << 4;
		tmp = tmp | (int8_t)7;
		out.push_back( tmp );

		return;
	}

	if (typed == true)
	{
		if (in.length() >= 15)
		{
			int8_t tmp = (int8_t)15;
			tmp = tmp << 4;
			tmp = tmp | (int8_t)7;
			out.push_back( tmp );

			make_typed_int(tmp_vector, in.length(), typed);
			out.insert( out.end(), tmp_vector.begin(), tmp_vector.end() );
		}
		else
		{
			int8_t tmp = (int8_t)in.length();
			tmp = tmp << 4;
			tmp = tmp | (int8_t)7;
			out.push_back( tmp );
		}
	}
	out.reserve(out.size()+in.size());
	copy(in.begin(), in.end(), back_inserter(out));
}

void entry::make_typed_int(vector<char> &out, const int &in, bool typed)
{
	vector<char> tmp_char;
	out.resize(0);

	int type;
	int8_t size_type = (int8_t)1;
	if (in < 127 and in >-127)
		type = 1;
	else if (in < 32767 and in>-32767)
		type = 2;
	else
		type = 3;

	make_int(tmp_char, in, type);

	if (typed == true)
	{
		size_type = size_type << 4;
		size_type = size_type | type;
		out.push_back(size_type);
	}
	out.insert(out.end(), tmp_char.begin(), tmp_char.end());
}

void entry::make_typed_string_vector( vector<char> &out, const vector<string> &in, int number )
{
	vector<char> tmp_char;
	int max_val = 0;
	int8_t size_type;
	out.resize(0);

	if (number == -1)
	{
		for (unsigned int ui=0; ui<in.size(); ui++)
		{
			if ((int)in[ui].size() > max_val)
				max_val = in[ui].size();
		}
	}
	else
		max_val = number;

	if (max_val < 15)
	{
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)7;

		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)7;

		out.push_back( size_type );

		make_typed_int(tmp_char, max_val, true);
		out.insert( out.end(), tmp_char.begin(), tmp_char.end() );
	}

	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		for (unsigned int uj=0; (int)uj<max_val; uj++)
		{
			if (in[ui] == ".")
				out.push_back( '\0' );
			else if (uj<in[ui].size())
				out.push_back( in[ui][uj] );
			else
				out.push_back( '\0' );
		}
	}
}

void entry::make_typed_GT_vector(vector<char> &out, vector<string> &in )
{
	vector<char> tmp_vector;
	int8_t size_type;
	int max_ploidy = 0;
	out.resize(0);

	max_ploidy = *max_element(ploidy.begin(), ploidy.end());

	if (max_ploidy < 15)
	{
		size_type = (int8_t)max_ploidy;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back( size_type );

		make_typed_int(tmp_vector, max_ploidy, true);
		out.insert( out.end(), tmp_vector.begin(), tmp_vector.end() );
		tmp_vector.resize(0);
	}

	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		encode_genotype(tmp_vector, in[ui], max_ploidy);
		out.insert( out.end(), tmp_vector.begin(), tmp_vector.end() );
		tmp_vector.resize(0);
	}
	out.insert( out.end(), tmp_vector.begin(), tmp_vector.end() );
}

void entry::encode_genotype(vector<char> &out, string &in, int exp_size)
{
	int8_t tmp_int;
	int8_t phased = 0;
	out.resize(exp_size);
	int idx = 0;

	for (unsigned int ui=0; ui<in.length(); ui++)
	{
		if (in[ui] =='|')
			phased = 1;
		else if (in[ui] == '/')
			phased = 0;
		else
		{
			if(in[ui] != '.')
			{
				tmp_int = header::str2int( in.substr(ui,1) );
				tmp_int++;
				tmp_int = tmp_int << 1;
				tmp_int = tmp_int | phased;
			}
			else
				tmp_int = 0x80;

			out[idx] = (int8_t)tmp_int;
			idx++;
		}
	}
	while (idx<exp_size)
	{
		out[idx] = (int8_t)0x81;
		idx++;
	}
}

void entry::make_typed_int_vector(vector<char> &out, const string &in, int number )
{
	vector<char> tmp_char;
	vector<int> tmp_ints;
	vector<string> split_string;
	int converted, type;
	int8_t size_type;
	unsigned int max = 0;
	unsigned int max_val = 0;
	out.resize(0);

	if (in == " " or in == "." or in == "")
	{
		size_type = (int8_t)0;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back( size_type );
		return;
	}

	header::tokenize(in, ',', split_string);
	if (number == -1)
	{
		if (split_string.size() > max_val)
			max_val = split_string.size();
	}
	else
		max_val = number;

	for (unsigned int ui=0; ui<max_val; ui++)
	{
		if ( ui<split_string.size() )
		{
			converted = header::str2int(split_string[ui], 0x80000000);

			if ((abs(converted) > (int)max) and ( converted != (int)0x80000000))
				max = abs(converted);
		}
		else
			converted = 0x80000001;

		tmp_ints.push_back( converted );
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;

	if (max_val < 15)
	{
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );

		make_typed_int(tmp_char, max_val, true);
		out.insert( out.end(), tmp_char.begin(), tmp_char.begin() );
	}

	for (unsigned int ui=0; ui<tmp_ints.size(); ui++)
	{
		make_int(tmp_char, tmp_ints[ui], type);
		out.insert( out.end(), tmp_char.begin(), tmp_char.end() );
	}
}

void entry::make_typed_int_vector(vector<char> &out, const vector<string> &in, int number )
{
	vector<char> tmp_char;
	vector<int> tmp_ints;
	vector<string> split_string;
	int converted, type;
	int8_t size_type;
	unsigned int max = 0;
	unsigned int max_val = 0;
	out.resize(0);

	if (number == -1)
	{
		unsigned int tmp_int = 0;
		for (unsigned int ui=0; ui<in.size(); ui++)
		{
			tmp_int = count(in[ui].begin(), in[ui].end(), ',');

			if (tmp_int > max_val)
				max_val = tmp_int;
		}
		max_val++;
	}
	else
		max_val = number;

	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		header::tokenize(in[ui], ',', split_string);
		for (unsigned int uj=0; uj<max_val; uj++)
		{
			if ( uj<split_string.size() )
			{
				converted = header::str2int(split_string[uj], 0x80000000);

				if ((abs(converted) > (int)max) and (converted != (int)0x80000000))
					max = abs(converted);
			}
			else
				converted = 0x80000001;

			tmp_ints.push_back( converted );
		}
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;
	if (max_val < 15)
	{
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );

		make_typed_int(tmp_char, max_val, true);
		out.insert( out.end(), tmp_char.begin(), tmp_char.begin() );
	}

	for (unsigned int ui=0; ui<tmp_ints.size(); ui++)
	{
		make_int(tmp_char, tmp_ints[ui], type);
		out.insert( out.end(), tmp_char.begin(), tmp_char.end() );
	}
}

void entry::make_typed_int_vector(vector<char> &out, const vector<int> &in )
{
	vector<char> tmp_char;
	int type;
	int8_t size_type;
	unsigned int max = 0;
	out.resize(0);

	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		if ((abs(in[ui]) > (int)max) and ( (int8_t)in[ui] != (int8_t)0x80))
			max = abs(in[ui]);
	}

	if (max < 127)
		type = 1;
	else if (max < 32767)
		type = 2;
	else
		type = 3;

	if (in.size() < 15)
	{
		size_type = (int8_t)in.size();
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)type;

		out.push_back( size_type );

		make_typed_int(tmp_char, in.size(), true);
		out.insert( out.end(), tmp_char.begin(), tmp_char.begin() );
	}

	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		make_int(tmp_char, in[ui], type);
		out.insert( out.end(), tmp_char.begin(), tmp_char.end() );
	}
}

void entry::make_int(vector<char> &out, const int &in, int type)
{
	out.resize(0);
	if (type == 1)
	{
		int8_t tmp_int;
		if (in == (int)0x80000000 || in >= 128)
			tmp_int = (int8_t)0x80;
		else if (in == (int)0x80000001)
			tmp_int = (int8_t)0x81;
		else
			tmp_int = (int8_t)in;
		out.push_back( (int8_t)tmp_int);
	}
	else if (type == 2)
	{
		int16_t tmp_int;

		if (in == (int)0x80000000 || in >= 32768)
			tmp_int = 0x8000;
		else if (in == (int)0x80000001)
			tmp_int = (int8_t)0x8001;
		else
			tmp_int = (int16_t)in;

		int8_t split;
		for(unsigned int ui=0; ui<2; ui++)
		{
			split = tmp_int & (int16_t)0x00FF;//0000000011111111
			out.push_back(split);
			tmp_int = tmp_int >> 8;
		}
	}
	else
	{
		int32_t tmp_int;
		tmp_int = (int32_t)in;

		int8_t split;
		for(unsigned int ui=0; ui<4; ui++)
		{
			split = tmp_int & (int32_t)0x0000FF;
			out.push_back( (int8_t)split);
			tmp_int = tmp_int >> 8;
		}
	}
}

void entry::make_typed_float_vector(vector<char> &out, const string &in, int number )
{
	vector<string> split_string;
	int8_t size_type;
	int max_val = 0;
	out.resize(0);

	if (in == " " or in == "." or in == "")
	{
		size_type = (int8_t)0;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)1;

		out.push_back( size_type );
		return;
	}

	header::tokenize(in, ',', split_string);
	if (number == -1)
		max_val = split_string.size();
	else
		max_val = number;

	if ( max_val < 15 )
	{
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back( size_type );

		vector<char> size_vector;
		make_typed_int(size_vector, max_val, true );
		out.insert(out.end(), size_vector.begin(), size_vector.end());
	}

	float value;
	char missing[4] = {static_cast<char>(0x01), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};
	char end[4] = {static_cast<char>(0x02), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};

	for(unsigned int ui=0; (int)ui<max_val; ui++)
	{
		if (ui < split_string.size() )
			value = (float)header::str2double(split_string[ui], 0x7F800001);
		else
			value = 0x7F800002;

		char *p = (char *)&value;

		for(unsigned int uj=0; uj<sizeof(value); uj++)
			if (value == (float)0x7F800001)
				out.push_back( missing[uj] );
			else if (value == (float)0x7F800002)
				out.push_back( end[uj] );
			else
				out.push_back( p[uj] );
	}
}

void entry::make_typed_float_vector(vector<char> &out, const vector<string> &in, int number )
{
	vector<string> split_string;
	int8_t size_type;
	unsigned int max_val = 0;
	out.resize(0);

	if (number == -1)
	{
		unsigned int tmp_int = 0;
		for (unsigned int ui=0; ui<in.size(); ui++)
		{
			tmp_int = count(in[ui].begin(), in[ui].end(), ',');
			if (tmp_int > max_val)
				max_val = tmp_int;
		}
		max_val++;
	}
	else
		max_val = number;

	if ( max_val < 15 )
	{
		size_type = (int8_t)max_val;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back( size_type );
	}
	else
	{
		size_type = (int8_t)15;
		size_type = size_type << 4;
		size_type = size_type | (int8_t)5;
		out.push_back( size_type );

		vector<char> size_vector;
		make_typed_int(size_vector, max_val, true );
		out.insert(out.end(), size_vector.begin(), size_vector.end());
	}

	float value;
	char missing[4] = {static_cast<char>(0x01), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};
	char end[4] = {static_cast<char>(0x02), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};
	for (unsigned int ui=0; ui<in.size(); ui++)
	{
		header::tokenize(in[ui], ',', split_string);
		for(unsigned int uj=0; uj<max_val; uj++)
		{
			if (uj < split_string.size() )
				value = (float)header::str2double(split_string[uj], 0x7F800001);
			else
				value = (float)0x7F800002;

			char *p = (char *)&value;

			for(unsigned int uk=0; uk<sizeof(value); uk++)
			{
				if (value == float(0x7F800001))
					out.push_back( missing[uk] );
				else if (value == float(0x7F800002))
					out.push_back( end[uk] );
				else
					out.push_back( p[uk] );
			}
		}
	}
}

void entry::make_type_size(vector<char> &out, const unsigned int &type, const unsigned int &size)
{
	uint8_t byte;
	vector<char> tmp_vector;
	tmp_vector.resize(0);
	out.resize(0);

	if (size < 15)
	{
		byte = size;
		byte = byte << 4;
	}
	else
	{
		byte = (uint8_t)15;
		make_typed_int(tmp_vector, size, true);
	}

	byte = byte | (uint8_t)type;
	out.push_back(byte);
	out.insert(out.end(), tmp_vector.begin(), tmp_vector.end());
}

string entry::get_typed_string(unsigned int * line_position, const vector<char>& line)
{
	unsigned int size, type;
	string out;

	get_type( line_position, line, type, size );
	if (type != 7)
	{
		LOG.printLOG("Error: Expected type 7 for string. Found type " + header::int2str(type) + ".\n");
	}

	char * tmp = new char[size];
	memcpy(tmp, &line[*line_position], size*sizeof(char));
	*line_position += size;
	out = string( tmp, size );

	if (out == "" or out == " ")
		out = ".";

	delete [] tmp;
	return out;
}

int entry::get_typed_int(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size)
{
	int out;

	get_type( line_position, line, type, size );

	if (size > 1)
	{
		LOG.printLOG("Error: Int vector when expected only a single Integer value.\n" );
		exit(0);
	}

	if (type == 1)
	{
		int8_t tmp;
		tmp = *reinterpret_cast<const int8_t*>(&line[*line_position]);
		*line_position += sizeof(tmp);
		out = tmp;
	}
	else if (type == 2)
	{
		int16_t tmp;
		tmp = *reinterpret_cast<const int16_t*>(&line[*line_position]);
		*line_position += sizeof(tmp);
		out = tmp;
	}
	else if (type == 3)
	{
		int32_t tmp;
		tmp = *reinterpret_cast<const int32_t*>(&line[*line_position]);
		*line_position += sizeof(tmp);
		out = tmp;
	}
	else
	{
		LOG.printLOG("Error: Invalid type for integer size.\n");
		exit(0);
	}
	return out;
}

vector<int> entry::get_int_vector(unsigned int * line_position, const vector<char>& line)
{
	unsigned int size, type;
	get_type( line_position, line, type, size );
	vector<int> out(size);

	if (type == 0)
	{
		return out;
	}
	else if (type == 1)
	{
		int8_t tmp;
		for (unsigned int ui=0; ui<size; ui++)
		{
			tmp = *reinterpret_cast<const int8_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);
			out[ui] = tmp;
		}
	}
	else if (type == 2)
	{
		int16_t tmp;
		for (unsigned int ui=0; ui<size; ui++)
		{
			tmp = *reinterpret_cast<const int16_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);
			out[ui] = tmp;
		}
	}
	else if (type == 3)
	{
		int32_t tmp;
		for (unsigned int ui=0; ui<size; ui++)
		{
			tmp = *reinterpret_cast<const int32_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);
			out[ui] = tmp;
		}
	}
	else
	{
		LOG.printLOG("Error: Invalid type for integer size.\n");
		exit(0);
	}
	return out;
}

void entry::get_type(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size)
{
	uint8_t byte = *reinterpret_cast<const uint8_t*>(&line[*line_position]);
	*line_position += sizeof(byte);
	size = byte >> 4;
	type = (byte & (uint8_t)15);

	if (size == 15)
	{
		int type2;
		byte = *reinterpret_cast<const uint8_t*>(&line[*line_position]);
		*line_position += sizeof(byte);

		type2 = (byte & (uint8_t)15);
		if (type2 == 1)
		{
			int8_t tmp;
			tmp = *reinterpret_cast<const int8_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);
			size = (unsigned int)tmp;
		}
		else if (type2 == 2)
		{
			int16_t tmp;
			tmp = *reinterpret_cast<const int16_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);

			size = (int)tmp;
		}
		else if (type2 == 3)
		{
			int32_t tmp;
			tmp = *reinterpret_cast<const int32_t*>(&line[*line_position]);
			*line_position += sizeof(tmp);
			size = (unsigned int)tmp;
		}
		else
		{
			LOG.printLOG("Error: Invalid type for integer size.\n");
			exit(0);
		}
	}
}

void entry::skip_section(unsigned int *line_position, const vector<char> &line)
{
	unsigned int type, size;
	get_type(line_position, line, type, size);

	if ( (type == 1) || (type == 7) )
		*line_position += sizeof(int8_t)*size;
	else if (type == 2)
		*line_position += sizeof(int16_t)*size;
	else if ( (type == 3) || (type == 5) )
		*line_position += sizeof(int32_t)*size;
}

bool entry::check_missing(unsigned int line_position, const unsigned int type, const vector<char> &line)
{
	static char missing_float[4] = {static_cast<char>(0x01), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};
	static char missing_int1 = static_cast<char>(0x80);
	static char missing_int2[2] = {static_cast<char>(0x00), static_cast<char>(0x80)};
	static char missing_int3[4] = {static_cast<char>(0x00), static_cast<char>(0x00), static_cast<char>(0x00), static_cast<char>(0x80)};

	char test_char;
	bool missing = true;
	if (type==1)
	{
		test_char = *reinterpret_cast<const char*>(&line[line_position]);
		missing = (test_char == missing_int1);
	}
	else if (type==2)
	{
		for (unsigned int ui=0; ui<sizeof(int16_t); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != missing_int2[ui])
			{
				missing = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==3)
	{
		for (unsigned int ui=0; ui<sizeof(int32_t); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != missing_int3[ui])
			{
				missing = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==5)
	{
		for (unsigned int ui=0; ui<sizeof(float); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != missing_float[ui])
			{
				missing = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==7)
		missing = false;

	return missing;
}

bool entry::check_end(unsigned int line_position, const unsigned int type, const vector<char> &line)
{
	static char end_float[4] = {static_cast<char>(0x02), static_cast<char>(0x00), static_cast<char>(0x80), static_cast<char>(0x7F)};
	static char end_int1 = static_cast<char>(0x81);
	static char end_int2[2] = {static_cast<char>(0x01), static_cast<char>(0x80)};
	static char end_int3[4] = {static_cast<char>(0x01), static_cast<char>(0x00), static_cast<char>(0x00), static_cast<char>(0x80)};

	char test_char;
	bool end = true;
	if (type==1)
	{
		test_char = *reinterpret_cast<const char*>(&line[line_position]);
		end = (test_char == end_int1);
	}
	else if (type==2)
	{
		for (unsigned int ui=0; ui<sizeof(int16_t); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != end_int2[ui])
			{
				end = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==3)
	{
		for (unsigned int ui=0; ui<sizeof(int32_t); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != end_int3[ui])
			{
				end = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==5)
	{
		for (unsigned int ui=0; ui<sizeof(float); ui++)
		{
			test_char = *reinterpret_cast<const char*>(&line[line_position]);
			if (test_char != end_float[ui])
			{
				end = false;
				break;
			}
			line_position += sizeof(char);
		}
	}
	else if (type==7)
		end = false;

	return end;
}

void entry::get_number(uint32_t &out, unsigned int *line_position, const vector<char>& line)
{
	memcpy(&out, &line[*line_position], sizeof(out));
	*line_position += sizeof(out);
}
