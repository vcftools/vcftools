/*
	entry_setters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "entry.h"

void entry::add_ALT_allele(const string &in)
{
	if (in != ".")
	{
		if (find(ALT.begin(), ALT.end(),in) == ALT.end())
		{
			ALT.push_back(in);
		}
		else
		{
			LOG.warning(" Duplicate alternate alleles found at " + CHROM + ":" + LOG.int2str(POS));
			ALT.push_back(in);
		}
	}
	parsed_ALT = true;
}

void entry::set_indv_DEPTH(unsigned int indv, int in)
{
	parsed_DP[indv] = true;
	if (in == -1)
	{
		if (!DEPTH.empty())
			DEPTH[indv] = -1;
		return;
	}
	if (DEPTH.empty())
		DEPTH.resize(N_indv, -1);

	DEPTH[indv] = in;
}
