/*
 * log.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 66 $)
 */

#include "output_log.h"

output_log::output_log()
{
	output_to_screen = true;
	output_to_file = true;
}

void output_log::open(bool stout, bool sterr, const string &filename_prefix)
{
	if (stout)
	{
		output_to_screen = false;
		output_to_file = true;
	}

	if (sterr)
	{
		output_to_screen = true;
		output_to_file = false;
	}

	if (output_to_file)
		LOG.open((filename_prefix + ".log").c_str(), ios::out);
}

void output_log::close()
{
	LOG.close();
}

void output_log::printLOG(string s)
{
	if (output_to_file)
		LOG << s; LOG.flush();
	if (output_to_screen)
		cerr << s; cerr.flush();
}

void output_log::error(string err_msg, int error_code)
{
	printLOG("Error: " + err_msg + "\n");
	exit(error_code);
}


void output_log::error(string err_msg, double value1, double value2, int error_code)
{
	printLOG("Error: " + err_msg + "\n");
	stringstream ss;
	ss << "Value1=" << value1 << " Value2=" << value2 << endl;
	printLOG(ss.str());
	exit(error_code);
}

void output_log::warning(string err_msg)
{
	printLOG(err_msg + "\n");
}

void output_log::one_off_warning(string err_msg)
{
	static set<string> previous_warnings;
	if (previous_warnings.find(err_msg) == previous_warnings.end())
	{
		printLOG(err_msg + "\n");
		previous_warnings.insert(err_msg);
	}
}

string output_log::int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string output_log::longint2str(long int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

string output_log::dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

string output_log::dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << setiosflags( ios::fixed );
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}
