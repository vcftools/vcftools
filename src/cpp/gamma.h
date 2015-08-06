#ifndef GAMMA_H
#define GAMMA_H


#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>

using namespace std;

double gammln(double xx);
double gcf(double a, double x, double &gln);
double gser(double a, double x, double &gln);
double gammp(double a, double x);
double gammq(double a, double x);

#endif
