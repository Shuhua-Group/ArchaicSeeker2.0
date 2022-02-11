/*
 * Utils.hpp
 *
 *  Created on: Aug 26, 2015
 *      Author: young
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <algorithm>
#include <vector>
#include <string>
#include "ParamExp.hpp"

//const double kMinP = 0.05; //minimum proportion for a wave
//const double kCriticalValue = 5.991; //chi-squared value for p-value 0.05 with degree of freedom 2
const std::string kVersion = "1.0.2";
const std::string kProgramName = "MultiWaver";

/*
 * perform permutation over a sequence
 */
template<class T>
std::vector<std::vector<T> > perm(std::vector<T> &seq)
{
	std::vector<std::vector<T> > result;
	int size = seq.size();
	std::sort(seq.begin(), seq.end());
	do
	{
		if (seq.at(size - 1) == seq.at(size - 2))
			continue;
		int rsize = result.size();
		if (rsize > 0 && result.at(rsize - 1).at(size - 2) == seq.at(size - 1))
			continue;
		result.push_back(seq);
	} while (std::next_permutation(seq.begin(), seq.end()));
	return result;
}

/*
 * calculate the critical value of chi-squared distribution with significance level alpha
 */
double cv_chisq(int df, double alpha = 0.05);

/*
 * perform EM to find the best number of exponential distribution and estimate
 * corresponding parameters
 */
ParamExp findOptPar(const std::vector<double> &observ, int maxIter, double aProp, \
		double cValue = 5.991, double epsilon = 0.000001, double minP = 0.001, \
		bool simple = false);
ParamExp findOptPar(const std::vector<double> &observ, int K, int maxIter, double aProp, \
		double cValue = 5.991, double epsilon = 0.000001, double minP = 0.001, \
		bool simple = false);


/*
 * solve for proportion when considering lower bound
*/
void solveTrueProp(ParamExp &par, double lower = 0);

/*
 * print help information
 */
void help();

#endif /* UTILS_HPP_ */
