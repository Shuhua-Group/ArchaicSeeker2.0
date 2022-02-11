/*
 * Utils.cpp
 *
 *  Created on: Aug 26, 2015
 *      Author: young
 */
#include <cmath>
#include <iostream>
#include "EMExp.hpp"
#include "Utils.hpp"
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

double cv_chisq(int df, double alpha)
{
	boost::math::chi_squared dist(df);
	return boost::math::quantile(dist, 1 - alpha);
}

ParamExp findOptPar(const vector<double> &observ, int K, int maxIter, double ancestryProp, \
		double criticalValue, double epsilon, double minP, bool simple)
{
//	bool findOpt = false;
//	int k = 1;
//	ParamExp parPrev(K);
//	EMExp em(parPrev, observ);
//	em.iterateClassify(maxIter, epsilon);
//	double llkPrev = em.getLik();
//	parPrev = em.getPar();
//	if (simple)
//	{
//		parPrev.sortByLambda();
//		return parPrev;
//	}
//	while (!findOpt)
//	{
//		k++;
		ParamExp parCur(K);
		EMExp em(parCur, observ);
		em.setPar(parCur);
		//em.updateLik();
		em.iterate(maxIter, epsilon);
		//cout << "K=" << k << "; Likelihood=" << setprecision(8) << em.getLik() << "; ";
		//em.getPar().print();
//		double llkCur = em.getLik();
//		if (2 * (llkCur - llkPrev) < criticalValue)
//		{
			//cout << "Optimal K " << k - 1 << endl;
//			findOpt = true;
//		}
//		else
//		{
			/*
			 * check survival proportion for each wave,
			 * if any one less than minP, stop to increase to larger K
			 */
/*			ParamExp tmpPar = em.getPar();
			tmpPar.sortByLambda();
			double tempSum = 0;
			double temp[k];
			for (int i = 0; i < k; ++i)
			{
				temp[i] = tmpPar.getProp(i) / tmpPar.getLambda(i);
				tempSum += temp[i] / temp[0];
			}
			tempSum = ancestryProp / tempSum;
			for (int i = 0; i < k; ++i)
			{
				if (tempSum * temp[i] / temp[0] < minP)
				{
					findOpt = true;
					break;
				}
			}
			if (!findOpt)
			{
				llkPrev = llkCur;
				parPrev = em.getPar();
			}
		}
	}
	//parPrev.print();
	 */
//	parPrev.sortByLambda();
	em.getPar().sortByLambda();
	return em.getPar();
}

ParamExp findOptPar(const vector<double> &observ, int maxIter, double ancestryProp, \
		double criticalValue, double epsilon, double minP, bool simple)
{
	bool findOpt = false;
	int k = 1;
	ParamExp parPrev(k);
	EMExp em(parPrev, observ);
	em.iterate(maxIter, epsilon);
	double llkPrev = em.getLik();
	parPrev = em.getPar();
	if (simple)
	{
		parPrev.sortByLambda();
		return parPrev;
	}
	while (!findOpt)
	{
		k++;
		ParamExp parCur(k);
		em.setPar(parCur);
		//em.updateLik();
		em.iterate(maxIter, epsilon);
		//cout << "K=" << k << "; Likelihood=" << setprecision(8) << em.getLik() << "; ";
		//em.getPar().print();
		double llkCur = em.getLik();
		if (2 * (llkCur - llkPrev) < criticalValue)
		{
			//cout << "Optimal K " << k - 1 << endl;
			findOpt = true;
		}
		else
		{
			/*
			 * check survival proportion for each wave,
			 * if any one less than minP, stop to increase to larger K
			 */
			ParamExp tmpPar = em.getPar();
			tmpPar.sortByLambda();
			double tempSum = 0;
			double temp[k];
			for (int i = 0; i < k; ++i)
			{
				temp[i] = tmpPar.getProp(i) / tmpPar.getLambda(i);
				tempSum += temp[i] / temp[0];
			}
			tempSum = ancestryProp / tempSum;
			for (int i = 0; i < k; ++i)
			{
				if (tempSum * temp[i] / temp[0] < minP)
				{
					findOpt = true;
					break;
				}
			}
			if (!findOpt)
			{
				llkPrev = llkCur;
				parPrev = em.getPar();
			}
		}
	}
	//parPrev.print();
	parPrev.sortByLambda();
	return parPrev;
}

void solveTrueProp(ParamExp &par, double lower)
{
	int numOfWave = par.getK();
	double temp[numOfWave];
	temp[0] = 1.0;
	double tempSum = temp[0];
	for (int i = 1; i < numOfWave; ++i)
	{
		temp[i] = par.getProp(i) * exp((par.getLambda(i) - par.getLambda(0)) * lower) / par.getProp(0);
		tempSum += temp[i];
	}
	par.setProp(0, 1.0 / tempSum);
	for(int i = 1; i < numOfWave; ++i)
	{
		par.setProp(i, par.getProp(0) * temp[i]);
	}
}

void help()
{
	cout << kProgramName << " v" << kVersion << endl;
	cout << kProgramName << " is designed to scan for multiple waves of admixture events and estimated corresponding parameters." << endl;
	cout << "General usage: " << kProgramName << " <arguments> [options]" << endl;
	cout << "Arguments" << endl;
	cout << "\t-i/--input\t<string>\tInput of the ancestral tracks [required]" << endl;
	cout << "\t-l/--lower\t[double]\tLower bound to discard short tracks [optional, default 0]" << endl;
	cout << "\t-b/--bootstrap\t[integer]\tNumber of bootstrapping [optional, default 100]" << endl;
	cout << "\t-a/--alpha\t[double]\tSignificance level to reject null hypothesis in LRT [optional, default 0.001]" << endl;
	cout << "\t-e/--epsilon\t[double]\tEpsilon to check whether a parameter converge or not [optional, default 1.0e-6] " << endl;
	cout << "\t-p/--minProp\t[double]\tMinimum survival proportion for a wave at the final generation [optional, default 0.05]" << endl;
	cout << "\t-m/--maxIter\t[integer]\tMaximum number of iterations to scan for waves of admixture events [optional, default 10000]" << endl;
	cout << "\t-t/--thread\t[integer]\tNumber of threads [optional, default 1]" << endl;
	cout << "\t-o/--output\t<string>\tPrefix of output [required]" << endl;
	cout << "Options" << endl;
	cout << "\t-h/--help\tPrint help message, default is OFF" << endl;
}

