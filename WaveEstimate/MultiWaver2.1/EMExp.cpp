/*
 * EMExp.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#include <cmath>
#include <iostream>
#include <cstring>
#include "EMExp.hpp"

double sum(double * data, int size)
{
	double tmp = 0;
	for (int i = 0; i < size; ++i)
	{
		tmp += data[i];
	}
	return tmp;
}

double modernWaver(unsigned int n, const double *x, double *grad, void *data)
{
	waverData *d = (waverData * ) data;
	double u (x[0]);
	double expLb(exp( - u * d->lb ) ), expUb(exp( - u * d->ub ) );
	double a ( log( u / ( expLb - expUb ) ) );
	int nobserv(d->observ.size());
	double b (0);
	const std::vector<double>& curObserv(d->observ);
	for(int i = 0 ; i < nobserv; ++i)
		b += curObserv[i];
	if(grad)
	{
		grad[0] = 1.0 / u + (d->lb * expLb - d->ub * expUb) / (expLb - expUb) - b / nobserv;
	}
//	std::cout << grad[0] << std::endl;
	double ret(a * nobserv - u * b);
	return ret;
}

double archaicWaver(unsigned int n, const double *x, double *grad, void *data)
{
	//x[0 .. n / 2] u
	//x[n / 2 .. n] w
	int nwave(n/2);
	const double *u(x) , *w(x + nwave);
	double *ugrad, *wgrad;
	if(grad)
	{
		ugrad = grad;
		wgrad = grad + nwave;
	}
	waverData *d = (waverData * ) data;
	double lb(d->lb), ub(d->ub);
	double denom(0);
	double dev_euc[nwave];
	double dev_ceuc[nwave];
	for(int i = 0 ; i < nwave ; ++i)
	{
		double euc1, euc2;
		euc1 = exp(-u[i] * lb);
		euc2 = exp(-u[i] * ub);
		dev_euc[i] = euc1 - euc2;
		denom += w[i] * dev_euc[i];
		dev_ceuc[i] = lb * euc1 - ub * euc2;
	}
	const std::vector<double>& curObserv(d->observ);
	int nobserv(curObserv.size());
	double ret(0);
	double eux[nwave], ueux[nwave];
	if(grad)
	{
		memset(grad, 0, sizeof(double) * nwave * 2);
		for(int i = 0 ; i < nobserv; ++i)
		{
			double x = curObserv[i];
			double sum_wueux(0);
			for(int j = 0 ; j < nwave ; ++j)
			{
				eux[j] = exp(-u[j] * x);
				ueux[j] = eux[j] * u[j];
				sum_wueux += ueux[j] * w[j];
			}
			ret += log(sum_wueux / denom);
			for(int j = 0 ; j < nwave; ++j)
			{
				wgrad[j] += ueux[j] / sum_wueux - dev_euc[j] / denom;
				ugrad[j] += w[j] * ( eux[j] - u[j] * ueux[j] ) / sum_wueux \
						+ w[j] * dev_ceuc[j] / denom;
			}
		}
	}
	else
		for(int i = 0 ; i < nobserv; ++i)
		{
			double x = curObserv[i];
			double sum_wueux(0);
			for(int j = 0 ; j < nwave ; ++j)
				sum_wueux += exp(-u[j] * x) * u[j] * w[j];
			ret += log(sum_wueux / denom);
		}
//	std::cout << "Iter" << std::endl;
//	for(int i = 0 ; i < nwave ; ++i)
//		std::cout << '\t' << u[i] << '\t' << w[i] << std::endl;
//	for(int i = 0 ; i < nwave ; ++i)
//		std::cout << '\t' << ugrad[i] << '\t' << wgrad[i] << std::endl;
//	std::cout << ret << std::endl;
	return ret;
}

double archaicWaverConstraint(unsigned int n, const double *x, double *grad, void *data)
{
	double ret(0);
	unsigned int nwave(n / 2);
	for(unsigned int i = nwave ; i < n ; ++i)
		ret += x[i];
	if(grad)
	{
		for(unsigned int i = 0 ; i < nwave ; ++i)
			grad[i] = 0;
		for(unsigned int i = nwave ; i < n; ++i)
			grad[i] = 1;
	}
	ret -= 1;
	return ret;
}

using namespace std;

//constructor
EMExp::EMExp(const ParamExp &par, const vector<double> &observ) :
		par(par), observ(observ)
{
	updateLik();
}

double EMExp::getLik() const
{
	return lik;
}

ParamExp & EMExp::getPar()
{
	return par;
}

void EMExp::setPar(const ParamExp &par)
{
	this->par = par;
}
//update likelihood
/*
 * Log-lik = sum(log(sum(m_j*l_j*exp(-l_j*x_i))))
 */
void EMExp::updateLik()
{
	double tmp = 0;
	for (size_t i = 0; i < observ.size(); ++i)
	{
		int k = par.getK();
		double fval = 0;
		for (int j = 0; j < k; ++j)
		{
			double m = par.getProp(j);
			double l = par.getLambda(j);
			fval += m * l * exp(-l * observ.at(i));
		}
		tmp += log(fval);
	}
	lik = tmp;
}

//EM iteration, converge or reach max iteration time will terminate
void EMExp::iterate(int maxIter, double epsilon)
{
	int it = 0; //iteration number
	int kval = par.getK();
	int size = observ.size();
	double *nlambda; //new lambda
	double *nprop; //new proportion
	double **pval;
	double **pxval;
	nlambda = new double[kval];
	nprop = new double[kval];
	pval = new double *[kval];
	pxval = new double *[kval];
	for (int i = 0; i < kval; ++i)
	{
		pval[i] = new double[size];
		pxval[i] = new double[size];
	}
	while (++it < maxIter)
	{
		//E-step
		for (int i = 0; i < size; ++i)
		{
			double denormator = 0;
			for (int j = 0; j < kval; ++j)
			{
				double mj = par.getProp(j);
				double lj = par.getLambda(j);
				pval[j][i] = mj * lj * exp(-lj * observ.at(i));
				denormator += pval[j][i];
			}
			for (int j = 0; j < kval; ++j)
			{
				pval[j][i] /= denormator;
				pxval[j][i] = pval[j][i] * observ.at(i);
			}
		}
		//M-step
		for (int i = 0; i < kval; ++i)
		{
			double sump = sum(pval[i], size);
			nprop[i] = sump / size;
			nlambda[i] = sump / sum(pxval[i], size);
		}
		ParamExp updatedPar(kval, nlambda, nprop);
		/*
		 * check converge or not, if converge, jump out of loop
		 */
		if (par.isConverge(updatedPar, epsilon))
		{
			par = updatedPar;
			updateLik();
			break;
		}
		else
		{
			par = updatedPar;
			updateLik();
		}
		//cout << "Iteration " << ++it << " --> llk: " << getLik() << "; ";
		//getPar().print();
	}
//	if (it >= maxIter)
//	{
//		cerr << "Warning: Max iteration reached before convergence, set a larger number" << endl;
//	}
	//clean stuff
	delete[] nlambda;
	delete[] nprop;
	for (int i = 0; i < kval; ++i)
	{
		delete[] pval[i];
		delete[] pxval[i];
	}
	delete[] pval;
	delete[] pxval;
}

void EMExp::iterateClassify(int maxIter, double epsilon)
{
	int it = 0; //iteration number
	int kval = par.getK();
	int size = observ.size();
	double *nlambda; //new lambda
	double *nprop; //new proportion
	double *sumLen;
	int *count;
	int *lab;
	nlambda = new double[kval];
	nprop = new double[kval];
	sumLen = new double[kval];
	count = new int[kval];
	lab = new int[size];
	double sum(0);
	for(int i = 0 ; i < size; ++i)
		sum += observ[i];
	while(++it < maxIter)
	{
		for(int i = 0 ; i < size ; ++i)
		{
			double x = observ[i];
			lab[i] = 0;
			double lambda = par.getLambda(0);
			double maxLk = lambda * exp(-lambda * x);
			for(int j = 1 ; j < kval; ++j)
			{
				double lk;
				lambda = par.getLambda(j);
				lk = lambda * exp(-lambda * x);
				if(lk > maxLk)
				{
					lab[i] = j;
					maxLk = lk;
				}
			}
		}
		memset(sumLen, 0, sizeof(double) * kval);
		memset(count, 0, sizeof(int) * kval);
		for(int i = 0 ; i < size ; ++i)
		{
			sumLen[lab[i]] += observ[i];
			++count[lab[i]];
		}
		for(int i = 0 ; i < kval; ++i)
		{
			nlambda[i] = count[i] / sumLen[i];
			nprop[i] = sumLen[i] / sum;
//			std::cout << '\t' << i << '\t' << nlambda[i] << '\t' << nprop[i] << std::endl;
		}

		ParamExp updatedPar(kval, nlambda, nprop);
		/*
		 * check converge or not, if converge, jump out of loop
		 */
		if (par.isConverge(updatedPar, epsilon))
		{
			par = updatedPar;
			updateLik();
			break;
		}
		else
		{
			par = updatedPar;
			updateLik();
		}
		//cout << "Iteration " << ++it << " --> llk: " << getLik() << "; ";
		//getPar().print();
	}
	if (it >= maxIter)
	{
		cerr << "Warning: Max iteration reached before convergence, set a larger number" << endl;
	}
	//clean stuff
	delete[] nlambda;
	delete[] nprop;
	delete[] sumLen;
	delete[] count;
	delete[] lab;
}

EMExp::~EMExp()
{
}

