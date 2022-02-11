/*
 * EMExp.hpp
 *
 *  Created on: May 26, 2015
 *      Author: young
 */
/******************************************************************************
 * @brief A class perform EM
 * Class EMExp to perform EM algorithm in parameters estimation with hidden variable
 * The EMExp class has three attributes:
 * 1) parameters to be estimated;
 * 2) a sequence of observations, which used to estimate the parameters;
 * 3) and the likelihood corresponding to the observations and parameters.

 * The method iterate is used to perform EM algorithm, including two steps:
 * E-Step:
 * 1) calculate p_j=prob(z_i=j|x_i, theta_t)
 * =prob(x_i|z_i=j,theta_t)*prob(z_i=j,theta_t)/(sum_j{from 1 to K}(.)
 * denotes the numerator
 * 2) calculate p_j*xi ; j=1,...,K
 * theta_t: parameters at the t-th iteration
 * M-Step:
 * 1) update prop, refer as m, m_j=sum(p_j)/sum(p_1+p_2+...+p_K)=sum(p_j)/n
 * 2) update lambda, lambda_j=sum(p_j)/sum(p_j*x_i)
 * note: update the parameters for t+1 times iteration
 ******************************************************************************/
#ifndef EMEXP_HPP_
#define EMEXP_HPP_

#include <vector>
#include "ParamExp.hpp"

/*
 * summation over data, with certain size
 */
double sum(double *data, int size);

double modernWaver(unsigned int n, const double *x, double *grad, void *data);

double archaicWaver(unsigned int n, const double *x, double *grad, void *data);

double archaicWaverConstraint(unsigned n, const double *x, double *grad, void *data);

class waverData
{
public:
	waverData(const std::vector<double>& _observ, double _lb, double _ub) : \
			lb(_lb), ub(_ub), observ(_observ) {};
	double lb, ub;
	const std::vector<double>& observ;
};

class EMExp
{
public:

	/*
	 * @brief constructor
	 * @param par Parameter for EM algorithm
	 * @param observ Observations of data used to estimate parameter of EM
	 */
	EMExp(const ParamExp &par, const std::vector<double> &observ);

    /*
	 * @brief get the log likelihood
	 * @return the log likelihood of current data and parameter
	 */
	double getLik() const;

    /*
	 * @brief get parameter of EM
	 * @return the reference of EM parameter
	 */
	ParamExp & getPar();

    /*
	 * @brief set a new parameter for EM
	 * @param new parameter to be set
	 */
	void setPar(const ParamExp &par);

    /*
	 * @brief update the log likelihood, usually when parameter is updated
	 */
	void updateLik();

    /*
	 * @brief perform EM iteration
	 * @param maxIter max number of iteration, the EM iteration terminate either
	 * parameter is converged, or reach the max number of iteration
	 */
	void iterate(int maxIter, double epsilon);

	void iterateClassify(int maxIter, double epsilon);

    /*
	 * destructor
	 */
	virtual ~EMExp();

private:
    /*
     * log likelihood value
     */
	double lik;

    /*
     * Parameter of EM
     */
	ParamExp par;

    /*
     * a vector of observations
     */
	std::vector<double> observ;
};

#endif /* EMEXP_HPP_ */
