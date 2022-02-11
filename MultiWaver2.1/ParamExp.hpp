/*
 * ParamExp.hpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */
/******************************************************************************
 * @brief Parameters of EM
 * Class ParamExp to store the parameters of EM
 * Here assume the observations are from the combination of K sets of numbers
 * which following exponential distributions, with parameters lambda[j], and the
 * mix proportions prop[j], the purpose of EM is used to estimate the parameters
 * lambda[j] and mix proportions[j], where j in [1,2,...,K]
 *****************************************************************************/

#ifndef PARAMEXP_HPP_
#define PARAMEXP_HPP_

//const double kDelta = 0.000001; //converge condition

# include <vector>

class ParamExp
{
public:

	/*
	 * @brief Constructor with given value K
	 * @param K number of exponential distribution
	 */
	ParamExp(int K = 1);

	ParamExp(int K, std::vector<double> observ);

    /*
	 * @brief Constructor with given value K, lambdas and proportions
	 * @param K number of exponential distribution
	 * @param lambda initial value of lambda
	 * @param prop initial value of proportions
	 */
	ParamExp(int K, double *lambda, double *prop);

    /*
	 * @brief Copy constructor
	 * @param rhs	old parameter
	 */
	ParamExp(const ParamExp &rhs);

    /*
	 * @brief overloading of operator =
	 * @param rhs	old parameter
	 */
	ParamExp & operator=(const ParamExp &rhs);

    /*
	 * @brief get the number of exponential distributions
	 * @return integer number
	 */
	int getK() const;

    /*
	 * @brief get the value of lambda with given index
	 * @param index the exponential distribution
	 * @return the value of lambda
	 */
	double getLambda(int index) const;

    /*
	 * @brief get the value of proportion with given index
	 * @param index the exponential distribution
	 * @return the value of proportion
	 */
	double getProp(int index) const;
//	void setK(int);
//	void setLambda(int index, double lambda);

    /*
	 * set the Proportion at index
	 */
	void setProp(int index, double prop);

    /*
	 * @brief check if the parameter is converged
	 * @param par old parameter
	 * @param epsilon epsilon to check whether a parameter converge or not
	 * @return true if is converged
	 */
	bool isConverge(const ParamExp &par, double epsilon = 0.000001);

    /*
	 * @brief sort the parameter by lambda in ascend order
	 */
	void sortByLambda();

    /*
	 * @brief print the value of parameters
	 */
	void print();

    /*
	 * @brief destructor
	 */
	virtual ~ParamExp();

private:

    /*
	 * number of exponential distributions
	 */
	int K;

    /*
	 * parameter lambdas for exponential distributions
	 */
	double *lambda;

    /*
	 * mix proportions for exponential distributions
	 */
	double *prop;
};

#endif /* PARAMEXP_HPP_ */
