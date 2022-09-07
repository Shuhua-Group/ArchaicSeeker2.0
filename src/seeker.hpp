/*
 * seeker.hpp
 *
 *  Created on: Oct 10, 2018
 *      Author: yuankai
 */

#ifndef SEEKER_HPP_
#define SEEKER_HPP_

# include "data.hpp"
# include "gzfstream.hpp"

# include <iostream>

const long MAXINTER = 100000;

class archaicSeekerPar
{
public:
	archaicSeekerPar() : alpha(0.02), introT(2000), emit(0.99), nthreads(1), \
			outPrefix(""), seekModel("") {}
	double alpha, introT, emit;
	int nthreads;
	std::string outPrefix, seekModel;
};

void archaicSeekerHMM(archaicSeekerPar& ASpar, const genomePar& gpar);

#endif /* SEEKER_HPP_ */
