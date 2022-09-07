/*
 * ArchaicSeeker.cpp
 *
 *  Created on: Oct 10, 2018
 *      Author: yuankai
 */

# include "data.hpp"
# include "seeker.hpp"

# include <unistd.h>
# include <getopt.h>
# include <fstream>
# include <cstdlib>
# include <cassert>

static struct option archaicSeekerHMMOpt[] =
{
	{"alpha",		required_argument,	NULL,	'a'},
	{"introT",		required_argument,	NULL,	'T'},
	{"emit",		required_argument,	NULL,	'e'},
	{"out",			required_argument,	NULL,	'o'},
	{"vcf",			required_argument,	NULL,	'v'},
	{"pop",			required_argument,	NULL,	'p'},
	{"remap",		required_argument,	NULL,	'r'},
	{"model",		required_argument,	NULL,	'm'},
	{"outgroup",	required_argument,	NULL,	'X'},
	{"anc",			required_argument,	NULL,	'A'},
	{"help",		required_argument,	NULL,	'h'},
};

void archaicSeekerHMMHelp()
{
	std::cout << "ArchaicSeekerHMM v2.0" << std::endl;
	std::cout << " -a / --alpha [0.02]      Introgression proportion." << std::endl;
	std::cout << " -T / --introT [2000]     Introgression time(in generation)." << std::endl;
	std::cout << " -e / --emit [0.99]       Emission probability parameter." << std::endl;
	std::cout << " -o / --out               Output prefix." << std::endl;
	std::cout << " -v / --vcf               VCF configuration file" << std::endl;
	std::cout << " -p / --pop               Population annotation file." << std::endl;
	std::cout << " -r / --remap             Recombination configuration file." << std::endl;
	std::cout << " -X / --outgroup			Outgroup genome configuration file." << std::endl;
	std::cout << " -A / --anc				Ancestral state configuration file." << std::endl;
	std::cout << " -m / --model             Admixture model file." << std::endl;
	std::cout << " -h / --help              Print this help." << std::endl;
	exit(0);
}

int main(int argc, char **argv)
{
	archaicSeekerPar ASpar;
	genomePar gpar;
	if(argc <= 1)
		archaicSeekerHMMHelp();
	char optstring[] = "a:T:e:o:v:p:r:X:A:m:h";
	char opt;
	std::map<char, int> parCheck;
	while((opt = getopt_long_only(argc, argv, optstring, archaicSeekerHMMOpt, NULL)) != -1)
	{
		switch (opt)
		{
			case 'a': ASpar.alpha = atof(optarg); break;
			case 'T': ASpar.introT = atof(optarg); break;
			case 'e': ASpar.emit = atof(optarg); break;
			case 'o': ASpar.outPrefix = optarg; break;
			case 'v': gpar.vcfPar = optarg; break;
			case 'p': gpar.popPar = optarg; break;
			case 'r': gpar.remapPar = optarg; break;
			case 'X': gpar.outgroupPar = optarg; break;
			case 'A': gpar.ancestralPar = optarg; break;
			case 'm': ASpar.seekModel = optarg; break;
			case 'h': archaicSeekerHMMHelp(); break;
			case '?': std::cout << "Cannot identify parameter " << opt << std::endl;
						archaicSeekerHMMHelp();
		}
		++parCheck[opt];
	}
	assert(parCheck['a'] <= 1);
	assert(parCheck['T'] <= 1);
	assert(parCheck['e'] <= 1);
	assert(parCheck['m'] == 1);
	assert(parCheck['o'] == 1);
	assert(parCheck['v'] == 1);
	assert(parCheck['p'] == 1);
	assert(parCheck['r'] == 1);
	assert(parCheck['X'] == 1);
	assert(parCheck['A'] == 1);

	std::ifstream fpi(ASpar.seekModel.c_str());
	getline(fpi, ASpar.seekModel);
//	omp_set_num_threads(ASpar.nthreads);
	gpar.load();
	archaicSeekerHMM(ASpar, gpar);
}










