/*
 * MultiWaveInfer.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */
#include <map>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Utils.hpp"

#include "omp.h"

using namespace std;


int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "Need more arguments than provided, please use -h/--help to get help" << endl;
		exit(1);
	}
	string filename = "";
	double lower = 0;
	double alpha = 0.001;
	double epsilon = 0.000001;
	double minP = 0.0001;
	int maxIter = 10000;
	int nbootstrap(100);
	int nthread(1);
	string outPrefix("");
	for (int i = 1; i < argc; ++i)
	{
		string arg(argv[i]);
		if (arg == "-h" || arg == "--help")
		{
			help();
			exit(0);
		}
		else if (arg == "-i" || arg == "--input")
		{
			filename = string(argv[++i]);
		}
		else if (arg == "-l" || arg == "--lower")
		{
			lower = atof(argv[++i]);
		}
		else if (arg == "-m" || arg == "--maxIter")
		{
			maxIter = atoi(argv[++i]);
		}
		else if (arg == "-a" || arg == "--alpha")
		{
			alpha = atof(argv[++i]);
		}
		else if (arg == "-e" || arg == "--epsilon")
		{
			epsilon = atof(argv[++i]);
		}
		else if (arg == "-p" || arg == "--minProp")
		{
			minP = atof(argv[++i]);
		}
		else if (arg == "-b" || arg == "--bootstrap")
		{
			nbootstrap = atol(argv[++i]);
		}
		else if (arg == "-o" || arg == "--output")
		{
			outPrefix = argv[++i];
		}
		else if (arg == "-t" || arg == "--thread")
		{
			nthread = atol(argv[++i]);
		}
		else
		{
			cerr << "unrecognizable argument found, please check again!" << endl;
			exit(1);
		}
	}
	if (filename.size() == 0)
	{
		cerr << "Input file name required, please check help" << endl;
		exit(1);
	}
	omp_set_num_threads(nthread);
	srand(time(NULL));
	vector<double> observ;
	//Reading test data from file named as test.dat
	ifstream fin(filename.c_str());
	if (!fin.is_open())
	{
		cerr << "Can't open file " << filename << endl;
		exit(1);
	}
	vector<double> segs;
	map<string, double> sumLengths;
	vector<string> labels;
	labels.push_back("Archaic");
	labels.push_back("Modern");
	double totalLength = 0;
	string line;
	while (getline(fin, line))
	{
		istringstream ss(line);
		double start, end;
		string label;
		ss >> start >> end >> label;
		double len = end - start;

		if (label == "Archaic")
		{
			if(len > lower)
				segs.push_back(len - lower);
		}
		else if (label != "Modern")
		{
			cerr << "Cannot identify " << label << endl;
			return 1;
		}
		sumLengths[label] += len;
		totalLength += len;
	}
	fin.close();

    if(segs.size() <= 0)
    {
    	cerr << "No effective segment length input." << endl;
    	return 1;
    }

	map<string, double> mixtureProps; //S_k
	ParamExp optPar;
	double criticalValue = cv_chisq(2, alpha);
	for (int i = 0 ; i < 2 ; ++i)
		mixtureProps[labels[i]] = sumLengths[labels[i]] / totalLength;

	optPar = findOptPar(segs, maxIter, mixtureProps.at("Archaic"), criticalValue, epsilon, minP, false);
	solveTrueProp(optPar, lower);

	int totalNumOfWaves = optPar.getK();
	vector<int> popOrder;
	double tempSum = 0;
	double temp[totalNumOfWaves];
	for (int i = 0 ; i < totalNumOfWaves; ++i)
	{
		popOrder.push_back(i);
		temp[i] = optPar.getProp(i) / optPar.getLambda(i);
		tempSum += temp[i] / temp[0];
	}
	tempSum = mixtureProps.at("Archaic") / tempSum;

	double mInOrder[totalNumOfWaves];
	for (int i = 0; i < totalNumOfWaves; ++i)
		mInOrder[i] = tempSum * temp[i] / temp[0];

	double alphaInOrder[totalNumOfWaves];
	double multiplier = 1.0;
	alphaInOrder[0] = mInOrder[0];
	for (int i = 1; i < totalNumOfWaves; ++i)
	{
		multiplier *= (1 - alphaInOrder[i - 1]);
		alphaInOrder[i] = mInOrder[i] / multiplier;
	}

	double hInOrder[totalNumOfWaves];
	hInOrder[totalNumOfWaves - 1] = alphaInOrder[totalNumOfWaves - 1];
	for(int i = totalNumOfWaves - 2; i >= 0; --i)
		hInOrder[i] = hInOrder[i + 1] * ( 1 - alphaInOrder[i]) + alphaInOrder[i];
//	hInOrder[totalNumOfWaves - 1] = hInOrder[totalNumOfWaves - 2];


	double admixTime[totalNumOfWaves];
	for (int i = 0; i < totalNumOfWaves; ++i)
	{
		double tempSum = 0;
		for (int k = 0; k < i; ++k)
		{
			tempSum += (1 - hInOrder[k]) * admixTime[k];
		}
		double rate = optPar.getLambda(i) - tempSum;
		admixTime[i] = rate / (1 - hInOrder[i]);
	}

	for (int i = 1; i < totalNumOfWaves; ++i)
	{
		admixTime[i] += admixTime[i - 1];
	}
//	cout << '\t';
//	for(int i = 0; i < totalNumOfWaves; ++i)
//		cout << admixTime[i] << '(' << alphaInOrder[i] << "),";
//	cout << endl;


	vector<vector<double> > bootstrapAdmixTime, bootstrapAlpha;
	bootstrapAdmixTime.resize(nbootstrap);
	bootstrapAlpha.resize(nbootstrap);
	string path;
	path = outPrefix + ".bootstrap";
	ofstream fpoboot(path.c_str());
	fpoboot << "NumBootstrap\tNumArchaic\tTime" << endl;

#pragma omp parallel for
	for(int iterBoot = 0 ; iterBoot < nbootstrap; ++iterBoot)
	{
		vector<double> curSegs;
		int nseg(segs.size());
		curSegs.resize(nseg);
		for(int i = 0 ; i < nseg ; ++i)
			curSegs[i] = segs[rand() % nseg];
		ParamExp curOptPar;

		curOptPar = findOptPar(curSegs, maxIter, mixtureProps.at("Archaic"), criticalValue, epsilon, minP, false);
		solveTrueProp(curOptPar, lower);

		int curTotalNumOfWaves = curOptPar.getK();
		vector<int> curPopOrder;
		double curTempSum = 0;
		double curTemp[curTotalNumOfWaves];
		for (int i = 0 ; i < curTotalNumOfWaves; ++i)
		{
			curPopOrder.push_back(i);
			curTemp[i] = curOptPar.getProp(i) / curOptPar.getLambda(i);
			curTempSum += curTemp[i] / curTemp[0];
		}
		curTempSum = mixtureProps.at("Archaic") / curTempSum;

		double mInOrder[curTotalNumOfWaves];
		for (int i = 0; i < curTotalNumOfWaves; ++i)
			mInOrder[i] = curTempSum * curTemp[i] / curTemp[0];

		vector<double> & alphaInOrder = bootstrapAlpha[iterBoot];
#pragma omp critical
		alphaInOrder.resize(curTotalNumOfWaves, 0);
		double multiplier = 1.0;
		alphaInOrder[0] = mInOrder[0];
		for (int i = 1; i < curTotalNumOfWaves; ++i)
		{
			multiplier *= (1 - alphaInOrder[i - 1]);
			alphaInOrder[i] = mInOrder[i] / multiplier;
		}

		double hInOrder[curTotalNumOfWaves];
		hInOrder[curTotalNumOfWaves - 1] = alphaInOrder[curTotalNumOfWaves - 1];
		for(int i = curTotalNumOfWaves - 2; i >= 0; --i)
			hInOrder[i] = hInOrder[i + 1] * ( 1 - alphaInOrder[i]) + alphaInOrder[i];
	//	hInOrder[totalNumOfWaves - 1] = hInOrder[totalNumOfWaves - 2];


		vector<double> & admixTime = bootstrapAdmixTime[iterBoot];
#pragma omp critical
		admixTime.resize(curTotalNumOfWaves, 0);
		for (int i = 0; i < curTotalNumOfWaves; ++i)
		{
			double tempSum = 0;
			for (int k = 0; k < i; ++k)
			{
				tempSum += (1 - hInOrder[k]) * admixTime[k];
			}
			double rate = curOptPar.getLambda(i) - tempSum;
			admixTime[i] = rate / (1 - hInOrder[i]);
		}

		for (int i = 1; i < curTotalNumOfWaves; ++i)
		{
			admixTime[i] += admixTime[i - 1];
		}
#pragma omp critical
		{
			fpoboot << iterBoot << "\t" << curTotalNumOfWaves << "\t";
			for(int i = 0 ; i < curTotalNumOfWaves - 1; ++i)
				fpoboot << admixTime[i] << "(" << alphaInOrder[i] << "),";
			fpoboot << admixTime[curTotalNumOfWaves - 1] << "(" << alphaInOrder[curTotalNumOfWaves - 1] << ")" << endl;
		}
	}
	map<int, int> wavesCount;
	map<int, vector<vector<double> > > adTimeSum, adProSum;
	for(int i = 0 ; i < nbootstrap; ++i)
	{
		int nwave = bootstrapAdmixTime[i].size();
		if(adTimeSum.count(nwave) == 0)
		{
			adTimeSum[nwave].resize(nwave);
			adProSum[nwave].resize(nwave);
		}
		++wavesCount[nwave];
		for(int j = 0 ; j < nwave ; ++j)
		{
			adTimeSum[nwave][j].push_back(bootstrapAdmixTime[i][j]);
			adProSum[nwave][j].push_back(bootstrapAlpha[i][j]);
		}
	}
	path = outPrefix + ".sum";
	ofstream fposum(path.c_str());
	fposum << "Dataset\tSupportRatio\tSupportNum\tNumOfArchaic\tTime(Generation)" << endl;
	fposum << "CompleteData\t-\t-\t" << totalNumOfWaves << "\t";
	for(int i = 0 ; i < totalNumOfWaves - 1; ++i)
		fposum << admixTime[i] << "(" << alphaInOrder[i] << "),";
	fposum << admixTime[totalNumOfWaves - 1] << "(" << alphaInOrder[totalNumOfWaves - 1] << ")" << endl;
	for(map<int, int>::iterator it = wavesCount.begin(); it != wavesCount.end(); ++it)
	{
		int n(it->second);
		int nwave(it->first);
		fposum << "BootstrapData\t" << (double) n / nbootstrap * 100 << "%\t" << n << '\t' << nwave << "\t";
		for(int i = 0 ; i < nwave - 1; ++i)
		{
			vector<double>& curAdTime(adTimeSum.at(nwave).at(i));
			sort(curAdTime.begin(), curAdTime.end());
			vector<double>& curPro(adProSum.at(nwave).at(i));
			sort(curPro.begin(), curPro.end());
			fposum << curAdTime[n * 0.025] << "~" << curAdTime[n * 0.975] << \
					"(" << curPro[n * 0.025] << "~" << curPro[n * 0.975] << "),";
		}
		vector<double>& curAdTime(adTimeSum.at(nwave).at(nwave - 1));
		sort(curAdTime.begin(), curAdTime.end());
		vector<double>& curPro(adProSum.at(nwave).at(nwave - 1));
		sort(curPro.begin(), curPro.end());
		fposum << curAdTime[n * 0.025] << "~" << curAdTime[n * 0.975] << \
				"(" << curPro[n * 0.025] << "~" << curPro[n * 0.975] << ")" << endl;
	}

	return 0;
}
