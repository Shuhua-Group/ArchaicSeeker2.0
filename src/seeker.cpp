/*
 * seeker.cpp
 *
 *  Created on: Oct 10, 2018
 *      Author: yuankai
 */

# include "seeker.hpp"
# include "data.hpp"
# include "matching.hpp"

# include <unistd.h>
# include <cmath>
# include <fstream>
# include <cassert>
# include <cctype>

# include <nlopt.h>

class transMat
{
public:
	void setP(double a, double t, double d)
	{
		double e = exp( - t * d );
		if(d < 0)
			e = 0;
		/*
		 * pro[0] Modern  -> Modern
		 * pro[1] Archaic -> Modern
		 * pro[2] Modern  -> Archaic
		 * pro[3] Archaic -> Archaic
		 */
		pro[0] = 1 - a * ( 1 - e );
		pro[1] = ( 1 - a ) * ( 1 - e );
		pro[2] = a * ( 1 - e );
		pro[3] = e + a * ( 1 - e );
	}
	long double pro[4];
};

class stateHMM
{
public:
	void setInit(double a, double e0, double e1)
	{
		pro[0] = (1 - a) * e0;
		pro[1] = a * e1;
		c = 0;
		check();
	}
	void trans(const stateHMM& pre, const transMat& transP, double e0, double e1)
	{
		pro[0] = pre.pro[0] * transP.pro[0] + pre.pro[1] * transP.pro[1];
		pro[0] *= e0;
		pro[1] = pre.pro[0] * transP.pro[2] + pre.pro[1] * transP.pro[3];
		pro[1] *= e1;
		c = pre.c;
		check();
	}
	int bestState() const
	{
		if(pro[0] < pro[1])
			return 1;
		else
			return 0;
	}
	long double pro[2];
	long c;
	static long double scale ;
private:
	void check()
	{
		if(pro[0] < scale)
		{
			pro[0] /= scale;
			pro[1] /= scale;
			++c;
		}
	}
	//0 Modern, 1 Archaic
};

long double stateHMM::scale = 1.0e-150 ;

void stateForward(const std::vector<transMat>& transP, std::vector<stateHMM>& sites, \
		const std::vector<uint8_t>& state, double emit[2][4])
{
	int nsite = sites.size();
	double* e0 = emit[0];
	double* e1 = emit[1];
	for(int i = 1; i < nsite; ++i)
		sites[i].trans(sites.at(i - 1), transP.at(i), e0[state[i]], e1[state[i]]);
}

void stateViterbi(const std::vector<stateHMM>& sites, std::vector<uint8_t>& bestState)
{
	int nsite = sites.size();
	bestState.clear();
	bestState.resize(nsite, 0);
	for(int i = 0; i < nsite ; ++i)
		bestState[i] = sites[i].bestState();
}

class bed
{
public:
	bed() : pstart(0), pend(0), start(0), end(0), t(0)
	{
		diff = new double [nleaf];
		memset(diff, 0, sizeof(double) * nleaf);
		memset(state, 0, sizeof(int) * 4);
	}
	~bed() { delete []diff; }
	bed(const bed& rhs)
	{
		if(this != &rhs)
		{
			pstart = rhs.pstart;
			pend = rhs.pend;
			start = rhs.start;
			end = rhs.end;
			t = rhs.t;
			memcpy(state, rhs.state, sizeof(int) * 4);
			diff = new double [nleaf];
			memcpy(diff, rhs.diff, sizeof(double) * nleaf);
		}
	}
	void clear()
	{
		memset(diff, 0, sizeof(double) * nleaf);
		memset(state, 0, sizeof(int) * 4);
	}
	int state[4];
	long pstart, pend;
	double start, end, t;
	double *diff;

	static int nleaf;
};

int bed::nleaf;

void archaicSeekerHMM(archaicSeekerPar& ASpar, const genomePar& gpar)
{
	genome data(gpar);
	int ncontig(data.comContigs.size());
	int nind(data.comIDs.size());

	matchModel seekM(ASpar.seekModel);
	int nleaf(seekM.leafLabs.size());
	int nlab(seekM.labels.size());
	bed::nleaf = nleaf;
	std::vector<std::string>& matLabels(seekM.labels);
	std::map<std::string, int> popCorrespond;
	std::vector<int> archaicPops;
	for(int i = 0 ; i < nleaf ; ++i)
		popCorrespond[seekM.leafLabs[i]] = i;
	for(int i = 0 ; i < nlab; ++i)
	{
		if(gpar.ASpopCheck.count(matLabels[i]) && gpar.ASpopCheck.at(matLabels[i]) == 0)
			archaicPops.push_back(i);
	}
	int narchaicPops(archaicPops.size());

	std::vector<int> indexPop;
	indexPop.resize(nind);
	for(int i = 0 ; i < nind; ++i)
	{
		assert(popCorrespond.count(data.comPop[i]));
		indexPop[i] = popCorrespond[data.comPop[i]];
	}

	std::vector<std::vector<std::vector<double> > > altAf, refAf;
	altAf.resize(ncontig);
	refAf.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		altAf[i].resize(nleaf);
		refAf[i].resize(nleaf);
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nleaf; ++j)
		{
			altAf[i][j].resize(nsite, 0);
			refAf[i][j].resize(nsite, 0);
		}
	}

	std::vector<std::vector<char> > outGroup;
	outGroup.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::vector<char>& curOutGroup(outGroup[i]);
		const std::vector<long>& curPos(data.comPos[i]);
		const std::vector<char>& curRefAllele(data.comRefAllele[i]), \
				& curAltAllele(data.comAltAllele[i]);
		int nsite(curPos.size());
		curOutGroup.resize(nsite);
		assert(-1 != access(gpar.outgroupFastaPaths.at(data.comContigs[i]).c_str(), R_OK));
		gzifstream fpog(gpar.outgroupFastaPaths.at(data.comContigs[i]));
		std::string faInfo;
		getline(fpog, faInfo);
		int p(0);
		int total(0);
		while(p < nsite)
		{
			while(total < curPos[p])
			{
				getline(fpog, faInfo);
				if(!fpog)
					break;
				total += faInfo.size();
			}
			if(!fpog)
			{
				do
				{
					curOutGroup[p] = -1;
					++p;
				}
				while(p < nsite);
				break;
			}
			char allele(toupper(faInfo[curPos[p] - total + faInfo.size() - 1]));
			if(allele == curRefAllele[p])
				curOutGroup[p] = 0;
			else if(allele == curAltAllele[p])
				curOutGroup[p] = 1;
			else if(allele == 'A' || allele == 'T' || allele == 'G' || allele == 'C')
				curOutGroup[p] = 2;
			else
				curOutGroup[p] = -1;
			++p;
		}

		std::vector<std::vector<int> > altCount, sumCount;
		altCount.resize(nleaf);
		sumCount.resize(nleaf);
		for(int j = 0 ; j < nleaf; ++j)
		{
			altCount[j].resize(nsite, 0);
			sumCount[j].resize(nsite, 0);
		}
		for(int j = 0 ; j < nind;  ++j)
		{
			const hapInfo& curHap1(data.haplotypes[i][j * 2]), \
					& curHap2(data.haplotypes[i][j * 2 + 1]);
			std::vector<int>& curAlt(altCount[indexPop[j]]), \
					& curSum(sumCount[indexPop[j]]);
			if(data.phased[j])
			{
				for(int n = 0 ; n < 2; ++n)
				{
					const hapInfo& curH(data.haplotypes[i][j * 2 + n]);
					for(int k = 0 ; k < nsite; ++k)
					{
						if(curH[k])
							++curAlt[k];
						++curSum[k];
					}
				}
			}
			else
			{
				for(int k = 0 ; k < nsite; ++k)
				{
					if(curHap1[k] && !curHap2[k])
					{
						++curAlt[k];
						curSum[k] += 2;
					}
					else if(!curHap1[k] && curHap2[k])
					{
						curAlt[k] += 2;
						curSum[k] += 2;
					}
					else if(!(curHap1[k] || curHap2[k]))
						curSum[k] += 2;
				}
			}
		}
		for(int j = 0; j < nleaf; ++j)
		{
			std::vector<int>& curAlt(altCount[j]), & curSum(sumCount[j]);
			std::vector<double>& curAltAf(altAf[i][j]), & curRefAf(refAf[i][j]);
			for(int k = 0 ; k < nsite; ++k)
			{
				if(curSum[k])
				{
					curAltAf[k] = (double)curAlt[k] / curSum[k];
					curRefAf[k] = 1 - curAltAf[k];
				}
			}
		}
	}
	std::map<int, int> outgroupCheck;
	for(int i = 0 ; i < ncontig; ++i)
	{
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nsite; ++j)
			++outgroupCheck[outGroup[i][j]];
	}
//	std::cout << "Outgroup:" << std::endl;
//	for(int i = -1; i <=2; ++i)
//		std::cout << '\t' << i << '\t' << outgroupCheck[i] << std::endl;

	std::vector<std::vector<double> > pairDiff;
	pairDiff.resize(nleaf);
	for(int i = 0 ; i < nleaf; ++i)
		pairDiff[i].resize(nleaf, 0);
	for(int i = 0 ; i < ncontig; ++i)
	{
		const std::vector<std::vector<double> >& curRefAf(refAf[i]), & curAltAf(altAf[i]);
		const std::vector<char>& curOutGroup(outGroup[i]);
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nsite; ++j)
		{
			if(curOutGroup[j] >= 0)
			{
				for(int m = 0; m < nleaf; ++m)
					for(int n = 0; n < m ; ++n)
						pairDiff[m][n] += curRefAf[m][j] * curAltAf[n][j] + \
						curAltAf[m][j] * curRefAf[n][j];
				if(curOutGroup[j] == 0)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += curAltAf[m][j];
				else if(curOutGroup[j] == 1)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += curRefAf[m][j];
				else if(curOutGroup[j] == 2)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += 1;
			}
		}
	}
	std::string path;
//	path = ASpar.outPrefix + "_contigs.txt";
//	std::ofstream fpcontig(path.c_str());
//	fpcontig << "Contig\tTotalLength" << std::endl;
	const std::vector<std::string>& contigs = data.comContigs;
	long sumRegion(0);
	std::map<std::string, long> contigLen;
	std::map<std::string, long> hmmProp;
	for(int i = 0 ; i < ncontig; ++i)
	{
		long curContigLen(0);
		const std::vector<long>& curPos(data.comPos[i]);
		int curNsite(curPos.size());
		for(int j = 1; j < curNsite; ++j)
		{
			long len(curPos[j] - curPos[j - 1]);
			if(len < MAXINTER)
				curContigLen += len;
		}
		contigLen[contigs[i]] = curContigLen;
		sumRegion += curContigLen;
//		fpcontig << contigs[i] << '\t' << curContigLen << std::endl;
	}

	seekM.modelCalibrate(pairDiff, sumRegion);
//	path = ASpar.outPrefix + "_calibratedModel.txt";
//	std::ofstream fpcm(path.c_str());
//	for(int i = 0 ; i < nleaf; ++i)
//		for(int j = 0 ; j < i ; ++j)
//			fpcm << seekM.leafLabs[i] << '\t' << seekM.leafLabs[j] << '\t' << \
//					pairDiff[i][j] << std::endl;
//	for(int i = 0 ; i < nlab ; ++i)
//		fpcm << seekM.labels[i] << '\t' << seekM.calibratedLen[i] << std::endl;

	std::vector<int>& asPopLabel = data.comASPop;


	double emit[2][4];
	/*
	 * emit[0][] Modern
	 * emit[1][] Archaic
	 * emit[][0-3] stands for observed state 2, 1, 3, 4
	 */

	emit[0][1] = ASpar.emit * ( 1 - ASpar.emit ); //ASpar.emit * ( 1 - ASpar.emit )
	emit[0][0] = ( 1 - ASpar.emit ) * ( 1 - ASpar.emit );
	emit[0][2] = ASpar.emit * ASpar.emit;
	emit[0][3] = ASpar.emit * ( 1 - ASpar.emit );//emmit[0][1];
	emit[1][1] = emit[0][3];
	emit[1][0] = emit[0][2];
	emit[1][2] = emit[0][0];
	emit[1][3] = emit[0][1];

	std::vector<std::vector<std::vector<uint8_t> > > observS;
	observS.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		observS[i].resize(nind * 2);
//		for(int j = 0 ; j < nind * 2 ; ++j)
//			observS[i][j].resize(4);
	}
	for(int i = 0 ; i < ncontig; ++i)
	{
		hapInfo afrT, afrF, archT, archF, stateCheck;
		std::vector<long>& curComPos(data.comPos[i]);
		int nsite(curComPos.size());
		afrT.resize(nsite);
		afrT.set();
		afrF.resize(nsite);
		afrF.reset();
		archT.resize(nsite);
		archT.set();
		archF.resize(nsite);
		archF.reset();
		std::vector<hapInfo>& curHap = data.haplotypes[i];;
		for(int j = 0 ; j < nind; ++j)
		{
			if(asPopLabel[j] == 0)
			{
				stateCheck = curHap[j * 2] & curHap[j * 2 + 1];
				stateCheck.flip();
				archT &= curHap[j * 2 + 1] & stateCheck;
				archF |= ( curHap[j * 2] | curHap[j * 2 + 1] ) & stateCheck;
			}
			else if(asPopLabel[j] == 1)
			{
				stateCheck = curHap[j * 2] & curHap[j * 2 + 1];
				stateCheck.flip();
				afrT &= curHap[j * 2 + 1] & stateCheck;
				afrF |= ( curHap[j * 2] | curHap[j * 2 + 1] ) & stateCheck;
			}
		}
		archF.flip();
		afrF.flip();
		hapInfo afrAnc, archAnc;
		hapInfo& curAnc(data.ancestralState[i]);
		hapInfo curAncN(~curAnc);
		afrAnc = ( afrT & curAnc ) | ( afrF & curAncN );
		archAnc = ( archT & curAnc ) | ( archF & curAncN );
		hapInfo afrAncN, archAncN;
		afrAncN = ~afrAnc;
		archAncN = ~archAnc;

		for(int j = 0; j < nind * 2; ++j)
		{
			if(asPopLabel[j / 2] == 2)
			{
				std::vector<uint8_t>& obs(observS[i][j]);
				obs.resize(nsite);
				hapInfo s[4];
				hapInfo& cur(curHap[j]);
				hapInfo curDer;
				curDer = cur ^ curAnc;
				//s[0] state2, s[1] state1, s[2] state3, s[3] state4
				s[1] = curDer & afrAnc & archAnc;
				s[0] = curDer & archAncN & afrAnc;
				s[2] = curDer & afrAncN & archAnc;
				s[3] = ~(s[0] | s[1] | s[2]);
				for(int k = 0 ; k < nsite ; ++k)
					for(int l = 0 ; l < 4 ; ++l)
						if(s[l][k])
						{
							obs[k] = l;
							break;
						}
			}
		}
	}
/*
	for(int i = 0 ; i < ncontig; ++i)
	{
		hapInfo afrT, afrF, archT, archF, stateCheck;
		std::vector<long>& curComPos(data.comPos[i]);
		int nsite(curComPos.size());
		afrT.resize(nsite);
		afrT.set();
		afrF.resize(nsite);
		afrF.reset();
		archT.resize(nsite);
		archT.set();
		archF.resize(nsite);
		archF.reset();
		std::vector<hapInfo>& curHap = data.haplotypes[i];;
		for(int j = 0 ; j < nind; ++j)
		{
			if(asPopLabel[j] == 0)
			{
				stateCheck = curHap[j * 2] & curHap[j * 2 + 1];
				stateCheck.flip();
				archT &= curHap[j * 2 + 1] & stateCheck;
				archF |= ( curHap[j * 2] | curHap[j * 2 + 1] ) & stateCheck;
			}
			else if(asPopLabel[j] == 1)
			{
				stateCheck = curHap[j * 2] & curHap[j * 2 + 1];
				stateCheck.flip();
				afrT &= curHap[j * 2 + 1] & stateCheck;
				afrF |= ( curHap[j * 2] | curHap[j * 2 + 1] ) & stateCheck;
			}
		}
		archF.flip();
		afrF.flip();
		hapInfo afrTarchT, afrFarchF, afrTarchNT, afrFarchNF, \
				afrNTarchT, afrNFarchF, afrNTarchNT, afrNFarchNF;
		afrTarchT = afrT & archT;
		afrFarchF = afrF & archF;
		afrTarchNT = afrT & (~archT);
		afrFarchNF = afrF & (~archF);
		afrNTarchT = (~afrT) & archT;
		afrNFarchF = (~afrF) & archF;
		afrNTarchNT = (~afrT) & (~archT);
		afrNFarchNF = (~afrF) & (~archF);
		for(int j = 0; j < nind * 2; ++j)
		{
			if(asPopLabel[j / 2] == 2)
			{
				std::vector<uint8_t>& obs(observS[i][j]);
				obs.resize(nsite);
				hapInfo s[4];
				hapInfo& cur(curHap[j]);
				hapInfo curN(~cur);
				//s[0] state2, s[1] state1, s[2] state3, s[3] state4
				s[1] = ( afrTarchT & curN ) | ( afrFarchF & cur );
				s[0] = ( afrTarchNT & curN ) | ( afrFarchNF & cur );
				s[2] = ( afrNTarchT & curN ) | ( afrNFarchF & cur );
				s[3] = ( afrNTarchNT & curN ) | ( afrNFarchNF & cur );
				for(int k = 0 ; k < nsite ; ++k)
					for(int l = 0 ; l < 4 ; ++l)
						if(s[l][k])
						{
							obs[k] = l;
							break;
						}
			}
		}
	}
*/
	double sumGdis(0);
	for(int i = 0 ; i < ncontig; ++i)
		sumGdis += (data.gdisData[i].back() - data.gdisData[i].front()) / 100;
	int testInd(0);
	for(unsigned int i = 0 ; i < asPopLabel.size(); ++i)
		if(asPopLabel[i] == 2)
			++testInd;
	indexPop.resize(nind);
	for(int i = 0 ; i < nind; ++i)
	{
		assert(popCorrespond.count(data.comPop[i]));
		indexPop[i] = popCorrespond[data.comPop[i]];
	}

	std::set<long> introTcheck;
	for(int iter = 0 ; iter < 100 ; ++iter)
	{
		std::cout << "Start iter: " << iter << std::endl;
		std::cout << "\tAlpha: " << ASpar.alpha << "\n\tintroT: " << ASpar.introT << std::endl;
//		std::cout << "\tEmmit: " << std::endl;
//		for(int i = 0 ; i < 2 ; ++i)
//		{
//			for(int j = 0 ; j < 4 ; ++j)
//				std::cout << "\t\t[" << i << "][" << j << "] " << emit[i][j] << std::endl;
//		}

		double minMod(-log(0.975) / ASpar.alpha / ASpar.introT * 100);
		double minArch(-log(0.975) / (1 - ASpar.alpha) / ASpar.introT * 100);
		std::vector<std::vector<std::vector<bed> > > candidate;
		candidate.resize(ncontig);
		for(int i = 0 ; i < ncontig; ++i)
		{
			std::vector<long>& curComPos(data.comPos[i]);
			int nsite(curComPos.size());
			const std::vector<std::vector<double> >& curAltAf(altAf[i]);
			const std::vector<std::vector<double> >& curRefAf(refAf[i]);
			std::vector<transMat> transP;
			transP.resize(nsite);
			std::vector<double>& curGdis = data.gdisData[i];
			candidate[i].resize(nind * 2);
			for(int j = 1; j < nsite; ++j)
				transP[j].setP(ASpar.alpha, ASpar.introT, ( curGdis[j] - curGdis[j - 1] ) / 100);
			for(int j = 0; j < nind * 2; ++j)
			{
				if(asPopLabel[j / 2] == 2)
				{
					std::vector<bed>& curCandidate(candidate[i][j]);
					const std::vector<uint8_t>& s(observS[i][j]);
					const hapInfo& curTest(data.haplotypes[i][j]);
					std::vector<stateHMM> sites;
					sites.resize(nsite);
					sites.front().setInit(ASpar.alpha, emit[0][s[0]], emit[0][s[0]]);

					stateForward(transP, sites, s, emit);
					std::vector<uint8_t> states;
					stateViterbi(sites, states);
					bool preArchaic(false);
					bed tbed;
					tbed.pstart = curComPos[0];
					tbed.pend = curComPos[0];
					tbed.start = curGdis[0];
					tbed.end = curGdis[0];

					++tbed.state[s[0]];
					if(curTest[0])
					{
						for(int n = 0 ; n < nleaf; ++n)
							tbed.diff[n] += curRefAf[n][0];
					}
					else
					{
						for(int n = 0 ; n < nleaf ; ++n)
							tbed.diff[n] += curAltAf[n][0];
					}
					if(states[0])
					{
						tbed.t = 1;
						preArchaic = true;
					}
					else
					{
						tbed.t = -1;
						preArchaic = false;
					}
					for(int k = 1 ; k < nsite; ++k)
					{
						double gdis = curGdis[k];
						long pos = curComPos[k];
						if(states[k])
						{
							if(!preArchaic)
							{
								tbed.end = gdis;
								tbed.pend = pos;
								curCandidate.push_back(tbed);

								tbed.clear();
								tbed.t = 1;
								tbed.start = gdis;
								tbed.pstart = pos;
								preArchaic = true;
							}
							tbed.end = gdis;
							tbed.pend = pos;
							++tbed.state[s[k]];
							if(curTest[k])
							{
								for(int n = 0 ; n < nleaf; ++n)
									tbed.diff[n] += curRefAf[n][k];
							}
							else
							{
								for(int n = 0 ; n < nleaf ; ++n)
									tbed.diff[n] += curAltAf[n][k];
							}
						}
						else
						{
							if(preArchaic)
							{
								curCandidate.push_back(tbed);
								tbed.clear();
								tbed.t = -1;
								tbed.start = curGdis[k - 1];
								tbed.pstart = curComPos[k - 1];
								preArchaic = false;
							}
							tbed.end = gdis;
							tbed.pend = pos;
							++tbed.state[s[k]];
							if(tbed.end - tbed.start >= minMod)
							{
								if(curTest[k])
								{
									for(int n = 0 ; n < nleaf; ++n)
										tbed.diff[n] += curRefAf[n][k];
								}
								else
								{
									for(int n = 0 ; n < nleaf ; ++n)
										tbed.diff[n] += curAltAf[n][k];
								}
							}
						}
					}
					if(preArchaic)
						curCandidate.push_back(tbed);
				}
			}
		}
		std::cout << "HMM done" << std::endl;
		double diff[nleaf], llk[nlab], mat[nlab];
		int obsAcc[4], obsUpd[4];
		std::vector<std::vector<double> > segs;
		segs.resize(nlab);
		std::vector<std::vector<int> > obsCount;
		obsCount.resize(nlab);
		for(int i = 0 ; i < nlab; ++i)
			obsCount[i].resize(4, 0);
		std::vector<int> modCount;
		modCount.resize(4, 0);
		for(int i = 0 ; i < ncontig; ++i)
			for(int j = 0 ; j < nind * 2; ++j)
			{
				std::vector<bed>& curCandidate(candidate[i][j]);
				int nseg(curCandidate.size());
				int testPop = indexPop[j / 2];
				for(int k = 0; k < nseg; ++k)
				{
					const bed& startBed(curCandidate[k]);
					double segStart, segEnd;
					long segPstart;
					if(startBed.t > 0)
					{
						int preMatch = seekM.segMatchLK (startBed.diff, \
								startBed.pend - startBed.pstart, testPop, mat, llk);
						++k;
						segStart = startBed.start;
						segPstart = startBed.pstart;
						segEnd = startBed.end;
						for(int l = 0 ; l < nleaf; ++l)
							diff[l] = startBed.diff[l];
						for(int l = 0 ; l < 4 ; ++l)
							obsAcc[l] = startBed.state[l];
						memset(obsUpd, 0, sizeof(int) * 4);
						while(k < nseg)
						{
							const bed& curBed(curCandidate[k]);
							for(int l = 0 ; l < nleaf ; ++l)
								diff[l] += curBed.diff[l];
							for(int l = 0 ; l < 4 ; ++l)
								obsUpd[l] += curBed.state[l];
							if(curBed.t > 0)
							{
								int match = seekM.segMatchLK(curBed.diff, \
										curBed.pend - curBed.pstart, testPop, mat, llk);
								if(match == preMatch)
								{
									match = seekM.segMatchLK(diff, \
											curBed.pend - segPstart, testPop, mat, llk);
									if(match == preMatch)
									{
										segEnd = curBed.end;
										for(int l = 0 ; l < 4 ; ++l)
											obsAcc[l] += obsUpd[l];
										memset(obsUpd, 0, sizeof(int) * 4);
									}
									else
									{
										segs[preMatch].push_back(segEnd - segStart);
										for(int l = 0 ; l < 4 ; ++l)
										{
											obsCount[preMatch][l] += obsAcc[l];
											modCount[l] += obsUpd[l] - curBed.state[l];
										}
										--k;
										break;
									}
								}
								else
								{
									segs[preMatch].push_back(segEnd - segStart);
									for(int l = 0 ; l < 4 ; ++l)
									{
										obsCount[preMatch][l] += obsAcc[l];
										modCount[l] += obsUpd[l] - curBed.state[l];
									}
									--k;
									break;
								}
							}
							else if(curBed.end - curBed.start > minMod)
							{
								segs[preMatch].push_back(segEnd - segStart);
								for(int l = 0 ; l < 4 ; ++l)
								{
									obsCount[preMatch][l] += obsAcc[l];
									modCount[l] += curBed.state[l];
								}
								break;
							}
							++k;
						}
						if(k == nseg)
						{
							segs[preMatch].push_back(segEnd - segStart);
							for(int l = 0 ; l < 4 ; ++l)
								obsCount[preMatch][l] += obsAcc[l];
						}
					}
				}
			}
		std::cout << "Matching done" << std::endl;
		double archaicLen(0);
		int narchaic(0);
		double temit[2][4];
		memset(temit[1], 0, sizeof(double) * 4);

//		for(int i = 0 ; i < nlab ; ++i)
//			std::cout << i << '\t' << matLabels[i] << '\t' <<  segs[i].size() << std::endl;
		for(int i = 0 ; i < narchaicPops; ++i)
		{
			int nseg(segs[archaicPops[i]].size());
			for(int j = 0 ; j < nseg ; ++j)
			{
				double s(segs[archaicPops[i]][j]);
				if(s > minArch)
				{
					archaicLen += s;
					++narchaic;
				}
			}
			for(int j = 0 ; j < 4 ; ++j)
				temit[1][j] += obsCount[archaicPops[i]][j];
		}
		double sumemit(0);
		for(int i = 0 ; i < 3 ; ++i)
		{
			temit[0][i] = modCount[i];
			sumemit += temit[0][i];
		}
		temit[0][3] = sumemit * (1 - ASpar.emit);
		sumemit += temit[0][3];
		for(int i = 0 ; i < 4 ; ++i)
			emit[0][i] = temit[0][i] / sumemit;
		if(emit[0][0] < 1e-10)
			emit[0][0] = 1e-10;
		sumemit = 0;
		for(int i = 0 ; i < 2 ; ++i)
			sumemit += temit[1][i];
		temit[1][2] = sumemit * (1 - ASpar.emit) * (1 - ASpar.emit);
		temit[1][3] = sumemit * (1 - ASpar.emit);
		sumemit += temit[1][2];
		sumemit += temit[1][3];
		for(int i = 0 ; i < 4 ; ++i)
			emit[1][i] = temit[1][i] / sumemit;

		double averLen(archaicLen / narchaic / 100);
		double talpha, tintroT;
		talpha = archaicLen / 100 / sumGdis / testInd / 2;
		tintroT = 1.0 / (1 - talpha) / averLen;

		long tcheck((long)(tintroT * 100));

		if(introTcheck.count(tcheck) || abs(tintroT - ASpar.introT) < 0.01)
		{
			ASpar.introT = tintroT;
			ASpar.alpha = talpha;
			std::cout << "Alpha and introT converged after " << iter << " iteration." << std::endl;
			std::cout << "\tAlpha: \t" << ASpar.alpha << std::endl;
			std::cout << "\tintroT: \t" << ASpar.introT << std::endl;
			std::cout << "\tEmmit:" << std::endl;
			for(int i = 0 ; i < 2 ; ++i)
			{
				for(int j = 0 ; j < 4 ; ++j)
					std::cout << "\t\t[" << i << "][" << j << "] " << emit[i][j] << std::endl;
			}
			break;
		}
		else
		{
			ASpar.introT = tintroT;
			ASpar.alpha = talpha;
			introTcheck.insert(tcheck);
		}
	}
	{
		//output
		std::string path;
//		path = ASpar.outPrefix + "_hmm.seg";
//		std::ofstream fpohmmseg(path.c_str());
//		fpohmmseg << "ID\tContig\tStart(bp)\tEnd(bp)\tStart(cM)\tEnd(cM)" << std::endl;
//		path = ASpar.outPrefix + "_hmm.sum";
//		std::ofstream fpohmmsum(path.c_str());
//		fpohmmsum << "ID\tArchaicProp(cM)" << std::endl;
		double minMod(-log(0.975) / ASpar.alpha / ASpar.introT * 100);
		double minArch(-log(0.975) / (1 - ASpar.alpha) / ASpar.introT * 100);
		//double minMod(0), minArch(0);
		std::vector<std::vector<std::vector<bed> > > candidate;
		candidate.resize(ncontig);
		double hmmLen[nind * 2];
		memset(hmmLen, 0, sizeof(double) * nind * 2);
		for(int i = 0 ; i < ncontig; ++i)
		{
			std::vector<long>& curComPos(data.comPos[i]);
			int nsite(curComPos.size());
			const std::vector<std::vector<double> >& curAltAf(altAf[i]);
			const std::vector<std::vector<double> >& curRefAf(refAf[i]);
			std::vector<transMat> transP;
			transP.resize(nsite);
			std::vector<double>& curGdis = data.gdisData[i];
			candidate[i].resize(nind * 2);
			for(int j = 1; j < nsite; ++j)
				transP[j].setP(ASpar.alpha, ASpar.introT, ( curGdis[j] - curGdis[j - 1] ) / 100);
			for(int j = 0; j < nind * 2; ++j)
			{
				if(asPopLabel[j / 2] == 2)
				{
					std::ostringstream curHapID;
					curHapID << data.comIDs[j / 2] << "_" << (j % 2) + 1;
					const std::string& curHapIDStr(curHapID.str());
					std::vector<bed>& curCandidate(candidate[i][j]);
					const std::vector<uint8_t>& s(observS[i][j]);
					const hapInfo& curTest(data.haplotypes[i][j]);
					std::vector<stateHMM> sites;
					sites.resize(nsite);
					for(int k = 0; k < 4 ; ++k)
						sites.front().setInit(ASpar.alpha, emit[0][s[0]], emit[0][s[0]]);

					stateForward(transP, sites, s, emit);
					std::vector<uint8_t> states;
					stateViterbi(sites, states);
					bool preArchaic(false);
					bed tbed;
					tbed.pstart = curComPos[0];
					tbed.pend = curComPos[0];
					tbed.start = curGdis[0];
					tbed.end = curGdis[0];
					if(curTest[0])
					{
						for(int n = 0 ; n < nleaf; ++n)
							tbed.diff[n] += curRefAf[n][0];
					}
					else
					{
						for(int n = 0 ; n < nleaf ; ++n)
							tbed.diff[n] += curAltAf[n][0];
					}
					if(states[0])
					{
						tbed.t = 1;
						preArchaic = true;
					}
					else
					{
						tbed.t = -1;
						preArchaic = false;
					}
					for(int k = 1 ; k < nsite; ++k)
					{
						double gdis = curGdis[k];
						long pos = curComPos[k];
						if(states[k])
						{
							if(!preArchaic)
							{
								tbed.end = gdis;
								tbed.pend = pos;
								curCandidate.push_back(tbed);

								tbed.clear();
								tbed.t = 1;
								tbed.start = gdis;
								tbed.pstart = pos;
								preArchaic = true;
							}
							tbed.end = gdis;
							tbed.pend = pos;
							if(curTest[k])
							{
								for(int n = 0 ; n < nleaf; ++n)
									tbed.diff[n] += curRefAf[n][k];
							}
							else
							{
								for(int n = 0 ; n < nleaf ; ++n)
									tbed.diff[n] += curAltAf[n][k];
							}
						}
						else
						{
							if(preArchaic)
							{
								curCandidate.push_back(tbed);
								hmmLen[j] += tbed.end - tbed.start;
								//fpohmmseg << curHapIDStr << "\t" << contigs[i] << '\t' \
										<< tbed.pstart << '\t' << tbed.pend << '\t' \
										<< tbed.start << '\t' << tbed.end << std::endl;
								tbed.clear();
								tbed.t = -1;
								tbed.start = curGdis[k - 1];
								tbed.pstart = curComPos[k - 1];
								preArchaic = false;
							}
							tbed.end = gdis;
							tbed.pend = pos;
							if(tbed.end - tbed.start >= minMod)
							{
								if(curTest[k])
								{
									for(int n = 0 ; n < nleaf; ++n)
										tbed.diff[n] += curRefAf[n][k];
								}
								else
								{
									for(int n = 0 ; n < nleaf ; ++n)
										tbed.diff[n] += curAltAf[n][k];
								}
							}
						}
					}
					if(preArchaic)
					{
						curCandidate.push_back(tbed);
						hmmLen[j] += tbed.end - tbed.start;
						//fpohmmseg << curHapIDStr << "\t" << contigs[i] << '\t' \
								<< tbed.pstart << '\t' << tbed.pend << '\t' \
								<< tbed.start << '\t' << tbed.end << std::endl;
					}
				}
			}
		}
//		for(int j = 0 ; j < nind ; ++j)
//		{
//			if(asPopLabel[j] == 2)
//			{
//				fpohmmsum << data.comIDs[j] << '\t' << \
						(double) (hmmLen[j * 2] + hmmLen[j * 2 + 1]) / sumGdis / 2 \
						<< "%" << std::endl;
//			}
//		}

		path = ASpar.outPrefix + ".seg";
		std::ofstream fpomseg(path.c_str());
		fpomseg << "ID\tContig\tStart(bp)\tEnd(bp)\tStart(cM)\tEnd(cM)\tBestMatchedPop\tBestMatchedTime" << std::endl;
		path = ASpar.outPrefix + ".sum";
		std::ofstream fpomsum(path.c_str());
		fpomsum << "ID";
		for(int i = 0; i < nlab; ++i)
			fpomsum << '\t' << seekM.labels[i] << "(cM)";
		fpomsum << std::endl;
		double diff[nleaf], llk[nlab], mat[nlab];

		std::vector<std::vector<double> > mergeLen;
		mergeLen.resize(nind * 2);
		for(int i = 0 ; i < nind * 2; ++i)
			mergeLen[i].resize(nlab, 0);

		for(int i = 0 ; i < ncontig; ++i)
			for(int j = 0 ; j < nind * 2; ++j)
			{
				std::ostringstream curHapID;
				curHapID << data.comIDs[j / 2] << "_" << (j % 2) + 1;
				const std::string& curHapIDStr(curHapID.str());
				std::vector<bed>& curCandidate(candidate[i][j]);
				int nseg(curCandidate.size());
				int testPop = indexPop[j / 2];
				for(int k = 0; k < nseg; ++k)
				{
					const bed& startBed(curCandidate[k]);
					double segStart, segEnd;
					long segPstart, segPend;
					double bestT;
					if(startBed.t > 0)
					{
						int preMatch = seekM.segMatchLK (startBed.diff, \
								startBed.pend - startBed.pstart, testPop, mat, llk);
						++k;
						segStart = startBed.start;
						segPstart = startBed.pstart;
						segEnd = startBed.end;
						segPend = startBed.pend;
						if(preMatch >= 0)
							bestT = mat[preMatch];
						else
							bestT = 0;
						for(int l = 0 ; l < nleaf; ++l)
							diff[l] = startBed.diff[l];
						while(k < nseg)
						{
							const bed& curBed(curCandidate[k]);
							for(int l = 0 ; l < nleaf ; ++l)
								diff[l] += curBed.diff[l];
							if(curBed.t > 0)
							{
								int match = seekM.segMatchLK(curBed.diff, \
										curBed.pend - curBed.pstart, testPop, mat, llk);
								if(match == preMatch)
								{
									match = seekM.segMatchLK(diff, \
											curBed.pend - segPstart, testPop, mat, llk);
									if(match == preMatch)
									{
										segEnd = curBed.end;
										segPend = curBed.pend;
										if(preMatch >= 0)
											bestT = mat[preMatch];									}
									else
									{
										if(segEnd - segStart > minArch)
										{
											mergeLen[j][preMatch] += segEnd - segStart;
											fpomseg << curHapIDStr << '\t' << contigs[i] << '\t' \
												<< segPstart << '\t' << segPend << '\t' \
												<< segStart << '\t' << segEnd << '\t' \
												<< seekM.labels[preMatch] << '\t' \
												<< bestT << std::endl;
										}
										--k;
										break;
									}
								}
								else
								{
									if(segEnd - segStart > minArch)
									{
										mergeLen[j][preMatch] += segEnd - segStart;
										fpomseg << curHapIDStr << '\t' << contigs[i] << '\t' \
											<< segPstart << '\t' << segPend << '\t' \
											<< segStart << '\t' << segEnd << '\t' \
											<< seekM.labels[preMatch] << '\t' \
											<< bestT << std::endl;
									}
									--k;
									break;
								}
							}
							else if(curBed.end - curBed.start > minMod)
							{
								if(segEnd - segStart > minArch)
								{
									mergeLen[j][preMatch] += segEnd - segStart;
									fpomseg << curHapIDStr << '\t' << contigs[i] << '\t' \
										<< segPstart << '\t' << segPend << '\t' \
										<< segStart << '\t' << segEnd << '\t' \
										<< seekM.labels[preMatch] << '\t' \
										<< bestT << std::endl;
								}
								break;
							}
							++k;
						}
						if(k == nseg)
						{
							if(segEnd - segStart > minArch)
							{
								mergeLen[j][preMatch] += segEnd - segStart;
								fpomseg << curHapIDStr << '\t' << contigs[i] << '\t' \
									<< segPstart << '\t' << segPend << '\t' \
									<< segStart << '\t' << segEnd << '\t' \
									<< seekM.labels[preMatch] << '\t' \
									<< bestT << std::endl;
							}
						}
					}
				}
			}
		for(int i = 0 ; i < nind ; ++i)
		{
			if(asPopLabel[i] == 2)
			{
				fpomsum << data.comIDs[i];
				for(int j = 0 ; j < nlab; ++j)
				{
					fpomsum << '\t' << \
							(mergeLen[i * 2][j] + mergeLen[i * 2 + 1][j]) / sumGdis / 2 \
							<< '%';
				}
				fpomsum << std::endl;
			}
		}
		std::cout << "Matching done" << std::endl;
	}
}





/*
	std::vector<std::vector<std::vector<bed> > > region, regionPointer; //[region][haplotype]
	region.resize(ncontig);
	regionPointer.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		region[i].resize(nind * 2);
		regionPointer[i].resize(nind * 2);
		std::vector<long>& curComPos(data.comPos[i]);
		int nsite(curComPos.size());
//		std::vector<hapInfo>& curHap = data.haplotypes[i];
		std::vector<transMat> transP;
		transP.resize(nsite);
		std::vector<double>& curGdis = data.gdisData[i];
		for(int j = 1; j < nsite; ++j)
			transP[j].setP(ASpar.alpha, ASpar.introT, ( curGdis[j] - curGdis[j - 1] ) / 100);

		for(int j = 0; j < nind * 2; ++j)
		{
			if(asPopLabel[j / 2] == 2)
			{
				std::vector<bed>& curRegion(region[i][j]);
				std::vector<bed>& curRegionPointer(regionPointer[i][j]);
				const std::vector<hapInfo>& s(observS[i][j]);

				std::vector<stateHMM> sites;
				sites.resize(nsite);
				for(int k = 0; k < 4 ; ++k)
					if(s[k][0])
					//if(s[k].count(0, SEARCHMK))
					{
						sites.front().setInit(ASpar.alpha, emmit[0][k], emmit[1][k]);
						break;
					}
				stateForward(transP, sites, s, emmit);
				//std::vector<char> states;
				hapInfo states;
				stateViterbi(sites, states);
				long start(0), end(0), itStart(0), itEnd(0);
				bool preArchaic(false);
				for(int k = 0 ; k < nsite; ++k)
				{
					long pos = curComPos[k];
					if(!states[k])
					{
						if(preArchaic)
						{
							if(end != start)
							{
								bed t;
								t.start = start;
								t.end = end;
								curRegion.push_back(t);
								t.start = itStart;
								t.end = itEnd;
								curRegionPointer.push_back(t);
							}
						}
						preArchaic = false;
					}
					else if(states[k])
					{
						if(!preArchaic)
						{
							start = pos;
							itStart = k;
							preArchaic = true;
						}
						else if(preArchaic && pos - curComPos[k - 1] > MAXINTER)
						{
							bed t;
							t.start = start;
							t.end = end;
							curRegion.push_back(t);
							t.start = itStart;
							t.end = itEnd;
							curRegionPointer.push_back(t);
							start = pos;
							itStart = k;
						}
						end = pos;
						itEnd = k;
					}
				}
				if(preArchaic && end != start)
				{
					bed t;
					t.start = start;
					t.end = end;
					curRegion.push_back(t);
					t.start = itStart;
					t.end = itEnd;
					curRegionPointer.push_back(t);
				}
			}
		}
	}

	std::string path;
	path = ASpar.outPrefix + "_hmm.seg";
	std::ofstream fpseg(path.c_str());
	fpseg << "ID\tContig\tStart\tEnd" << std::endl;
	path = ASpar.outPrefix + "_hmm.sum";
	std::ofstream fpsum(path.c_str());
	fpsum << "ID\tArchaicProp" << std::endl;
	path = ASpar.outPrefix + "_contigs.txt";
	std::ofstream fpcontig(path.c_str());
	fpcontig << "Contig\tTotalLength" << std::endl;
	const std::vector<std::string>& contigs = data.comContigs;
	long sumRegion(0);
	std::map<std::string, long> contigLen;
	std::map<std::string, long> hmmProp;
	for(int i = 0 ; i < ncontig; ++i)
	{
		long curContigLen(0);
		const std::vector<long>& curPos(data.comPos[i]);
		int curNsite(curPos.size());
		for(int j = 1; j < curNsite; ++j)
		{
			long len(curPos[j] - curPos[j - 1]);
			if(len < MAXINTER)
				curContigLen += len;
		}
		contigLen[contigs[i]] = curContigLen;
		sumRegion += curContigLen;
		fpcontig << contigs[i] << '\t' << curContigLen << std::endl;
	}
	for(int i = 0 ; i < ncontig; ++i)
	{
		for(int j = 0 ; j < nind ; ++j)
		{
			if(asPopLabel[j] == 2)
			{
				std::string& curID(data.comIDs[j]);
				for(int k = 1; k <= 2 ; ++k)
				{
					long length(0);
					std::vector<bed>& curSeg = region[i][j * 2 + k - 1];
					int nseg = curSeg.size();
					std::ostringstream curHapID;
					curHapID << curID << "_" << k;
					std::string curHapIDStr(curHapID.str());
					if(nseg == 0)
					{
						hmmProp[curHapIDStr] = 0;
						continue;
					}
					for(int l = 0; l < nseg; ++l)
					{
						fpseg << curHapIDStr << '\t' << contigs[i] << '\t' << \
								curSeg[l].start << '\t' << curSeg[l].end << std::endl;
						length += curSeg[l].end - curSeg[l].start;
						hmmProp[curHapIDStr] += curSeg[l].end - curSeg[l].start;
					}
				}
			}
		}
	}
	for(std::map<std::string, long>::iterator it = hmmProp.begin(); it != hmmProp.end(); ++it)
		fpsum << it->first << '\t' << (double) it->second / sumRegion << std::endl;

	matchModel seekM(ASpar.seekModel);
	int nleaf(seekM.leafLabs.size());
	int nlab(seekM.labels.size());
	std::vector<std::string>& matLabels = seekM.labels;
	std::map<std::string, int> popCorrespond;
	for(int i = 0 ; i < nleaf ; ++i)
		popCorrespond[seekM.leafLabs[i]] = i;
	std::vector<int> indexPop;
	indexPop.resize(nind);
	for(int i = 0 ; i < nind; ++i)
	{
		assert(popCorrespond.count(data.comPop[i]));
		indexPop[i] = popCorrespond[data.comPop[i]];
	}

	std::vector<std::vector<std::vector<double> > > altAf, refAf;
	altAf.resize(ncontig);
	refAf.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		altAf[i].resize(nleaf);
		refAf[i].resize(nleaf);
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nleaf; ++j)
		{
			altAf[i][j].resize(nsite, 0);
			refAf[i][j].resize(nsite, 0);
		}
	}
	std::vector<std::vector<char> > outGroup;
	outGroup.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::vector<char>& curOutGroup(outGroup[i]);
		const std::vector<long>& curPos(data.comPos[i]);
		const std::vector<char>& curRefAllele(data.comRefAllele[i]), \
				& curAltAllele(data.comAltAllele[i]);
		int nsite(curPos.size());
		curOutGroup.resize(nsite);
		gzifstream fpog(gpar.outgroupFastaPaths.at(data.comContigs[i]));
		std::string faInfo;
		getline(fpog, faInfo);
		getline(fpog, faInfo);
		int block(0);
		for(int j = 0 ; j < nsite ; ++j)
		{
			long p(curPos[j]);
			int pBlock((p - 1) / 50);
			while(pBlock > block)
			{
				getline(fpog, faInfo);
				++block;
			}
			assert(pBlock == block);
			int pSite((p - 1) % 50);
			char allele(toupper(faInfo[pSite]));
			if(allele == curRefAllele[j])
				curOutGroup[j] = 0;
			else if(allele == curAltAllele[j])
				curOutGroup[j] = 1;
			else if(allele == 'A' || allele == 'T' || allele == 'G' || allele == 'C')
				curOutGroup[j] = 2;
			else
				curOutGroup[j] = -1;
		}
		std::vector<std::vector<int> > altCount, sumCount;
		altCount.resize(nleaf);
		sumCount.resize(nleaf);
		for(int j = 0 ; j < nleaf; ++j)
		{
			altCount[j].resize(nsite, 0);
			sumCount[j].resize(nsite, 0);
		}
		for(int j = 0 ; j < nind;  ++j)
		{
			const hapInfo& curHap1(data.haplotypes[i][j * 2]), \
					& curHap2(data.haplotypes[i][j * 2 + 1]);
			std::vector<int>& curAlt(altCount[indexPop[j]]), \
					& curSum(sumCount[indexPop[j]]);
			if(data.phased[j])
			{
				for(int n = 0 ; n < 2; ++n)
				{
					const std::vector<bed>& curP(regionPointer[i][j * 2 + n]);
					const hapInfo& curH(data.haplotypes[i][j * 2 + n]);
					int p(0), curNseg(curP.size());
					for(int k = 0 ; k < nsite; ++k)
					{
						if(p != curNseg)
						{
							if(curP[p].start == k)
							{
								k = curP[p].end;
								++p;
								continue;
							}
						}
						if(curH[k])
							++curAlt[k];
						++curSum[k];
					}
				}
			}
			else
			{
				for(int k = 0 ; k < nsite; ++k)
				{
					if(curHap1[k] && !curHap2[k])
					{
						++curAlt[k];
						curSum[k] += 2;
					}
					else if(!curHap1[k] && curHap2[k])
					{
						curAlt[k] += 2;
						curSum[k] += 2;
					}
					else if(!(curHap1[k] || curHap2[k]))
						curSum[k] += 2;
				}
			}
		}
		for(int j = 0; j < nleaf; ++j)
		{
			std::vector<int>& curAlt(altCount[j]), & curSum(sumCount[j]);
			std::vector<double>& curAltAf(altAf[i][j]), & curRefAf(refAf[i][j]);
			for(int k = 0 ; k < nsite; ++k)
			{
				if(curSum[k])
				{
					curAltAf[k] = (double)curAlt[k] / curSum[k];
					curRefAf[k] = 1 - curAltAf[k];
				}
			}
		}
	}
	std::map<int, int> outgroupCheck;
	for(int i = 0 ; i < ncontig; ++i)
	{
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nsite; ++j)
			++outgroupCheck[outGroup[i][j]];
	}
//	std::cout << "Outgroup:" << std::endl;
//	for(int i = -1; i <=2; ++i)
//		std::cout << '\t' << i << '\t' << outgroupCheck[i] << std::endl;
	std::vector<std::vector<double> > pairDiff;
	pairDiff.resize(nleaf);
	for(int i = 0 ; i < nleaf; ++i)
		pairDiff[i].resize(nleaf, 0);
	for(int i = 0 ; i < ncontig; ++i)
	{
		const std::vector<std::vector<double> >& curRefAf(refAf[i]), & curAltAf(altAf[i]);
		const std::vector<char>& curOutGroup(outGroup[i]);
		int nsite(data.comPos[i].size());
		for(int j = 0 ; j < nsite; ++j)
		{
			if(curOutGroup[j] >= 0)
			{
				for(int m = 0; m < nleaf; ++m)
					for(int n = 0; n < m ; ++n)
						pairDiff[m][n] += curRefAf[m][j] * curAltAf[n][j] + \
						curAltAf[m][j] * curRefAf[n][j];
				if(curOutGroup[j] == 0)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += curAltAf[m][j];
				else if(curOutGroup[j] == 1)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += curRefAf[m][j];
				else if(curOutGroup[j] == 2)
					for(int m = 0 ; m < nleaf; ++m)
						pairDiff[m][m] += 1;
			}
		}
	}
	seekM.modelCalibrate(pairDiff, sumRegion);
	path = ASpar.outPrefix + "_calibratedModel.txt";
	std::ofstream fpcm(path.c_str());
	for(int i = 0 ; i < nleaf; ++i)
		for(int j = 0 ; j < i ; ++j)
			fpcm << seekM.leafLabs[i] << '\t' << seekM.leafLabs[j] << '\t' << \
					pairDiff[i][j] << std::endl;
	for(int i = 0 ; i < nlab ; ++i)
		fpcm << seekM.labels[i] << '\t' << seekM.calibratedLen[i] << std::endl;

	double minMod(-log(0.975) / ASpar.alpha / ASpar.introT * 100);
	double minArch(-log(0.975) / (1 - ASpar.alpha) / ASpar.introT * 100);

	std::vector< std::map<std::string, std::vector<long> > > matProp;
	matProp.resize(ncontig);
	indexPop.resize(nind);
	for(int i = 0 ; i < nind; ++i)
	{
		assert(popCorrespond.count(data.comPop[i]));
		indexPop[i] = popCorrespond[data.comPop[i]];
	}
	std::vector<std::vector<std::vector<std::string> > > matchRe;
	matchRe.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
		matchRe[i].resize(nind * 2);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::map<std::string, std::vector<long> >& curMatProp(matProp[i]);
		const std::vector<double>& curGdis(data.gdisData[i]);
		const std::vector<long>& curPos(data.comPos[i]);
		const std::vector<std::vector<double> >& curAltAf(altAf[i]), & curRefAf(refAf[i]);
		for(int j = 0; j < nind * 2; ++j)
		{
			if(asPopLabel[j / 2] == 2)
			{
				int testPop = indexPop[j / 2];
				const hapInfo& curTest(data.haplotypes[i][j]);
				std::vector<bed>& curRegionPointer(regionPointer[i][j]);
				std::ostringstream curHapID;
				curHapID << data.comIDs[j / 2] << "_" << (j % 2) + 1;
				std::string curHapIDStr(curHapID.str());
				if(curMatProp[curHapIDStr].size() == 0)
					curMatProp[curHapIDStr].resize(nlab, 0);
				int nseg(curRegionPointer.size());
				if(nseg == 0)
					continue;
				std::vector<int> hmmMatch;
				hmmMatch.resize(nseg, 0);
				for(int k = 0 ; k < nseg; ++k)
				{
					double diff[nleaf], mat[nlab], llk[nlab];
					memset(diff, 0, sizeof(double) * nleaf);
					int start(curRegionPointer[k].start), end(curRegionPointer[k].end);
					for(int l = start; l <= end; ++l)
					{
						if(curTest[l])
						{
							for(int n = 0 ; n < nleaf; ++n)
								diff[n] += curRefAf[n][l];
						}
						else
						{
							for(int n = 0 ; n < nleaf; ++n)
								diff[n] += curAltAf[n][l];
						}
					}
					double curLen(curPos[end] - curPos[start]);
					hmmMatch[k] = seekM.segMatchLK(diff, curLen, testPop, mat, llk);
					curRegionPointer[k].t = mat[hmmMatch[k]];
				}
				bed t(curRegionPointer[0]);
				int preMatch(hmmMatch[0]);
				std::vector<std::string>& curResults(matchRe[i][j]);
				for(int k = 1 ; k < nseg; ++k)
				{
					//0.128
					if(curGdis[curRegionPointer[k].start] - curGdis[t.end] > minMod || \
							(t.end - curRegionPointer[k].start <= 1 && \
							curPos[t.end] - curPos[curRegionPointer[k].start] > MAXINTER))
					{
						if(curGdis[t.end] - curGdis[t.start] > minArch)
						{
							std::ostringstream tre;
							tre << curHapIDStr << '\t' << contigs[i] << '\t' << \
									curPos[t.start] << '\t' << curPos[t.end] << '\t';
							if(preMatch == -1)
								tre << "Unknown\t-";
							else
							{
								tre << matLabels[preMatch] << '\t' << t.t;
								curMatProp[curHapIDStr][preMatch] += curPos[t.end] - \
										curPos[t.start];
							}
							curResults.push_back(tre.str());
						}
						t.start = curRegionPointer[k].start;
						t.end = curRegionPointer[k].end;
						t.t = curRegionPointer[k].t;
						preMatch = hmmMatch[k];
					}
					else
					{
						if(preMatch == hmmMatch[k])
						{
							double diff[nleaf], mat[nlab], llk[nlab];
							memset(diff, 0, sizeof(double) * nleaf);
							for(int l = t.start; l <= curRegionPointer[k].end; ++l)
							{
								if(curTest[l])
								{
									for(int n = 0 ; n < nleaf; ++n)
										diff[n] += curRefAf[n][l];
								}
								else
								{
									for(int n = 0 ; n < nleaf; ++n)
										diff[n] += curAltAf[n][l];
								}
							}
							double curLen(curPos[curRegionPointer[k].end] - curPos[t.start]);
							int best = seekM.segMatchLK(diff, curLen, testPop, mat, llk);
							if(best == preMatch)
							{
								t.end = curRegionPointer[k].end;
								t.t = mat[best];
							}
							else
							{
								if(curGdis[t.end] - curGdis[t.start] > minArch)
								{
									std::ostringstream tre;
									tre << curHapIDStr << '\t' << contigs[i] << '\t' << \
											curPos[t.start] << '\t' << curPos[t.end] << '\t';
									if(preMatch == -1)
										tre << "Unknown\t-";
									else
									{
										tre << matLabels[preMatch] << '\t' << t.t;
										curMatProp[curHapIDStr][preMatch] += curPos[t.end] - \
												curPos[t.start];
									}
									curResults.push_back(tre.str());
								}
								t.start = curRegionPointer[k].start;
								t.end = curRegionPointer[k].end;
								t.t = curRegionPointer[k].t;
								preMatch = hmmMatch[k];
							}
						}
						else
						{
							if(curGdis[t.end] - curGdis[t.start] > minArch)
							{
								std::ostringstream tre;
								tre << curHapIDStr << '\t' << contigs[i] << '\t' << \
										curPos[t.start] << '\t' << curPos[t.end] << '\t';
								if(preMatch == -1)
									tre << "Unknown\t-";
								else
								{
									tre << matLabels[preMatch] << '\t' << t.t;
									curMatProp[curHapIDStr][preMatch] += curPos[t.end] - \
											curPos[t.start];
								}
								curResults.push_back(tre.str());
							}
							t.start = curRegionPointer[k].start;
							t.end = curRegionPointer[k].end;
							t.t = curRegionPointer[k].t;
							preMatch = hmmMatch[k];
						}
					}
				}
				if(curGdis[t.end] - curGdis[t.start] > minArch)
				{
					std::ostringstream tre;
					tre << curHapIDStr << '\t' << contigs[i] << '\t' << curPos[t.start] << \
							'\t' << curPos[t.end] << '\t';
					if(preMatch == -1)
						tre << "Unknown\t-";
					else
					{
						tre << matLabels[preMatch] << '\t' << t.t;
						curMatProp[curHapIDStr][preMatch] += curPos[t.end] - curPos[t.start];
					}
					curResults.push_back(tre.str());
				}
			}
		}
	}
	path = ASpar.outPrefix + "_mergeLK.seg";
	std::ofstream fpmerge(path.c_str());
	fpmerge << "ID\tContig\tStart\tEnd\tBestMatchedPop\tBestMatchedTime" << std::endl;
	for(int i = 0 ; i < ncontig; ++i)
		for(int j = 0 ; j < nind * 2; ++j)
			for(unsigned int k = 0 ; k < matchRe[i][j].size(); ++k)
				fpmerge << matchRe[i][j][k] << std::endl;
	path = ASpar.outPrefix + "_mergeLK.sum";
	std::ofstream fpmergeSum(path.c_str());
	fpmergeSum << "ID";
	for(int i = 0; i < nlab; ++i)
		fpmergeSum << '\t' << seekM.labels[i];
	fpmergeSum << std::endl;
	for(std::map<std::string, std::vector<long> >::iterator it = matProp[0].begin(); \
			it != matProp[0].end(); ++it)
	{
		std::string id(it->first);
		fpmergeSum << id;
		double len[nlab];
		memset(len, 0, sizeof(double) * nlab);
		for(int i = 0 ; i < ncontig; ++i)
		{
			for(int j = 0 ; j < nlab; ++j)
				len[j] += matProp[i][id][j];
		}
		for(int i = 0 ; i < nlab ; ++i)
			fpmergeSum << '\t' << len[i] / sumRegion;
		fpmergeSum << std::endl;
	}
}


*/

