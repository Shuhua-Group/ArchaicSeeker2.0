/*
 * getArchaicMK.cpp
 *
 *  Created on: Mar 9, 2020
 *      Author: kyuan
 */

# include <iostream>
# include <fstream>
# include <sstream>
# include <map>
# include <set>
# include <vector>
# include <cstdlib>
# include <zlib.h>
# include "boost/dynamic_bitset.hpp"

using namespace std;

class reg
{
public:
	int chr;
	long start, end;
	boost::dynamic_bitset<> snp;
};

class pos
{
public:
	int chr;
	long p;
	string ref, alt;
	boost::dynamic_bitset<> snp;
};

//-1 pos < reg; 0 pos -> reg; 1 pos > reg
int cmp(const pos& p, const reg& r)
{
	if(p.chr < r.chr)
		return -1;
	if(p.chr > r.chr)
		return 1;
	if(p.p < r.start)
		return -1;
	if(p.p > r.end)
		return 1;
	return 0;
}

int cmp(const pos& p1, const pos& p2)
{
	if(p1.chr < p2.chr)
		return -1;
	if(p1.chr > p2.chr)
		return 1;
	if(p1.p < p2.p)
		return -1;
	if(p1.p > p2.p)
		return 1;
	return 0;
}

const int MAX = 4096 * 16;

int main(int argc, char **argv)
{
	gzFile fptest, fpas;
	fptest = gzopen(argv[1], "rb");
	fpas = gzopen(argv[2], "rb");
	ifstream fpind(argv[3]);
	string pop(argv[4]);
	ofstream fpo(argv[5]);

	set<string> ids;
	string id, p;
	while(fpind >> id >> p)
		if(p == pop)
			ids.insert(id);

	map<string, int> check;
	vector<int> sel;
	char buff[MAX];
	gzgets(fpas, buff, MAX);
	istringstream ashead(buff);
	string line;
	for(int i = 0 ; i < 9 ; ++i)
		ashead >> line;
	int n(0);
	while(ashead >> line)
		if(ids.count(line))
		{
			check[line] = n;
			++n;
			sel.push_back(1);
		}
		else
			sel.push_back(0);
	vector<reg> asReg;
	reg treg;
	treg.snp.resize(n * 2);
	while(gzgets(fpas, buff, MAX))
	{
		istringstream cl(buff);
		treg.snp.reset();
		cl >> treg.chr >> treg.start >> line >> line >> line >> line >> line >> line;
		treg.end = atol(line.substr(4).c_str());
		cl >> line;
		int m(0);
		for(int i = 0 ; i < sel.size() ; ++i)
		{
			cl >> line;
			if(sel[i])
			{
				if(line[0] == '1')
					treg.snp.set(m * 2);
				if(line[2] == '1')
					treg.snp.set(m * 2 + 1);
				++m;
			}
		}
		asReg.push_back(treg);
	}

	while(gzgets(fptest, buff, MAX))
		if(buff[1] != '#')
			break;
	istringstream head(buff);
	for(int i = 0 ; i < 9 ; ++i)
		head >> line;
	vector<int> index;
	int count(0);
	while(head >> id)
	{
		if(check.count(id))
			index.push_back(check[id]);
		else
			index.push_back(-1);
	}
	pos tpos;
	tpos.snp.resize(n * 2);
	int s(0);
	fpo << "CHR\tPos\tRef\tAlt\tAltArchaic\tNArchaic\tAltMod\tNMod\tSum" << endl;
	while(gzgets(fptest, buff, MAX))
	{
		istringstream cl(buff);
		cl >> tpos.chr >> tpos.p >> line >> tpos.ref >> tpos.alt >> line >> line >> line >> line;
		tpos.snp.reset();
		for(int i = 0 ; i < index.size(); ++i)
		{
			cl >> line;
			if(index[i] >= 0)
			{
				if(line[0] == '1')
					tpos.snp.set(index[i] * 2);
				if(line[2] == '1')
					tpos.snp.set(index[i] * 2 + 1);
			}
		}
		int c(cmp(tpos, asReg[s]));
		if(c == -1)
			continue;
		while(c == 1)
		{
			++s;
			if(s == asReg.size())
				break;
			if(cmp(tpos, asReg[s]) <= 0)
				break;
		}
		if(s == asReg.size())
			break;
		if(c == 0)
		{
				fpo << tpos.chr<<'\t' << tpos.p << '\t' << tpos.ref << '\t' << tpos.alt << '\t' << (tpos.snp & asReg[s].snp).count() << '\t' << asReg[s].snp.count() << '\t' << (tpos.snp & (~asReg[s].snp) ).count() << '\t' << (~asReg[s].snp).count() << '\t' << n * 2 << endl;
		}
	}
}


