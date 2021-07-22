/*
 * data.hpp
 *
 *  Created on: Sep 29, 2018
 *      Author: yuankai
 */

#ifndef DATA_HPP_
#define DATA_HPP_

# include <iostream>
# include <cstring>
# include <bitset>
# include <vector>
# include <map>
# include <sstream>

class hapInfo
{
public:
	hapInfo() : nBlock(0), nByte(0), capBlock(16), nInfo(0), mem(NULL)
			{ mem = new unsigned char[16]; 	memset(mem, 0, 16);}
	hapInfo(const hapInfo& rhs);
	hapInfo(int cap);
	~hapInfo() { delete[] mem; }

	void set(int n);
	void set();
	void reset(int n);
	void reset();
	void flip();
	int count() const;
	int count(int start, int end) const;
	int size() const;
	void resize(int s);
	void reserve(int s);

	hapInfo operator&(const hapInfo& rhs) const;
	void operator&=(const hapInfo& rhs);
	hapInfo operator|(const hapInfo& rhs) const;
	void operator|=(const hapInfo& rhs);
	hapInfo operator^(const hapInfo& rhs) const;
	void operator^=(const hapInfo& rhs);
	bool operator[](const int n) const;
	hapInfo& operator=(const hapInfo& rhs);
	hapInfo operator~() const;

	friend int hapDiff(const hapInfo& a, const hapInfo& b);
	friend int hapDiff(const hapInfo& a, const hapInfo& b, int start, int end);

#ifdef DEBUG
	void print()
	{
		std::cout << "nBlock:\t" << nBlock << std::endl;
		std::cout << "nByte:\t" << nByte << std::endl;
		std::cout << "capBlock:\t" << capBlock << std::endl;
		std::cout << "nInfo:\t" << nInfo << std::endl;
		for(int i = 0 ; i < capBlock; ++i)
		{
			std::cout << i << '\t';
		    std::cout << std::bitset<8>(mem[i]);
		    std::cout << std::endl;
		}
	}
	void simplePrint()
	{
		for(int i = 0 ; i < capBlock; ++i)
		    std::cout << std::bitset<8>(mem[i]) << '\t';
		std::cout << std::endl;
	}
	std::string to_string()
	{
		std::ostringstream s;
		for(int i = 0 ; i < nBlock; ++i)
		    s << std::bitset<8>(mem[i]) << '\t';
		if(nByte)
			s << std::bitset<8>(mem[nBlock]).to_string().substr(8-nByte);
		return s.str();
	}
#endif //DEBUG
private:
	//nBlock, nByte and nInfo are the current end of the haplotype; capBlock is the capability.
	int nBlock, nByte, capBlock, nInfo;

	unsigned char *mem;
};

int hapDiff(const hapInfo& a, const hapInfo& b);
int hapDiff(const hapInfo& a, const hapInfo& b, int start, int end);

class genomePar
{
public:
	genomePar(){};
	genomePar(const std::string& _vcfPar, const std::string& _popPar, \
			const std::string& _remapPar) : vcfPar(_vcfPar), popPar(_popPar), \
			remapPar(_remapPar) { vcfParLoad(); remapParLoad(); popParLoad(); }

	void load() {vcfParLoad(); remapParLoad(); popParLoad(); outgroupFastaPathLoad(); \
			ancestralFastaPathLoad(); }

	//Path to vcf par file, path to recombination par file and path to pop info file
	std::string vcfPar, popPar, remapPar;

	//Path to outgroup par file
	std::string outgroupPar;

	//Path to ancestral par file
	std::string ancestralPar;

	//vcf file paths, recombination map file paths, recombination map contigs
	std::vector<std::string> vcfPaths, remapPaths, remapContigs;

	//ID info for "Archaic"(0), "African"(1) and "Test"(2) population
		//corresponding population information
	std::vector<std::vector<std::string> > idInfos, popInfos;

	//From id to pop label
	std::map<std::string, std::string> popCheck;

	//From id to as pop label
	std::map<std::string, int> ASpopCheck;

	//Outgroup fasta
	std::map<std::string, std::string> outgroupFastaPaths;

	//Ancestral fasta
	std::map<std::string, std::string> ancestralFastaPaths;
private:
	void vcfParLoad();
	void remapParLoad();
	void popParLoad();
	void outgroupFastaPathLoad();
	void ancestralFastaPathLoad();
};

class genome
{
public:
	genome(const genomePar& _par) : par(_par)
	{
		remapLoad(par.remapPaths, par.remapContigs);
		vcfsLoad(par.vcfPaths);
		ancestralLoad(par.ancestralFastaPaths);
	};

	const genomePar& par;

	std::vector<std::string> comIDs, comPop, comContigs;
	std::vector<int> comASPop;
	std::vector<std::vector<long> > comPos;
	std::vector<std::vector<char> > comRefAllele, comAltAllele;
	std::vector<int> phased;
	//[contig] [id]
	std::vector<std::vector<hapInfo> > haplotypes;
	//T means Alt is ancestral allele, F is other
	std::vector<hapInfo> ancestralState;
	std::vector<std::vector<double> > gdisData;

private:
	void vcfsLoad(const std::vector<std::string>& paths);
	void remapLoad(const std::vector<std::string>& paths, \
			const std::vector<std::string>& contigs);
	void ancestralLoad(const std::map<std::string, std::string>& ancestralFastaPaths);

	std::vector<std::vector<double> > gdisMap;
	std::vector<std::vector<long> > posMap;
	std::map<std::string, int> contigMap;
};

#endif /* DATA_HPP_ */
