/*
 * data.cpp
 *
 *  Created on: Sep 29, 2018
 *      Author: yuankai
 */

# include "data.hpp"
# include "gzfstream.hpp"

# include <unistd.h>
# include <cassert>
# include <set>
# include <sstream>

const unsigned char bitTable[] = {1, 2, 4, 8, 16, 32, 64, 128};

const int byteCountTable[] =
{
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

hapInfo::hapInfo(const hapInfo& rhs) : nBlock(rhs.nBlock), nByte(rhs.nByte), \
		capBlock(rhs.capBlock), nInfo(rhs.nInfo)
{
	mem = new unsigned char[capBlock];
	memcpy(mem, rhs.mem, capBlock);
}

hapInfo::hapInfo(int s) : mem(NULL)
{
	nBlock = s >> 3;
	nByte = s % 8;
	nInfo = s;
	capBlock = nBlock + 1;
	mem = new unsigned char[capBlock];
	memset(mem, 0, capBlock);
}

void hapInfo::set(int n)
{
	assert(n < nInfo && n >= 0 );
	int block(n >> 3), bit(n % 8);
	mem[block] |= bitTable[bit];
}

void hapInfo::set()
{
	memset(mem, -1, nBlock);
	if(nByte)
		mem[nBlock] = -1;
}

void hapInfo::reset(int n)
{
	assert(n < nInfo && n >= 0);
	int block(n >> 3), bit(n % 8);
	mem[block] &= ~bitTable[bit];
}

void hapInfo::reset()
{
	memset(mem, 0, nBlock);
	if(nByte)
		mem[nBlock] = 0;
}

void hapInfo::flip()
{
	for(int i = 0 ; i < nBlock ; ++i)
		mem[i] = ~mem[i];
	if(nByte)
		mem[nBlock] = ~mem[nBlock];
}

int hapInfo::count() const
{
	int ret(0);
	for(int i = 0 ; i < nBlock; ++i)
		ret += byteCountTable[mem[i]];
	if(nByte)
	{
		unsigned char mask = -1;
		mask <<= (8 - nByte);
		mask >>= (8 - nByte);
		ret += byteCountTable[mem[nBlock] & mask];
	}
	return ret;
}

int hapInfo::count(int start, int end) const
{
	int ret(0);
	--end;
	if(start < 0)
		start = 0;
	if(end > nInfo)
		end = nInfo;
	int startBlock(start >> 3), startByte(start % 8);
	int endBlock(end >> 3), endByte(end % 8);
	if(startBlock != endBlock)
	{
		unsigned char mask(-1), tmp;
		mask >>= startByte;
		mask <<= startByte;
		tmp = mem[startBlock] & mask;
		ret += byteCountTable[tmp];
		mask = -1;
		mask <<= (7 - endByte);
		mask >>= (7 - endByte);
		tmp = mem[endBlock] & mask;
		ret += byteCountTable[tmp];
	}
	else
	{
		unsigned char mask(-1), tmp;
		mask >>= startByte;
		mask <<= startByte + (7 - endByte);
		mask >>= 7 - endByte;
		tmp = mem[startBlock] & mask;
		ret += byteCountTable[tmp];
	}
	for(int i = startBlock + 1; i < endBlock; ++i)
		ret += byteCountTable[mem[i]];
	return ret;
}

int hapInfo::size() const
{
	return nInfo;
}

void hapInfo::resize(int s)
{
	assert(s >= 0);
	int block = s >> 3;
	if(block > capBlock)
	{
		unsigned char *newMem;
		capBlock = block + 1;
		newMem = new unsigned char[capBlock];
		memset(newMem, 0, capBlock);
		if(nByte)
			memcpy(newMem, mem, nBlock + 1);
		else
			memcpy(newMem, mem, nBlock);
		delete []mem;
		mem = newMem;
	}
	nBlock = block;
	nByte = s % 8;
	nInfo = s;
}

void hapInfo::reserve(int s)
{
	assert(s >= 0);
	int block = s >> 3;
	if(block > capBlock)
	{
		unsigned char *newMem;
		capBlock = block + 1;
		newMem = new unsigned char[capBlock];
		memset(newMem, 0, capBlock);
		if(nByte)
			memcpy(newMem, mem, nBlock + 1);
		else
			memcpy(newMem, mem, nBlock);
		delete []mem;
		mem = newMem;
	}
}

hapInfo hapInfo::operator &(const hapInfo& rhs) const
{
	assert(nInfo == rhs.size());
	hapInfo ret(nInfo);
	for(int i = 0 ; i < nBlock ; ++i)
		ret.mem[i] = mem[i] & rhs.mem[i];
	if(nByte)
		ret.mem[nBlock] = mem[nBlock] & rhs.mem[nBlock];
	return ret;
}

void hapInfo::operator &=(const hapInfo& rhs)
{
	assert(nInfo == rhs.size());
	for(int i = 0 ; i < nBlock ; ++i)
		mem[i] &= rhs.mem[i];
	if(nByte)
		mem[nBlock] &= rhs.mem[nBlock];
}

hapInfo hapInfo::operator |(const hapInfo& rhs) const
{
	assert(nInfo == rhs.size());
	hapInfo ret(nInfo);
	for(int i = 0 ; i < nBlock ; ++i)
		ret.mem[i] = mem[i] | rhs.mem[i];
	if(nByte)
		ret.mem[nBlock] = mem[nBlock] | rhs.mem[nBlock];
	return ret;
}

void hapInfo::operator |=(const hapInfo& rhs)
{
	assert(nInfo == rhs.size());
	for(int i = 0 ; i < nBlock ; ++i)
		mem[i] |= rhs.mem[i];
	if(nByte)
		mem[nBlock] |= rhs.mem[nBlock];
}

hapInfo hapInfo::operator ^(const hapInfo& rhs) const
{
	assert(nInfo == rhs.size());
	hapInfo ret(nInfo);
	for(int i = 0 ; i < nBlock ; ++i)
		ret.mem[i] = mem[i] ^ rhs.mem[i];
	if(nByte)
		ret.mem[nBlock] = mem[nBlock] ^ rhs.mem[nBlock];
	return ret;
}

void hapInfo::operator ^=(const hapInfo& rhs)
{
	assert(nInfo == rhs.size());
	for(int i = 0 ; i < nBlock ; ++i)
		mem[i] ^= rhs.mem[i];
	if(nByte)
		mem[nBlock] ^= rhs.mem[nBlock];
}

bool hapInfo::operator [](int n) const
{
	assert(n < nInfo && n>=0);
	int block(n >> 3), byte(n % 8);
	return mem[block] & bitTable[byte];
}

hapInfo& hapInfo::operator =(const hapInfo& rhs)
{
	if(this == &rhs)
		return *this;
	if(capBlock < rhs.capBlock)
	{
		delete[] mem;
		mem = new unsigned char[rhs.capBlock];
		capBlock = rhs.capBlock;
	}
	memcpy(mem, rhs.mem, rhs.capBlock);
	nBlock = rhs.nBlock;
	nByte = rhs.nByte;
	nInfo = rhs.nInfo;
	return *this;
}

hapInfo hapInfo::operator~() const
{
	hapInfo ret(*this);
	ret.flip();
	return ret;
}

int hapDiff(const hapInfo& a, const hapInfo& b)
{
	assert(a.size() == b.size());
	int ret(0);
	unsigned char tmp;
	int block(a.nBlock);
	for(int i = 0 ; i < block; ++i)
	{
		tmp = a.mem[i] ^ b.mem[i];
		ret += byteCountTable[tmp];
	}
	if(a.nByte)
	{
		unsigned char mask(-1);
		mask <<= (8 - a.nByte);
		mask >>= (8 - a.nByte);
		tmp = ( a.mem[block] ^ b.mem[block] ) & mask;
		ret += byteCountTable[tmp];
	}
	return ret;
}

int hapDiff(const hapInfo& a, const hapInfo& b, int start, int end)
{
	assert(a.size() == b.size() && end > start && start >= 0);
	--end;
	int ret(0);
	int startBlock(start >> 3), startByte(start % 8);
	int endBlock(end >> 3), endByte(end % 8);
	if(startBlock != endBlock)
	{
		unsigned char mask(-1), tmp;
		mask >>= startByte;
		mask <<= startByte;
		tmp = ( a.mem[startBlock] ^ b.mem[startBlock] ) & mask;
		ret += byteCountTable[tmp];
		mask = -1;
		mask <<= (7 - endByte);
		mask >>= (7 - endByte);
		tmp = ( a.mem[endBlock] ^ b.mem[endBlock] ) & mask;
		ret += byteCountTable[tmp];
	}
	else
	{
		unsigned char mask(-1), tmp;
		mask >>= startByte;
		mask <<= startByte + (7 - endByte);
		mask >>= 7 - endByte;
		tmp = ( a.mem[startBlock] ^ b.mem[startBlock] ) & mask;
		ret += byteCountTable[tmp];
	}
	for(int i = startBlock + 1; i < endBlock; ++i)
	{
		unsigned char tmp;
		tmp = a.mem[i] ^ b.mem[i];
		ret += byteCountTable[tmp];
	}
	return ret;
}

void genomePar::vcfParLoad()
{
	assert(-1 != access(vcfPar.c_str(), R_OK));
	gzifstream fpi(vcfPar);
	std::string path;
	fpi >> path;
	assert(path == "vcf");
	while(fpi >> path)
		vcfPaths.push_back(path);
	assert(vcfPaths.size());
}

void genomePar::remapParLoad()
{
	assert(-1 != access(remapPar.c_str(), R_OK));
	gzifstream fpi(remapPar);
	std::string path, contig;
	fpi >> path >> contig;
	assert(path == "remap");
	assert(contig == "contig");
	while(fpi >> path >> contig)
	{
		remapPaths.push_back(path);
		remapContigs.push_back(contig);
	}
	assert(remapPaths.size() && remapContigs.size());
}

void genomePar::popParLoad()
{
	assert(-1 != access(popPar.c_str(), R_OK));
	gzifstream fpi(popPar);
	std::string id, pop, asPop;
	fpi >> id >> pop >> asPop;
	assert(id == "ID");
	assert(pop == "Pop");
	assert(asPop == "ArchaicSeekerPop");
	idInfos.resize(3);
	popInfos.resize(3);
	while(fpi >> id >> pop >> asPop)
	{
		int asLab(-1);
		if(asPop == "Archaic")
			asLab = 0;
		else if(asPop == "African")
			asLab = 1;
		else if(asPop == "Test")
			asLab = 2;
		assert(asLab >= 0);
		if(ASpopCheck.count(pop))
			assert(ASpopCheck[pop] == asLab);
		else
			ASpopCheck[pop] = asLab;
		assert(popCheck.count(id) == 0);
		popCheck[id] = pop;
		idInfos[asLab].push_back(id);
		popInfos[asLab].push_back(pop);
	}
	for(int i = 0 ; i < 3; ++i)
	{
		assert(idInfos[i].size());
		assert(popInfos[i].size());
	}
}

void genomePar::outgroupFastaPathLoad()
{
	assert(-1 != access(outgroupPar.c_str(), R_OK));
	gzifstream fpi(outgroupPar);
	std::string path, contig;
	fpi >> path >> contig;
	assert(path == "outgroup");
	assert(contig == "contig");
	while(fpi >> path >> contig)
		outgroupFastaPaths[contig] = path;
	assert(outgroupFastaPaths.size());
}

void genomePar::ancestralFastaPathLoad()
{
	assert(-1 != access(ancestralPar.c_str(), R_OK));
	gzifstream fpi(ancestralPar);
	std::string path, contig;
	fpi >> path >> contig;
	assert(path == "ancestor");
	assert(contig == "contig");
	while(fpi >> path >> contig)
		ancestralFastaPaths[contig] = path;
	assert(ancestralFastaPaths.size());
}

class comSNP
//Only used in vcfsLoad().
{
public:
	comSNP() : ref(' '), alt(' '), findT(false), findF(false), count(0), pos(0) {}
	char ref, alt;
	bool findT, findF;
	short count;
	long pos;
};

void genome::vcfsLoad(const std::vector<std::string>& paths)
{
	assert(paths.size());
	int npath = paths.size();
	std::map<std::string, std::map<std::string, int> > idContigFileCheck;
	std::map<std::string, std::vector<comSNP> > consistSNPCount;
	std::map<std::string, int> contigFileCount;
	for(int i = 0 ; i < npath; ++i)
	{
		std::cout << "Loading SNP information from " << paths[i] << std::endl;
		assert(-1 != access(paths[i].c_str(), R_OK));
		gzifstream fpi(paths[i]);
		std::string line;
		while(getline(fpi, line))
			if(line[1] != '#')
				break;
		std::istringstream head(line);
		std::vector<std::string> tids;
		std::string tid;
		for(int j = 0 ; j < 9 ; ++j)
			head >> tid;
		std::set<std::string> idRepCheck;
		while(head >> tid)
		{
			tids.push_back(tid);
			assert(idRepCheck.count(tid) == 0);
			assert(par.popCheck.count(tid));
			idRepCheck.insert(tid);
		}
		int nind(tids.size());
		long tpos, prePos(0);
		std::string tchr, preChr(""), ref, alt;
		std::set<std::string> checkContigs;
		int p(0);
		std::vector<comSNP>* curSNPs(NULL);
		bool newAdd(false);
		int curNsnp(0);
		while(getline(fpi, line))
		{
			std::istringstream cl(line);
			cl >> tchr >> tpos >> line >> ref >> alt >> line >> line >> line >> line;
			bool findT(false), findF(false);
			while(cl >> line && !(findT && findF))
			{
				if(line[0] == '0' || line[2] == '0')
					findF = true;
				if(line[0] == '1' || line[2] == '1')
					findT = true;
			}
			if(contigMap.count(tchr) == 0)
				continue;
			if(preChr != tchr)
			{
				assert(checkContigs.count(tchr) == 0);
				checkContigs.insert(tchr);
				if(consistSNPCount.count(tchr) == 0)
					newAdd = true;
				else
				{
					newAdd = false;
					curNsnp = consistSNPCount[tchr].size();
				}
				curSNPs = & (consistSNPCount[tchr]);
				p = 0;
			}
			else
				assert(tpos > prePos);
			prePos = tpos;
			preChr = tchr;
			if(newAdd)
			{
				if(ref.size() == 1 && alt.size() == 1)
				{
					comSNP t;
					t.ref = ref[0];
					t.alt = alt[0];
					t.count = 1;
					t.pos = tpos;
					t.findF = findF;
					t.findT = findT;
					curSNPs->push_back(t);
				}
			}
			else
			{
				if(p == curNsnp)
					continue;
				while(p < curNsnp && tpos > (*curSNPs)[p].pos)
					++p;
				if(p == curNsnp || tpos < (*curSNPs)[p].pos)
					continue;
				if(ref.size() == 1 && alt.size() == 1 && \
						ref[0] == (*curSNPs)[p].ref)
				{
					if(alt[0] == (*curSNPs)[p].alt || alt[0] == '.')
					{
						++(*curSNPs)[p].count;
						(*curSNPs)[p].findF |= findF;
						(*curSNPs)[p].findT |= findT;
					}
					else if((*curSNPs)[p].alt == '.')
					{
						++(*curSNPs)[p].count;
						(*curSNPs)[p].alt = alt[0];
						(*curSNPs)[p].findF |= findF;
						(*curSNPs)[p].findT |= findT;
					}
				}
			}
		}
		for(std::set<std::string>::iterator it = checkContigs.begin(); \
			it != checkContigs.end(); ++it)
		{
			++contigFileCount[*it];
			for(int j = 0; j < nind; ++j)
			{
				const std::string& curID = tids[j];
				assert(idContigFileCheck[curID].count(*it) == 0);
				idContigFileCheck[curID][*it] = i;
			}
		}
	}
#ifdef DEBUG_VCF_INTER
	std::cout << "idContigFileCheck [ID] [contig] [file]" << std::endl;
	for(std::map<std::string, std::map<std::string, int> >::iterator \
			it = idContigFileCheck.begin(); it != idContigFileCheck.end(); ++it)
	{
		std::cout << "\tID: " << it->first << std::endl;
		std::map<std::string, int>& contigFile = it->second;
		for(std::map<std::string, int>::iterator tit = \
				contigFile.begin(); tit != contigFile.end(); ++tit)
			std::cout << "\t\t" << tit->first << '\t' << tit->second << std::endl;
	}
	std::cout << "consistSNPCount [contig] [position] [comSNP]" << std::endl;
	for(std::map<std::string, std::vector<comSNP> >::iterator \
			it = consistSNPCount.begin(); it != consistSNPCount.end(); ++it)
	{
		std::cout << "\tID: " << it->first << std::endl;
		std::vector<comSNP>& snpCount = it->second;
		for(std::vector<comSNP>::iterator tit = \
				snpCount.begin(); tit != snpCount.end(); ++tit)
			std::cout << "\t\t" << tit->pos << '\t' << tit->ref << '\t' \
			<< tit->alt << '\t' << tit->count << std::endl;
	}
	std::cout << "contigFileCount [contig] [int]" << std::endl;
	for(std::map<std::string, int>::iterator it = contigFileCount.begin(); \
			it != contigFileCount.end(); ++it)
		std::cout << "\t" << it->first << '\t' << it->second << std::endl;
#endif//DEBUG_VCF_INTER

	int nind(idContigFileCheck.size());
	int ncontig(contigFileCount.size());
	std::vector<std::vector<char> > idContigCheckTable;
	std::vector<int> contigCountCheck;
	idContigCheckTable.resize(nind);
	for(int i = 0 ; i < nind ; ++i)
		idContigCheckTable[i].resize(ncontig, 0);
	contigCountCheck.resize(ncontig, 0); //should be = nind

	comIDs.resize(nind);
	comPop.resize(nind);
	comASPop.resize(nind);
	std::map<std::string, int> comIDCheckTable;
	int idc(0), contigc(0);
	for(std::map<std::string, std::map<std::string, int> >::iterator \
			it = idContigFileCheck.begin(); it != idContigFileCheck.end(); ++it)
	{
		contigc = 0;
		comIDs[idc] = it->first;
		comPop[idc] = par.popCheck.at(comIDs[idc]);
		comASPop[idc] = par.ASpopCheck.at(comPop[idc]);
		comIDCheckTable[it->first] = idc;
		for(std::map<std::string, int>::iterator tit = contigFileCount.begin(); \
				tit != contigFileCount.end(); ++tit)
		{
			if(it->second.count(tit->first))
			{
				idContigCheckTable[idc][contigc] = 1;
				++contigCountCheck[contigc];
			}
			++contigc;
		}
		++idc;
	}
	std::map<std::string, int> comContigCheckTable;
	contigc = 0;
	for(std::map<std::string, int>::iterator it = contigFileCount.begin(); \
			it != contigFileCount.end(); ++it)
	{
		if(contigCountCheck[contigc] == nind)
		{
			comContigCheckTable[it->first] = comContigs.size();
			comContigs.push_back(it->first);
		}
		++contigc;
	}
	ncontig = comContigs.size();
	comPos.resize(ncontig);
	comRefAllele.resize(ncontig);
	comAltAllele.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::vector<long>& curPos = comPos[i];
		std::vector<char>& curRefAllele = comRefAllele[i];
		std::vector<char>& curAltAllele = comAltAllele[i];
		std::vector<comSNP>& curContig = consistSNPCount[comContigs[i]];
		int nfile = contigFileCount[comContigs[i]];
		for(std::vector<comSNP>::iterator it = curContig.begin(); it != curContig.end(); \
				++it)
			if(it->count == nfile && it->findF && it->findT)
			{
				curPos.push_back(it->pos);
				curRefAllele.push_back(it->ref);
				curAltAllele.push_back(it->alt);
			}
	}
#ifdef DEBUG_VCF_INTER
	std::cout << "comIDs\tcomPop\tcomASPop" << std::endl;
	for(int i = 0 ; i < nind; ++i)
		std::cout << comIDs[i] << '\t' << comPop[i] << '\t' << comASPop[i] << std::endl;
	std::cout << "comConigs & comPos" << std::endl;
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::cout << comContigs[i] << std::endl;
		int npos = comPos[i].size();
		for(int j = 0 ; j < npos; ++j)
			std::cout << '\t' << comPos[i][j] << std::endl;
	}
	std::cout << "idContigCheckTable" << std::endl;
	for(std::map<std::string, int>::iterator tit = contigFileCount.begin(); \
			tit != contigFileCount.end(); ++tit)
		std::cout << '\t' << tit->first ;
	std::cout << std::endl;
	for(int i = 0 ; i < nind ; ++i)
	{
		std::cout << comIDs[i];
		for(unsigned int j = 0 ; j < idContigCheckTable[i].size(); ++j)
			std::cout << '\t' << (int)idContigCheckTable[i][j];
		std::cout << std::endl;
	}
#endif//DEBUG_VCF_INTER

	haplotypes.resize(ncontig);
	phased.resize(nind, 0);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::vector<hapInfo>& curHap = haplotypes[i];
		curHap.resize(nind * 2);
		int nmk(comPos[i].size());
		for(int j = 0 ; j < nind * 2; ++j)
		{
			curHap[j].resize(nmk);
			curHap[j].reset();
		}
	}
	for(int i = 0 ; i < nind; ++i)
		if(comASPop[i] == 2)
			phased[i] = 1;
		else
			phased[i] = 0;
	for(int i = 0 ; i < npath ; ++i)
	{
		std::cout << "Loading genotype information from " << paths[i] << std::endl;
		gzifstream fpi(paths[i]);
		std::string line;
		while(getline(fpi, line))
			if(line[1] != '#')
				break;
		std::istringstream head(line);
		std::vector<std::string> tids;
		std::string tid;
		for(int j = 0 ; j < 9 ; ++j)
			head >> tid;
		std::vector<int> check;
		int isphased(-1);
		while(head >> tid)
		{
			int index(comIDCheckTable[tid]);
			check.push_back(index);
			if(isphased == -1)
				isphased = phased[index];
			else
				assert(isphased == phased[index]);
		}
		int tnind = check.size();
		std::string pchr(""), tchr;
		long tpos;
		int p(0);
		std::vector<hapInfo>* curHap(NULL);
		std::vector<long>* curPos(NULL);
		int curNpos(1024);
		if(isphased)
			while(getline(fpi,line) && p < curNpos)
			{
				std::istringstream cl(line);
				cl >> tchr >> tpos;
				if(pchr != tchr || curHap == NULL)
				{
					if(comContigCheckTable.count(tchr) == 0)
						continue;
					p = 0;
					pchr = tchr;
					int chrIndex(comContigCheckTable[tchr]);
					curHap = &haplotypes[chrIndex];
					curPos = &comPos[chrIndex];
					curNpos = curPos->size();
				}
				if(tpos < (*curPos)[p])
					continue;
				assert(tpos == (*curPos)[p]);
				for(int j = 0 ; j < 7 ; ++j)
					cl >> line;
				for(int j = 0 ; j < tnind; ++j)
				{
					cl >> line;
					int idp(check[j] * 2);
					if(line[0] == '1')
						(*curHap)[idp].set(p);
					else
						assert(line[0] == '0');
					if(line[2] == '1')
						(*curHap)[idp + 1].set(p);
					else
						assert(line[2] == '0');
					assert(line[1] == '|');
				}
				++p;
			}
		else
			while(getline(fpi,line) && p < curNpos)
			{
				std::istringstream cl(line);
				cl >> tchr >> tpos;
				if(pchr != tchr || curHap == NULL)
				{
					if(comContigCheckTable.count(tchr) == 0)
						continue;
					p = 0;
					pchr = tchr;
					int chrIndex(comContigCheckTable[tchr]);
					curHap = &haplotypes[chrIndex];
					curPos = &comPos[chrIndex];
					curNpos = curPos->size();
				}
				if(tpos < (*curPos)[p])
					continue;
				assert(tpos == (*curPos)[p]);
				for(int j = 0 ; j < 7 ; ++j)
					cl >> line;
				for(int j = 0 ; j < tnind; ++j)
				{
					cl >> line;
					int n(0), s(2);
					int idp(check[j] * 2);
					if(line[0] == '1')
						++n;
					else if(line[0] == '.')
						--s;
					else
						assert(line[0] == '0');
					if(line[2] == '1')
						++n;
					else if(line[2] == '.')
						--s;
					else
						assert(line[2] == '0');
					if(s != 2)
					{
						(*curHap)[idp].set(p);
						(*curHap)[idp + 1].set(p);
					}
					else if(n == 1)
						(*curHap)[idp].set(p);
					else if(n == 2)
						(*curHap)[idp + 1].set(p);
					else
						assert(n == 0);
				}
				++p;
			}
	}
#ifdef DEBUG_VCF_INTER
	std::cout << "haplotypes" << std::endl;
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::cout << "\tcontig: " << comContigs[i] << std::endl;
		for(int j = 0; j < nind * 2; ++j)
		{
			std::cout << "\t\t" << j << " ";
			haplotypes[i][j].simplePrint();
		}
	}
#endif//DEBUG_VCF_INTER

	gdisData.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::vector<long>& curPosMap = posMap[contigMap[comContigs[i]]];
		std::vector<double>& curGdisMap = gdisMap[contigMap[comContigs[i]]];
		std::vector<long>& curComPos = comPos[i];
		std::vector<double>& curGdisData = gdisData[i];
		int ncomPos(curComPos.size()), nmap(curPosMap.size());
		curGdisData.resize(ncomPos);
		int pMap(0), pData(0);
		while(pData < ncomPos && curComPos[pData] <= curPosMap.front())
		{
			curGdisData[pData] = 0;
			++pData;
		}
		if(pData == ncomPos)
			continue;
		while(pData < ncomPos)
		{
			while(pMap < nmap && curPosMap[pMap] <= curComPos[pData])
				++pMap;
			if(pMap == nmap)
				break;
			curGdisData[pData] = ( curGdisMap[pMap] - curGdisMap[pMap - 1] ) / \
					( curPosMap[pMap] - curPosMap[pMap - 1] ) * \
					( curComPos[pData] - curPosMap[pMap - 1] ) + curGdisMap[pMap - 1];
			++pData;
		}
		while(pData < ncomPos)
		{
			curGdisData[pData] = curGdisMap.back();
			++pData;
		}
	}
#ifdef DEBUG_VCF_INTER
	std::cout << "gdisData" << std::endl;
	for(int i = 0 ; i < ncontig; ++i)
	{
		std::cout << "\tcontig: " << comContigs[i] << std::endl;
		std::cout << "\t map:" << std::endl;
		for(unsigned int j = 0 ; j < gdisMap[i].size(); ++j)
			std::cout << "\t\t" << posMap[i][j] << '\t' << gdisMap[i][j] << std::endl;
		std::cout << "\t data:" << std::endl;
		int nmk(gdisData[i].size());
		for(int j = 0 ; j < nmk; ++j)
			std::cout << "\t\t" << comPos[i][j] << '\t' << gdisData[i][j] << std::endl;
	}
#endif//DEBUG_VCF_INTER
}

void genome::remapLoad(const std::vector<std::string>& paths, \
		const std::vector<std::string>& contigs)
{
	assert(paths.size() == contigs.size());
	int npath = paths.size();
	gdisMap.resize(npath);
	posMap.resize(npath);
	for(int i = 0 ; i < npath; ++i)
	{
		assert(-1 != access(paths[i].c_str(), R_OK));
		gzifstream fpi(paths[i]);
		std::string line;
		getline(fpi, line);
		long tpos;
		double trate, tgdis;
		std::vector<long>& curPos = posMap[i];
		std::vector<double>& curGdis = gdisMap[i];
		while(fpi >> tpos >> trate >> tgdis)
		{
			curPos.push_back(tpos);
			curGdis.push_back(tgdis);
		}
		assert(curPos.size() && curGdis.size());
		contigMap[contigs[i]] = i;
	}
}

void genome::ancestralLoad(const std::map<std::string, std::string>& ancestralFastaPaths)
{
	int ncontig(comContigs.size());
	ancestralState.resize(ncontig);
	for(int i = 0 ; i < ncontig; ++i)
	{
		assert(ancestralFastaPaths.count(comContigs[i]));
		const std::string& curPath(ancestralFastaPaths.at(comContigs[i]));
		assert(-1 != access(curPath.c_str(), R_OK));
		gzifstream fpi(curPath);
		std::string line;
		getline(fpi, line);
		int p(0);
		const std::vector<long>& curPos(comPos[i]);
		const std::vector<char>& curAlt(comAltAllele[i]);
		int nsite(curPos.size());
		hapInfo& curAnc(ancestralState[i]);
		curAnc.resize(nsite);
		curAnc.reset();
		int total(0);
		while(p < nsite)
		{
			while(total < curPos[p])
			{
				getline(fpi,line);
				if(!fpi)
					break;
				total += line.size();
			}
			if(!fpi)
				break;
			if(curAlt[p] == line[curPos[p] - total + line.size() - 1])
				curAnc.set(p);
			++p;
		}
	}
}












