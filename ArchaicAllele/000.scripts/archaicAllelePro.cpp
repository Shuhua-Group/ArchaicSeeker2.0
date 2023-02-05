# include <iostream>
# include <fstream>
# include <vector>
# include <sstream>
# include <zlib.h>

using namespace std;

class snp
{
	public:
		long pos;
		char ref, alt, anc;
		int nafr, sumafr;
};

const int MAX = 4096 * 16;

double power(0.95);

int main(int argc, char **argv)
{
	gzFile fdt, fdafr;
	fdt = gzopen(argv[1], "rb");
	fdafr = gzopen(argv[2], "rb");
	ofstream fpo(argv[3]);
	vector<vector<snp> > data;
	data.resize(23);
	char buff[MAX];
	gzgets(fdafr, buff, MAX);
	snp tsnp;
	while(gzgets(fdafr, buff, MAX))
	{
		istringstream cl(buff);
		int chr;
		cl >> chr >> tsnp.pos >> tsnp.ref >> tsnp.alt >> tsnp.anc >> tsnp.nafr >> tsnp.sumafr;
		data[chr].push_back(tsnp);
	}
	fpo << "CHR\tPOS\tREF\tALT\tANC\tAncConfid\tAfrDerFreq\tArchaicDerFreq\tModDerFreq\tArchaicAlleleProb" << endl;
	gzgets(fdt, buff, MAX);
	vector<int> count;
	count.resize(23, 0);
	while(gzgets(fdt, buff, MAX))
	{
		istringstream cl(buff);
		int chr;
		long pos;
		char ref, alt;
		int nas, sumas, nmod, summod, sum;
		cl >> chr >> pos >> ref >> alt >> nas >> sumas >> nmod >> summod >> sum;
		if(sumas == 0 || nas == 0)
			continue;
		while(count[chr] < data[chr].size() && data[chr][count[chr]].pos < pos)
			++count[chr];
		if(data[chr][count[chr]].pos > pos || count[chr] == data[chr].size())
			continue;
		if(data[chr][count[chr]].pos == pos)
		{
			snp& cur(data[chr][count[chr]]);
			int state(0);
			if(cur.ref == ref && cur.alt == alt)
			{
				char anc(cur.anc);
				anc &= ~' ';
				if(anc == cur.ref)
					state = 0;
				else if(anc == cur.alt)
					state = 1;
				else
					state = -9;
				fpo << chr << '\t' << pos << '\t' << ref << '\t' << alt << '\t' << cur.anc << '\t' << state;
				double afrf, archf, modf, pro;
				afrf = (double) cur.nafr / cur.sumafr;
				archf = (double) nas / sumas;
				modf = (double) nmod / summod;
				if(state == 1)
				{
					afrf = 1 - afrf;
					archf = 1 - archf;
					modf = 1 - modf;
				}
				if(nmod == 0)
					modf = 0;
				pro = ( modf * ( 1 - power ) + 1 - modf ) * ( afrf * 0.003 + 1 - afrf );
				fpo << '\t' << afrf << '\t' << archf << '\t' << modf << '\t' << pro << endl;
			}
		}
	}
}
