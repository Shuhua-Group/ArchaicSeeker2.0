# include <iostream>
# include <fstream>
# include <vector>
# include <sstream>
# include <zlib.h>
# include <string>

using namespace std;

const int MAX = 4096 * 16;

int main(int argc, char **argv)
{
	gzFile fdaf;
	fdaf = gzopen(argv[1], "rb");
	string pre(argv[2]);
	ofstream fpo(argv[3]);
	vector<vector<string> > anc;
	anc.resize(23);
	char buff[MAX];
	for(int i = 1; i <= 22; ++i)
	{
		ostringstream path;
		path << pre << "_" << i << ".fa.gz";
		gzFile fdanc;
		fdanc = gzopen(path.str().c_str(), "rb");
		gzgets(fdanc, buff, MAX);
		while(gzgets(fdanc, buff, MAX))
		{
			string line(buff);
			anc[i].push_back(line);
//			cout << line.size() << endl;
		}
		gzclose(fdanc);
	}
	int n(anc[1].front().size() - 1);
	fpo << "CHR\tPOS\tREF\tALT\tANC\tN\tSum" << endl;
	gzgets(fdaf, buff, MAX);
	while(gzgets(fdaf, buff, MAX))
	{
		istringstream cl(buff);
		int chr;
		long pos;
		string ref, alt, num, sum;
		cl >> chr >> pos >> ref >> alt >> num >> sum;
		char a(anc[chr][(pos - 1) / n][(pos - 1) % n]);
		//char b(anc[chr][pos / n][pos % n]);

		fpo << chr << '\t' << pos << '\t' << ref << '\t' << alt << '\t' << a << '\t' << num << '\t' << sum << endl;
	}
}
