# include <iostream>
# include <vector>
# include <sstream>
# include <fstream>
# include <map>
# include <cstring>
# include <iomanip>
# include <set>
# include <zlib.h>

using namespace std;

const int MAX = 4096 * 256;

class gmap
{
	public:
		string gmapPath;
		gmap(string gmapPath): gmapPath(gmapPath)
		{
			pos.resize(23);
			gdis.resize(23);
			mem.resize(23);
			for(int i = 1 ; i <= 22 ; ++i)
			{
				ostringstream path;
				path << gmapPath << "/genetic_map_chr" << i << "_combined_b37.txt";
				ifstream fpi(path.str().c_str());
				vector<long> &curPos(pos[i]);
				vector<double> &curGdis(gdis[i]);
				curPos.clear();
				curGdis.clear();
				string line;
				fpi >> line >> line >> line;
				long tpos;
				double tgdis;
				while(fpi >> tpos >> line >> tgdis)
				{
					curPos.push_back(tpos);
					curGdis.push_back(tgdis);
				}
			}
		}
		double search(int chr, long p)
		{
			if(mem[chr].count(p))
				return mem[chr][p];
			const vector<long>& curPos(pos[chr]);
			const vector<double>& curGdis(gdis[chr]);
			int start(0), end(curPos.size() - 1);
			int mid;
			while(end - start > 1)
			{
				mid = (end + start) / 2;
				if(curPos[mid] > p)
					end = mid;
				else if(curPos[mid] < p)
					start = mid;
				else
				{
					mem[chr][p] = curGdis[mid];
					return curGdis[mid];
				}
			}
			while(start < curPos.size() && curPos[start] < p)
				++start;
			if(start == curPos.size())
			{
				mem[chr][p] = curGdis.back();
				return curGdis.back();
			}
			else if(start <= 1)
			{
				mem[chr][p] = curGdis.front();
				return curGdis.front();
			}
			--start;
			double r = curGdis[start] + ( p - curPos[start]	) * ( curGdis[start + 1] - curGdis[start] ) / ( curPos[start + 1] - curPos[start] );
			mem[chr][p] = r;
			return r;
		}
		vector<map<long,double> > mem;
		vector<vector<long> > pos;
		vector<vector<double> > gdis;
};

class bed
{
	public:
		int chr;
		long start, end;
		double gstart, gend;
		string id;
};

int main(int argc, char **argv)
{
	string gmapPath(argv[1]);
	gmap g(gmapPath);
	ifstream fpi(argv[2]);
	string out(argv[3]);
	set<string> check;
	for(int i = 3; i < argc ; ++i)
		check.insert(argv[i]);
	string line;
	getline(fpi,line);
	map<string, vector<vector<bed> > > segs;
	while(getline(fpi, line))
	{
		istringstream cl(line);
		string lab;
		bed tbed;
		cl >> tbed.id >> tbed.chr >> tbed.start >> tbed.end >> line >> line >> lab;
		if(!check.count(lab))
			continue;
		tbed.gstart = g.search(tbed.chr, tbed.start);
		tbed.gend = g.search(tbed.chr, tbed.end);
		segs[tbed.id].resize(23);
		segs[tbed.id][tbed.chr].push_back(tbed);
	}
	string path;
	path = out + ".seg";
	ofstream fpo(path.c_str());
	fpo << setprecision(20);

	for(map<string, vector<vector<bed> > >::iterator it = segs.begin(); it != segs.end(); ++it)
	{
		for(int i = 1 ; i <= 22 ; ++i)
		{
			const vector<bed> & curData(it->second.at(i));
			if(curData.size() ==0)
				continue;
			double gstart(g.gdis[i].front());
			long start(g.pos[i].front());
			for(int j = 0 ; j < curData.size(); ++j)
			{
				fpo << gstart / 100 << '\t' << curData[j].gstart / 100 << "\tModern" << endl;
				fpo << curData[j].gstart / 100 << '\t' << curData[j].gend / 100 << "\tArchaic" << endl;
				gstart = curData[j].gend;
				start = curData[j].end;
			}
			fpo << gstart / 100 << '\t' << g.gdis[i].back() / 100 << "\tModern" << endl;
		}
	}
}


