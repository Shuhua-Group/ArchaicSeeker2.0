# include <iostream>
# include <fstream>
# include <vector>
# include <sstream>
# include <map>
# include <set>

using namespace std;

class bed
{
	public:
		long start, end;
		double gstart, gend;
};

int main(int argc, char **argv)
{
	ifstream fpind(argv[1]);
	map<string, int> idCheck;
	vector<string> ids;
	string tid, tpop;
	while(fpind >> tid)
	{
		idCheck[tid + "_1"] = ids.size();
		ids.push_back(tid + "_1");
		idCheck[tid + "_2"] = ids.size();
		ids.push_back(tid + "_2");
	}
	int nind = ids.size();
	vector<map<string,vector<bed> > >data;
	data.resize(nind);
	ifstream fpi(argv[2]);
	string line;
	getline(fpi,line);
	set<string> labs;
	map<string, set<long> > pos;
	string chr;
	map<long,double> gmap;
	while(getline(fpi,line))
	{
		istringstream cl(line);
		bed t;
		string lab;
		cl >> tid >> chr >> t.start >> t.end >> t.gstart >> t.gend >> lab;
		if(lab == "YRI_Modern_Denisova_Neanderthal" || lab == "Denisova_Neanderthal")
			lab = "Others";
		else if(lab != "Denisova" && lab != "Neanderthal")
			continue;
		if(idCheck.count(tid) == 0)
			continue;
		labs.insert(lab);
		pos[lab].insert(t.start);
		pos[lab].insert(t.end);
		gmap[t.start] = t.gstart;
		gmap[t.end] = t.gend;
		data[idCheck[tid]][lab].push_back(t);
	}
	for(set<string>::iterator it = labs.begin(); it != labs.end(); ++it)
	{
		if(pos.count(*it) == 0)
			continue;
		vector<long> curPos(pos[*it].begin(), pos[*it].end());
		int npos = curPos.size();
		vector<vector<short> > mat;
		mat.resize(npos - 1);
		for(int i = 0 ; i < npos - 1; ++i)
			mat[i].resize(nind, 0);
		for(int i = 0 ; i < nind; ++i)
		{
			if(data[i].count(*it) == 0)
				continue;
			vector<bed>& curBed(data[i][*it]);
			int p(0);
			int nbed(curBed.size());
			for(int j = 0 ; j < npos - 1; ++j)
			{
				if(p == nbed)
					break;
				long start = curPos[j];
				long end = curPos[j + 1];
				if(curBed[p].start > end)
					continue;
				else if(curBed[p].start <= start && curBed[p].end >=end)
					mat[j][i] = 1;
				if(curBed[p].end == end)
					++p;
			}
		}
		string out(argv[3]);
		string path = out + "_" + *it + ".txt";
		ofstream fpo(path.c_str());
		fpo << "Contig\tStart(bp)\tEnd(bp)\tStart(cM)\tEnd(cM)";
		for(int i = 0 ; i < nind ; i += 2)
			fpo << '\t' << ids[i].substr(0,ids[i].size() - 2);
		fpo << endl;
		for(int i = 0 ; i < npos - 1; ++i)
		{
			fpo << chr << '\t' << curPos[i] << '\t' << curPos[i + 1] << '\t' << gmap[curPos[i]] << '\t' << gmap[curPos[i+1]];
			for(int j = 0 ; j < nind ; j += 2)
				fpo << '\t' << mat[i][j] << "|" << mat[i][j + 1];
			fpo << endl;
		}
	}
}
