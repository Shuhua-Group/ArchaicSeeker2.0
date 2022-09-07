# include <iostream>
# include <fstream>
# include <map>
# include <sstream>
# include <set>
# include <vector>

using namespace std;

int main(int argc, char **argv)
{
	bool all(false);
	if(argc >= 2)
		all = true;
	ifstream fpind;
	if(!all)
		fpind.open(argv[1]);
	set<string> ids;
	string id;
	vector<int> check;
	string line;
	cout << "#CHR\tPOS\tREF\tALT\tN\tSum" << endl;
	while(getline(cin, line))
		if(line[1] != '#')
			break;
	istringstream head(line);

	if(all)
	{
		while(fpind >> id)
			ids.insert(id);
		for(int i = 0 ; i < 9 ; ++i)
			head >> line;
		while(head >> id)
			if(ids.count(id))
				check.push_back(1);
			else
				check.push_back(0);
	}
	else
	{
		for(int i = 0 ; i < 9 ; ++i)
			head >> line;
		while(head >> id)
			check.push_back(1);
	}
	while(getline(cin, line))
	{
		istringstream cl(line);
		if(line[0] == '#')
			continue;
		string chr, pos, ref, alt;
		cl >> chr >> pos >> line >> ref >> alt >> line >> line >> line >> line;
		int n(0), sum(0);
		for(int i = 0 ; i < check.size(); ++i)
		{
			cl >> line;
			if(check[i])
			{
				if(line[0] == '1')
				{
					++n;
					++sum;
				}
				else if(line[0] == '0')
					++sum;
				if(line[2] == '1')
				{
					++n;
					++sum;
				}
				else if(line[2] == '0')
					++sum;
			}
		}
		cout << chr << '\t' << pos << '\t' << ref << '\t' << alt << '\t' << n << '\t' << sum << endl;
	}
}
