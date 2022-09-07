# include <iostream>
# include <vector>
# include <fstream>
# include <sstream>

using namespace std;

int main()
{
	cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	string line;
	getline(cin,line);
	istringstream head(line);
	head >> line >> line >> line >> line >> line;
	int nind(0);
	while(head >> line)
		cout << '\t' << line;
	cout << endl;
	while(getline(cin,line))
	{
		istringstream cl(line);
		string chr;
		long start, end;
		double gstart, gend;
		cl >> chr >> start >> end >> gstart >> gend;
		cout << chr << '\t' << start << "\t.\t.\t.\t.\tGenetDis:" << gstart << "-" << gend <<  "\tEND=" << end - 1 << "\tGT";
		while(cl >> line)
			cout << '\t' << line;
		cout << endl;
	}
}
