# include <iostream>
# include <fstream>
# include <vector>
# include <sstream>

using namespace std;

int main()
{
	string line, tmp;
	getline(cin,line);
	cout << line << endl;
	while(getline(cin,line))
	{
		if(line[0] == 'C')
			continue;
		istringstream cl(line);
		cl >> tmp >> tmp >> tmp >> tmp >> tmp;
		while(cl >> tmp)
			if(tmp[0] == '1' || tmp[2] == '1')
			{
				cout << line << endl;
				break;
			}
	}
}
