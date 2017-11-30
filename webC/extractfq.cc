#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <tr1/unordered_set>
// use <unordered_set> for more recent g++
using namespace std;

int main(int argc,char **argv)
{
    string line;
    tr1::unordered_set<string> names;
    ifstream in(argv[1],ios::in);
    if(!in.is_open()) return EXIT_FAILURE;
    while(getline(in,line,'\n'))
        names.insert(line);
    in.close();
    string name, seq, name2, qual;
    tr1::unordered_set<string>::iterator r;
    ifstream in2("/dev/stdin", ios::in);
    while(!names.empty())
    {
        if(!getline(in2,name,'\n')) break;
        if(!getline(in2,seq,'\n')) break;
        if(!getline(in2,name2,'\n')) break;
        if(!getline(in2,qual,'\n')) break;
        r=names.find(name);
        if(r==names.end()) continue;
        names.erase(r);
//names are unique we don't need this one anymore
        cout << name << "\n" << seq << "\n" << name2 << "\n" << qual << endl;
    }
    in2.close();
    return 0;
}
