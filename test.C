#include<vector>
#include <iostream>

using namespace std;
void foo(const std::vector<int>& v)
{
  std::vector< std::vector<int> > vectors;
  vectors.push_back(v);
  vectors.back()[0] =2;
  vectors.push_back(v);
}

int main()
{
  std::vector<int> v={1,2,3,4};
  cout <<v.size()<<endl;
  std::vector< std::vector<int> > vectors;
  vectors.push_back(v);
  cout <<v.size()<<endl;
  foo(v);
  cout <<v.size()<<endl;

}
