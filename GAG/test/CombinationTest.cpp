#include <iostream>
#include <vector>
#include <string>
#include <GAGPL/MATH/combination.h>

using namespace stdcomb;
using namespace std;

template<class BidIt>
void display(BidIt begin,BidIt end)
{
  for (BidIt it=begin;it!=end;++it)
      cout<<*it<<" ";
  cout<<endl;
}

int main()
{
  vector<int> ca;
  ca.push_back (1);
  ca.push_back (2);
  ca.push_back (3);
  ca.push_back (4);
  ca.push_back (5);
  ca.push_back (6);
  vector<int> cb;
  cb.push_back (1);
  cb.push_back (2);
  cb.push_back (3);
  cb.push_back (4);
   
  do
  {
    display(cb.begin(),cb.end());
  }
  while(next_combination(ca.begin (),ca.end (),cb.begin (),cb.end()) );
  
  cout<<"Complete!"<<endl;
  
  return 0;
}