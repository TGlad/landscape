#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>
#include <cmath>

struct Adj // lower triangular undirected adjacency matrix
{
  void resize(int n)
  {
    data.resize((n*(n+1))/2, 0);
  }
  int size() const
  {
    return (int)std::round((std::sqrt(8.0*(double)data.size() + 1.0) - 1.0)/2.0);
  }
  int operator()(int i, int j) const
  {
    return i >= j ? data[(i * (i + 1) / 2) + j] : data[(j * (j + 1) / 2) + i]; 
  }
  int &operator()(int i, int j)
  {
    return i >= j ? data[(i * (i + 1) / 2) + j] : data[(j * (j + 1) / 2) + i]; 
  }
  std::vector<int> data;    
};