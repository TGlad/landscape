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
/*  bool operator==(Adj &other)
  {
    if (hash == "")
    {
      hash = get_wl_hash(); // dynamic programming / caching
    }
    if (other.hash == "")
    {
      other.hash = other.get_wl_hash();
    }
    return hash == other.hash;
  }*/
  std::vector<int> data;
//  std::string hash;

  std::string get_wl_hash(int iterations = 3) 
  {
    int n = size();
    std::vector<std::string> labels(n, "1"); // Initial colors (all same)

    for (int iter = 0; iter < iterations; ++iter) 
    {
        std::vector<std::string> next_labels(n);
        
        for (int i = 0; i < n; ++i) 
        {
            std::vector<std::string> neighbor_labels;
            for (int j = 0; j < n; ++j) 
            {
                if ((*this)(i,j) > 0) 
                { // If edge exists
                    neighbor_labels.push_back(labels[j]);
                }
            }
            // Sort neighbor labels to ensure the representation is unique
            std::sort(neighbor_labels.begin(), neighbor_labels.end());

            // Create a new label based on: (own_label, {sorted_neighbor_labels})
            std::stringstream ss;
            ss << labels[i] << ",";
            for (const auto& s : neighbor_labels) ss << s;
            next_labels[i] = ss.str();
        }

        // Relabel with smaller strings to prevent memory explosion (compression)
        std::map<std::string, std::string> compression_map;
        int id = 0;
        std::vector<std::string> sorted_next = next_labels;
        std::sort(sorted_next.begin(), sorted_next.end());
        
        for (const auto& s : sorted_next) 
        {
            if (compression_map.find(s) == compression_map.end()) 
            {
                compression_map[s] = std::to_string(id++);
            }
        }

        for (int i = 0; i < n; ++i) 
        {
            labels[i] = compression_map[next_labels[i]];
        }
    }

    // The final graph signature is the sorted multiset of the final labels
    std::sort(labels.begin(), labels.end());
    std::string final_hash = "";
    for (const auto& l : labels) 
      final_hash += l + "|";
    return final_hash;
  }    
};