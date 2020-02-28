#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "../src/get_common_variants.h"

int main() {
  int num_success = 0;
  std::string files[] = {"../data/test.1.tsv", "../data/test.2.tsv", "../data/test.3.tsv"};
  std::ifstream fin;
  std::map<std::string, std::string> file_tab_delimited_str;
  std::map<std::string, unsigned int> counts;
  // Read all files
  for (int i = 0; i < 3; ++i) {
    fin.open(files[i]);
    read_variant_file(fin, i, counts, file_tab_delimited_str);
    fin.close();
  }
  std::map<std::string, std::string>::iterator it = file_tab_delimited_str.begin();
  // Test tab delimited strings for test.1.tsv.
  std::string keys[] = {
    "test\t320\tA\tT\t0",
    "test\t331\tT\tG\t0",
    "test\t331\tT\tG\t1",
    "test\t365\tA\tT\t0",
    "test\t365\tA\tT\t1",
    "test\t365\tA\tT\t2",
    "test\t42\tG\tT\t0",
    "test\t42\tG\tT\t1",
    "test\t42\tG\tT\t2"
  };
  std::string values[] = {
    "1\t1\t35\t1\t1\t46\t0.5\t2\t0.666667\tFALSE\t",
    "1\t1\t20\t1\t1\t21\t0.5\t2\t0.666667\tFALSE\t",
    "2\t2\t20\t2\t2\t21\t0.5\t2\t0.666667\tFALSE\t",
    "0\t0\t0\t1\t1\t27\t1\t1\t1\tFALSE\t",
    "0\t0\t0\t1\t1\t27\t1\t1\t1\tFALSE\t",
    "0\t0\t0\t1\t1\t22\t1\t1\t1\tFALSE\t",
    "0\t0\t0\t1\t0\t49\t1\t1\t1\tFALSE\t",
    "0\t0\t0\t1\t0\t49\t1\t1\t1\tFALSE\t",
    "0\t0\t0\t1\t0\t49\t1\t1\t1\tFALSE\t"
  };
  int ctr = 0;
  while(it!=file_tab_delimited_str.end()){
    if((it->first).compare(keys[ctr]) != 0 || (it->second).compare(values[ctr]) != 0){
      std::cout << it->first << ": " << it->second << std::endl;
      num_success = -1;
    }
    it++;
    ctr++;
  }
  // Test counts
  std::map<std::string, unsigned int>::iterator count_it = counts.begin();
  std::string count_keys[] = {"test\t320\tA\tT\t", "test\t331\tT\tG\t", "test\t365\tA\tT\t", "test\t42\tG\tT\t"};
  int count_values[] = {1, 2, 3, 3};
  ctr = 0;
  while(count_it!=counts.end()){
    if((count_it->first).compare(count_keys[ctr]) != 0 || count_it->second != count_values[ctr]){
      std::cout << count_it->first << ": " << count_it->second << std::endl;
      num_success = -1;
    }
    count_it++;
    ctr++;
  }  
  if(num_success == 0)
    return 0;
  return -1;
}
