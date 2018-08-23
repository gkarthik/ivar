#include<iostream>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int main() {
  int num_tests = 6;
  allele a1 = {
    "A",
    40,
    0,
    30
  };
  allele a2 = {
    "T",
    40,
    10,
    30
  };
  allele a3 = {
    "G",
    40,
    10,
    10
  };
  allele a4 = {
    "+T",
    30,
    10,
    10
  };
  allele a5 = {
    "+AG",
    30,
    10,
    10
  };
  allele a6 = {
    "+AC",
    30,
    4,
    10
  };
  allele a7 = {
    "+WT",
    30,
    2,
    10
  };
  int success = 0;
  allele arr[] = {a1,a2,a3,a4,a5};
  ret_t s;
  std::vector<allele> ad(arr, arr+sizeof(arr)/sizeof(allele));
  int size = sizeof(arr)/sizeof(allele);
  for(int i = 0;i<size;i++){
    ad.at(i) = arr[i];
  }
  s = get_consensus_allele(ad, 0, 0, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("DW") == 0) ? 1: 0;
  success += (s.q.compare("8+") == 0) ? 1 : 0;
  ad.push_back(a6);
  s = get_consensus_allele(ad, 0, 0, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("DAS") == 0) ? 1: 0;
  success += (s.q.compare("8++") == 0) ? 1 : 0;
  ad.push_back(a7);
  s = get_consensus_allele(ad, 0, 0, 'N');
  std::cout << s.nuc << ": " << s.q << std::endl;
  success += (s.nuc.compare("DAB") == 0) ? 1: 0;
  success += (s.q.compare("8++") == 0) ? 1 : 0;
  return (success == num_tests) ? 0 : -1;
}
