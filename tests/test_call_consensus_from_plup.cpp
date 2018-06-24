#include<iostream>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int main() {
  int num_tests = 2;
  allele a1 = {
    "A",
    30,
    0,
    30
  };
  allele a2 = {
    "T",
    30,
    10,
    30
  };
  allele a3 = {
    "G",
    30,
    10,
    10
  };
  allele a4 = {
    "GT",
    30,
    10,
    10
  };
  allele a5 = {
    "AAG",
    30,
    10,
    10
  };
  int success = 0;
  allele arr[] = {a1,a2,a3,a4,a5};
  std::vector<allele> ad(arr,arr+5);
  std::string c0 = get_consensus_allele(ad, 0);
  std::cout << c0 << std::endl;
  success += (c0.compare("DWG") == 0) ? 1: 0;
  std::string c20 = get_consensus_allele(ad, 20);
  std::cout << c0 << std::endl;
  success += (c20.compare("W") == 0) ? 1: 0;
  return (success == num_tests) ? 0 : -1;
}
