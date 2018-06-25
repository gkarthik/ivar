#include <iostream>

#include "../src/allele_functions.h"

int main(){
  int num_tests = 2;
  std::string b = "AAAATTTG+3ATGT-3ATG";
  int size = 9, success = 0, i;
  char ref = 'A';
  uint8_t _q[] = {30,30,20,23,20,15,20,40,20};
  std::string q = "";
  for (int i = 0; i < size; ++i) {
    _q[i] += 33;
    q += _q[i];
  }
  std::cout << q << std::endl;
  std::vector<allele> ad = update_allele_depth(ref, b, q, 10);
  print_allele_depths(ad);
  i = find_ref_in_allele(ad, 'T');
  success += (ad.at(i).mean_qual == 18) ? 1 : 0;
  i = find_ref_in_allele(ad, 'A');
  success += (ad.at(i).mean_qual == 25) ? 1 : 0;
  return (num_tests == success) ? 0 : -1;
}
