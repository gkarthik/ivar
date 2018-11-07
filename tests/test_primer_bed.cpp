#include <iostream>
#include <vector>

#include "../src/primer_bed.h"

int main(){
  int num_tests = 2;
  int success = 0;
  std::vector<primer>::iterator it;
  std::vector<primer> primers = populate_from_file("../data/test.bed");
  std::string primer_names[] = {"WNV_400_1_LEFT", "WNV_400_2_LEFT", "WNV_400_1_RIGHT", "WNV_400_2_RIGHT", "WNV_400_3_LEFT"};
  unsigned int primer_start[] = {8, 249, 359,658,569};
  int flag = 1;
  for(it = primers.begin(); it != primers.end(); ++it) {
    if(it->get_name() != primer_names[it-primers.begin()] || it->get_start() != primer_start[it-primers.begin()])
      flag = 0;
  }
  success += flag;
  populate_pair_indices(primers, "../data/pair_information.tsv");
  int pair_indices[] = {2,3,0,1,-1};
  flag = 1;
  for(it = primers.begin(); it != primers.end(); ++it) {
    if(it->get_pair_indice() != pair_indices[it-primers.begin()])
      flag = 0;
  }
  success += flag;
  return (num_tests == success) ? 0 : -1;
}
