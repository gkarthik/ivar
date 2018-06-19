#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>

#include "allele_functions.h"

const std::string gap="N";

void format_alleles(std::vector<allele> &ad){
  std::vector<allele>::iterator it = ad.begin();
  int prev_ind = 0;
  while(it < ad.end()){
    if(it->nuc[0] == '-'){	// Remove deletions
      it = ad.erase(it);
    } // else if (it->nuc[0] == '+'){
    //   it->nuc[0] = it->prev_base;
    //   prev_ind = find_ref_in_allele(ad, it->prev_base);
    //   if(prev_ind!=-1)
    // 	ad.at(prev_ind).depth -= it->depth;
    // }
    ++it;
  }
}

std::string get_consensus_allele(std::vector<allele> ad){
  if(ad.size()==0)
    return "N";
  format_alleles(ad);
  print_allele_depths(ad);
  if(ad.size() == 1)
    return (ad.at(0).nuc.compare("*") == 0) ? "" : ad.at(0).nuc;
  std::string cnuc = "";
  char n;
  uint32_t max_l = 0, max_depth = 0, tmp_depth = 0, cur_depth = 0, prev_depth = 0;
  int32_t gap_depth = 0;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc.length() > max_l){
      max_l = it->nuc.length();
    }
  }
  std::cout << max_l << std::endl;
  for (int i = 0; i < max_l; ++i){
    n = '*';
    max_depth = 0;
    tmp_depth = 0;
    cur_depth = 0;
    prev_depth = 0;
    std::vector<allele>::iterator it = ad.begin();
    gap_depth = 0;
    while(it!=ad.end()){
      if(!(i < it->nuc.length())){
	prev_depth += it->depth;
	it++;
	continue;
      }
      if(it->nuc[i] == '+'){
	it++;
	continue;
      }
      tmp_depth = it->depth;
      while(it<=ad.end() -1 && it->nuc[i] == (it+1)->nuc[i]){ // Third condition for cases with more than 1 base in insertion  && (i+1) < (it+1)->nuc.length()
	tmp_depth += (it+1)->depth;
	it++;
      }
      cur_depth += tmp_depth;
      if(tmp_depth > max_depth){
	n = it->nuc[i];
	max_depth = tmp_depth;
      } else if(tmp_depth == max_depth){
	n = gt2iupac(n, it->nuc[i]);
	std::cout << "Nuc: " << n << std::endl;
      }
      it++;
    }
    gap_depth = prev_depth - cur_depth;
    std::cout << max_depth << " " << prev_depth << " " << cur_depth <<  " " << gap_depth << " " << std::endl;
    if(n!='*' && max_depth >= gap_depth) // TODO: Check what to do when equal.
      cnuc += n;
  }
  std::cout << "Cns: " << cnuc << std::endl;
  return cnuc;
}

int call_consensus_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual){
  std::string line, cell;
  std::ofstream fout(out_file+".fa");
  fout << ">Consensus"<<std::endl;
  int ctr = 0, pos = 0, mdepth = 0;
  std::stringstream lineStream;
  char ref;
  std::string bases;
  std::string qualities;
  std::vector<allele> ad;
  while (std::getline(cin, line)){
    lineStream << line;
    ctr = 0;
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	ref = cell[0];
	break;
      case 3:
	mdepth = stoi(cell);
	break;
      case 4:
	bases = cell;
	break;
      case 5:
	qualities = cell;
	break;
      case 6:
	break;
      }
      ctr++;
    }
    std::cout << pos << std::endl;
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    fout << get_consensus_allele(ad);
    lineStream.clear();
  }
  fout.close();
  return 0;
}

// int main_old(int argc, char* argv[]) {
//   std::string line, cell;
//   std::string out_file = argv[1];
//   uint8_t min_qual = 20;
//   if(argc > 1){
//     min_qual = atoi(argv[2]);
//   }
//   std::ofstream fout(out_file+".fa");
//   fout << ">Consensus"<<std::endl;
//   int ctr = 0, pos = 0, mdepth = 0;
//   std::stringstream lineStream;
//   char ref;
//   std::string bases;
//   std::string qualities;
//   std::vector<allele> ad;
//   while (std::getline(std::cin, line)){
//     lineStream << line;
//     ctr = 0;
//     while(std::getline(lineStream,cell,'\t')){
//       switch(ctr){
//       case 0:
// 	break;
//       case 1:
// 	pos = stoi(cell);
// 	break;
//       case 2:
// 	ref = cell[0];
// 	break;
//       case 3:
// 	mdepth = stoi(cell);
// 	break;
//       case 4:
// 	bases = cell;
// 	break;
//       case 5:
// 	qualities = cell;
// 	break;
//       case 6:
// 	break;
//       }
//       ctr++;
//     }
//     ad = update_allele_depth(ref, bases, qualities, min_qual);
//     // print_allele_depths(ad);
//     fout << get_consensus_allele(ad);
//     // std::cout << std::endl << "Consensus: " << get_consensus_allele(ad) << std::endl;
//     lineStream.clear();
//   }
//   fout.close();
//   return 0;
// }
