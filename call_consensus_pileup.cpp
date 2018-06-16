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
    } else if (it->nuc[0] == '+'){
      it->nuc[0] = it->prev_base;
      prev_ind = find_ref_in_allele(ad, it->prev_base);
      if(prev_ind!=-1)
	ad.at(prev_ind).depth -= it->depth;
    }
    ++it;
  }
}

std::string get_consensus_allele(std::vector<allele> ad){
  if(ad.size()==0)
    return "N";
  print_allele_depths(ad);
  format_alleles(ad);
  print_allele_depths(ad);
  if(ad.size() == 1)
    return (ad.at(0).nuc.compare("*") == 0) ? "" : ad.at(0).nuc;
  std::string cnuc = "";
  char n;
  int max_l = 0, mdepth = 0, tdepth = 0;
  uint32_t gaps = 0;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc.length() > max_l){
      max_l = it->nuc.length();
    }
  }
  std::cout << max_l << std::endl;
  for (int i = 0; i < max_l; ++i){
    n = '*';
    mdepth = 0;
    tdepth = 0;
    std::vector<allele>::iterator it = ad.begin();
    gaps = 0;
    while(it!=ad.end()){
      if(!(i < it->nuc.length())){
	gaps += it->depth;
	it++;
	continue;
      }
      tdepth = it->depth;
      while(it!=ad.end() -1 && it->nuc[i] == (it+1)->nuc[i] && (i+1) < (it+1)->nuc.length()){ // Third condition for cases with more than 1 base in insertion
	tdepth += (it+1)->depth;
	it++;
      }
      if(tdepth > mdepth){
	n = it->nuc[i];
	mdepth = tdepth;
      } else if(tdepth == mdepth){
	n = gt2iupac(n, it->nuc[i]);
      }
      it++;
    }
    if(n!='*' && mdepth >= gaps) // TODO: Check what to do when equal.
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
