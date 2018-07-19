#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>

#include "allele_functions.h"
#include "call_consensus_pileup.h"

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

bool compare_allele_depth(const allele &a, const allele &b){
  return b.depth < a.depth;
}

ret_t get_consensus_allele(std::vector<allele> ad, uint8_t min_qual, double threshold){
  ret_t t;
  t.nuc="";
  t.q = min_qual+33;
  if(ad.size()==0)
    return t;
  format_alleles(ad);
  if(ad.size() == 1){
    t.nuc = (ad.at(0).nuc.compare("*") == 0) ? "" : ad.at(0).nuc;
    t.q = (ad.at(0).nuc.compare("*") == 0) ? min_qual + 33 : ad.at(0).mean_qual + 33;
    return t;
  }
  std::string cnuc = "";
  std::vector<allele> nuc_pos;
  allele tmp_a;
  char n;
  uint32_t max_l = 0, max_depth = 0, tmp_depth = 0, cur_depth = 0, total_max_depth = 0, gap_depth = 0, total_indel_depth = 0;
  uint8_t ambg_n = 1, ctr = 0;
  double q = 0, tq = 0, cur_threshold = 0;
  std::vector<allele>::iterator it;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc.length() > max_l){
      max_l = it->nuc.length();
    }
    if(it->nuc.length() == 1){
      total_max_depth += it->depth;
      total_indel_depth += it->depth;
      total_indel_depth = total_indel_depth - it->end; // For indels
    }
  }
  t.nuc = "";
  t.q = "";
  for (uint16_t i = 0; i < max_l; ++i){
    n = '*';
    q = 0;
    tq = 0;
    max_depth = 0;
    tmp_depth = 0;
    cur_depth = 0;
    // prev_depth = 0;
    ctr = 1;
    it = ad.begin();
    gap_depth = 0;
    nuc_pos.clear();
    while(it!=ad.end()){
      if(!(i < it->nuc.length())){
	it++;
	continue;
      }
      if(it->nuc[i] == '+' || it->nuc[i] == '-'){
	it++;
	continue;
      }
      tmp_depth = it->depth;
      tq = it->mean_qual;
      ctr = 1;
      while(it!=ad.end()-1 && i < (it+1)->nuc.length() && it->nuc[i] == (it+1)->nuc[i]){ // Third condition for cases with more than 1 base in insertion  && (i+1) < (it+1)->nuc.length()
	tmp_depth += (it+1)->depth;
	tq = ((tq * (ctr)) + (it+1)->mean_qual)/(ctr + 1);
	it++;
	ctr++;
      }
      cur_depth += tmp_depth;
      tmp_a.nuc = it->nuc[i];
      tmp_a.mean_qual = tq;
      tmp_a.reverse = 0;	// Reverse reads not important for consensus
      tmp_a.depth = tmp_depth;
      nuc_pos.push_back(tmp_a);
      // if(tmp_depth > max_depth){
      // 	n = it->nuc[i];
      // 	max_depth = tmp_depth;
      // 	q = tq;
      // 	ambg_n = 1;
      // } else if(tmp_depth == max_depth){
      // 	n = gt2iupac(n, it->nuc[i]);
      // 	q = ((q * ambg_n) + tq)/(ambg_n+1);
      // 	ambg_n += 1;
      // }
      it++;
    }
    // Sort nuc_pos by depth of alleles
    std::sort(nuc_pos.begin(), nuc_pos.end(), compare_allele_depth);
    it=nuc_pos.begin();
    n = it->nuc[0];
    q = it->mean_qual;
    max_depth = it->depth;
    it++;
    ambg_n = 1;
    if(i > 0){
      cur_threshold = threshold * (double) total_indel_depth;
    } else {
      cur_threshold =  threshold * (double)total_max_depth;
    }
    while(it!=nuc_pos.end() && (max_depth < cur_threshold || it->depth == (it -1)->depth)){ // Iterate till end or till cross threshold and no more allele with same depth
      n = gt2iupac(n, it->nuc[0]);
      q = ((q * ambg_n) + it->mean_qual)/(ambg_n+1);
      ambg_n += 1;
      max_depth += it->depth;
      it++;
    }
    if(max_depth < cur_threshold) // If depth still less than threshold
      n = 'N';
    if(i > 0)
      gap_depth = (total_indel_depth > cur_depth) ? total_indel_depth - cur_depth : 0;
    else
      gap_depth = 0;			  // For first position of allele
    if(n!='*' && max_depth >= gap_depth){ // TODO: Check what to do when equal.{
      t.nuc += n;
      q += 0.5;			// For rounding before converting to int
      t.q += (((uint8_t)q)+33);
    }
  }
  return t;
}

int call_consensus_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double threshold){
  std::string line, cell;
  std::ofstream fout(out_file+".fa");
  std::ofstream tmp_qout(out_file+".qual.txt");
  fout << ">Consensus" << "_" << threshold <<std::endl;
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
    ret_t t = get_consensus_allele(ad, min_qual, threshold);
    fout << t.nuc;
    tmp_qout << t.q;
    lineStream.clear();
    ad.clear();
  }
  tmp_qout.close();
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
