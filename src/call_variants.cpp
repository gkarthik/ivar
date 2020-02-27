#include<stdint.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>
#include<cmath>
#include<htslib/kfunc.h>
#include "htslib/kfunc.h"

#include "allele_functions.h"

const char gap='N';
const float sig_level = 0.01;

std::vector<allele>::iterator get_ref_allele(std::vector<allele> &ad, char ref){
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc[0] == ref)
      return it;
  }
  // print_allele_depths(ad);
  return ad.end();
}

double* get_frequency_depth(allele a, uint32_t pos_depth, uint32_t total_depth){ // pos_depth is ungapped depth after passing quality. total_depth is depth without quality filtering for indels.
  double *val = new double[2];
  if(a.nuc[0] == '+'){	// For insertions use depth discarding quality
    val[0] = a.depth/(double)total_depth;
    val[1] = total_depth;
    return val;
  }
  val[0] = a.depth/(double)pos_depth;
  val[1] = pos_depth;
  return val;
}

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold, uint8_t min_depth){
  std::string line, cell, bases, qualities, region;
  std::ofstream fout((out_file+".tsv").c_str());
  fout << "REGION"
    "\tPOS"
    "\tREF"
    "\tALT"
    "\tREF_DP"
    "\tREF_RV"
    "\tREF_QUAL"
    "\tALT_DP"
    "\tALT_RV"
    "\tALT_QUAL"
    "\tALT_FREQ"
    "\tTOTAL_DP"
    "\tPVAL"
    "\tPASS"
       << std::endl;
  int ctr = 0, pos = 0;
  uint32_t mdepth = 0, pdepth = 0; // mpdepth for mpileup depth and pdeth for ungapped depth at position
  double pval_left, pval_right, pval_twotailed, *freq_depth, err;
  std::stringstream lineStream;
  char ref;
  std::vector<allele> ad;
  std::vector<allele>::iterator ref_it;
  while (std::getline(cin, line)){
    lineStream << line;
    ctr = 0;
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	region = cell;
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
    if(mdepth < min_depth) {	// Check for minimum depth
      lineStream.clear();
      continue;
    }
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    if(ad.size() == 0){
      lineStream.clear();
      continue;
    }
    ref_it = get_ref_allele(ad, ref);
    if (ref_it == ad.end()) {	// If ref not present in reads.
      allele a;
      a.nuc = ref;
      a.depth = 0;
      a.reverse = 0;
      a.mean_qual = 0;
      ad.push_back(a);
      ref_it = ad.end() - 1;
    }
    // Get ungapped coverage
    pdepth = 0;
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if(it->nuc[0]=='*' || it->nuc[0] == '+' || it->nuc[0] == '-')
	continue;
      pdepth += it->depth;
    }
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if((*it == *ref_it) || it->nuc[0]=='*')
	continue;
      freq_depth = get_frequency_depth(*it, pdepth, mdepth);
      if(freq_depth[0] < min_threshold)
	continue;
      fout << region << "\t";
      fout << pos << "\t";
      fout << ref << "\t";
      fout << it->nuc << "\t";
      fout << ref_it->depth << "\t";
      fout << ref_it->reverse << "\t";
      fout << (uint16_t)ref_it->mean_qual << "\t";
      fout << it->depth << "\t";
      fout << it->reverse << "\t";
      fout << (uint16_t) it->mean_qual << "\t";
      fout << freq_depth[0] << "\t";
      fout << freq_depth[1] << "\t";
      /*
	    | Var   | Ref      |
	Exp | Error | Err free |
	Obs | AD    | RD       |
       */
      err = pow(10, ( -1 * (it->mean_qual)/10));
      kt_fisher_exact((err * mdepth), (1-err) * mdepth, it->depth, ref_it->depth, &pval_left, &pval_right, &pval_twotailed);
      fout << pval_left << "\t";
      if(pval_left <= sig_level){
	fout << "TRUE";
      } else {
	fout << "FALSE";
      }
      fout << std::endl;
    }
    lineStream.clear();
  }
  fout.close();
  return 0;
}
