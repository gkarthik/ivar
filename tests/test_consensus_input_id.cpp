#include<iostream>
#include<fstream>
#include<string>
#include "../src/call_consensus_pileup.h"
#include "../src/allele_functions.h"

int call_cns_check_outfile(std::string input_id, std::string prefix, std::string cns, char gap, bool call_min_depth, int min_depth){
  std::string path = "../data/test.gap.sorted.mpileup";
  std::string expctd_hdr = "";
  std::ifstream mplp(path);
  call_consensus_from_plup(mplp, input_id, prefix, 20, 0, min_depth, gap, call_min_depth);
  std::ifstream outFile(prefix+".fa");
  std::string l;
  getline(outFile, l);		// header
  if(input_id.empty()) {
    char *o = new char[prefix.length() + 1];
    strcpy(o, prefix.c_str());
    expctd_hdr = ">Consensus_" + std::string(basename(o)) + "_threshold_0_quality_20";
    
  } else {
    expctd_hdr = ">" + input_id;
  }

  return l.compare(expctd_hdr);
}

int main() {
  int num_success = 0;
  std::string c_ = "CTGCTGGGTCATGGGCCCATCATGATGGTCTTGGCGATTCTAGCCTTTTTGAGATTCACGGCAATCAAGCCATCACTGGGTCTCATCAATAGATGGGGTTCAGTGGGGAAAAAAGAGGCTATGGAAACAATAAAGAAGTTCAAGAAAGAT------------------------------AGGAAGGAGAAGAAGAGACGWGGCGCAGATACTAGTGTCGGAATTGTTGGMCTCCTGCTGACCACAGCTATGGMAGCGGAGGTCACKAGACGTGGGAGTGCATACTATATGTACTTGGACWGAAACGATGCKGGGGAGGCCATATCTTTTCCAACCACATTGGGGTTGAATAAGTG";
  num_success = call_cns_check_outfile("TESTID", "../data/test.gap", c_, '-', true, 0);
  std::cout << num_success << std::endl;
  num_success += call_cns_check_outfile("", "../data/test.gap", c_, 'N', true, 0);
  std::cout << num_success << std::endl;
  num_success += call_cns_check_outfile("TESTID2", "../data/test.gap", c_, 'N', false, 0);
  std::cout << num_success << std::endl;
  if(num_success == 0)
    return 0;
  return -1;
}