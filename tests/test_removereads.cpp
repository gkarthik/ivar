#include <iostream>
#include <fstream>
#include <vector>
#include "../src/remove_reads_from_amplicon.h"

#include "htslib/sam.h"

int main(int argc, char *argv[]){
  int num_success = 1, num_tests = 1;
  std::vector<std::string> amp;
  std::ifstream fin("../data/test.masked_primer_indices.txt");
  std::string s, region = "Consensus_ZI-27_threshold_0_quality_20";
  while(getline(fin, s, '\t' ) ){
    amp.push_back(s);
  }
  rmv_reads_from_amplicon("../data/test.trimmed.bam", region, "../data/test.trimmed.masked.bam", amp, "../data/test.bed", "ivar removereads -i ../data/test.trimmed.bam -p ../data/test.masked -t ../data/test.masked_primer_indices.txt -b ../data/test.bed ");
  std::string out_file = "../data/test.trimmed.masked.bam";
  samFile *in = hts_open(out_file.c_str(), "r");
  bam_hdr_t *header = sam_hdr_read(in);
  hts_itr_t *iter = NULL;
  hts_idx_t *idx = sam_index_load(in, out_file.c_str());
  iter  = sam_itr_querys(idx, header, region.c_str());
  bam1_t *aln = bam_init1();
  int ctr = 0;
  bool w;
  while(sam_itr_next(in, iter, aln) >= 0) {
    uint8_t* a = bam_aux_get(aln, "XA");
    w = true;
    if(a != 0){
      if(bam_aux2i(a) == 1 || bam_aux2i(a) == 3 || bam_aux2i(a) == 4){
	w = false;
      }
    }
    if(!w)
      break;
  }
  if(w)
    num_success += 1;
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  return (num_success == num_tests) ? 0 : 1;
}
