#include<iostream>
#include "../src/trim_primer_quality.h"
#include "../src/primer_bed.h"
#include "htslib/sam.h"

int main(){
  int success = 0;
  std::string bam = "../data/test.unmapped.sorted.bam";
  std::vector<primer> primers = populate_from_file("../data/test.bed");
  std::string region_;
  samFile *in = hts_open(bam.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  bam_hdr_t *header = sam_hdr_read(in);
  region_.assign(header->target_name[0]);
  std::string temp(header->text);
  hts_itr_t *iter = NULL;
  iter  = sam_itr_querys(idx, header, region_.c_str());
  bam1_t *aln = bam_init1();
  cigar_ t;
  uint32_t *cigar;
  int ctr = 0, primer_ctr = 0;
  int primer_indices[] = {5, 0, 0, 5, 5, 5, 5, 2, 2, 1};
  int cigar_flag[5][3] = {{BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP}, {BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP}, {BAM_CSOFT_CLIP, BAM_CMATCH}};
  int cigar_len[5][3] = {{11, 139}, {24, 11, 115}, {130, 14, 6}, {112, 14, 24}, {20, 130}};
  int condense_cigar_flag[5][3] = {{BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CMATCH, BAM_CSOFT_CLIP}, {BAM_CMATCH, BAM_CSOFT_CLIP}, {BAM_CSOFT_CLIP, BAM_CMATCH}};
  int condense_cigar_len[5][3] = {{11, 139}, {35, 115}, {130, 20}, {112, 38}, {20, 130}};
  uint8_t primer_indice = 0;
  while(sam_itr_next(in, iter, aln) >= 0) {
    primer_indice = get_overlapping_primer_indice(aln, primers);
    if(primer_indice != primer_indices[ctr]){
      success = -1;
      std::cout << "Primer indice wrong. Expected: " << primer_indices[ctr] << ". Got: " << primer_indice << std::endl;
    }
    if(primer_indice < primers.size()){
      if(bam_is_rev(aln)){
	t = primer_trim(aln, primers[primer_indice].get_start() + 1);
      } else {
	t = primer_trim(aln, primers[primer_indice].get_end() + 1);
      }
      replace_cigar(aln, t.nlength, t.cigar);
      cigar = bam_get_cigar(aln);
      for (int i = 0; i < t.nlength; ++i){
	if(((cigar[i]) & BAM_CIGAR_MASK) != cigar_flag[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar flag didn't match! Expected " << cigar_flag[primer_ctr][i]  << " " << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
	}
	if((((cigar[i]) >> BAM_CIGAR_SHIFT)) != cigar_len[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar flag didn't match! Expected " << cigar_len[primer_ctr][i]  << " " << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
	}
      }
      // Check condense
      t = condense_cigar(t.cigar, t.nlength);
      replace_cigar(aln, t.nlength, t.cigar);
      cigar = bam_get_cigar(aln);
      for (int i = 0; i < t.nlength; ++i){
	if(((cigar[i]) & BAM_CIGAR_MASK) != condense_cigar_flag[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar flag didn't match! Expected " << condense_cigar_flag[primer_ctr][i]  << " " << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
	}
	if((((cigar[i]) >> BAM_CIGAR_SHIFT)) != condense_cigar_len[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar flag didn't match after condense! Expected " << condense_cigar_len[primer_ctr][i]  << " " << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
	}
      }
      primer_ctr++;
    }
    ctr++;
  }
  return success;
}
