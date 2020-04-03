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
  if (idx == NULL) {
    if (sam_index_build2(bam.c_str(), 0, 0) < 0) {
      std::cerr << ("Unable to open BAM/SAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
      if (idx == NULL) {
        std::cerr << "Unable to create BAM/SAM index." << std::endl;
        return -1;
      }
    }
  }
  bam_hdr_t *header = sam_hdr_read(in);
  region_.assign(header->target_name[0]);
  std::string temp(header->text);
  hts_itr_t *iter = NULL;
  iter  = sam_itr_querys(idx, header, region_.c_str());
  bam1_t *aln = bam_init1();
  cigar_ t;
  uint32_t *cigar;
  int primer_ctr = 0;
  int primer_indices[] = {0, 0, 7, 7, 6};
  int cigar_flag[5][3] = {{BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CSOFT_CLIP, BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP}, {BAM_CMATCH, BAM_CSOFT_CLIP, BAM_CSOFT_CLIP}, {BAM_CSOFT_CLIP, BAM_CMATCH}};
  int cigar_len[5][3] = {{11, 139}, {24, 11, 115}, {121, 23, 6}, {103, 23, 24}, {23, 127}};
  int condense_cigar_flag[5][3] = {{BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CSOFT_CLIP, BAM_CMATCH}, {BAM_CMATCH, BAM_CSOFT_CLIP}, {BAM_CMATCH, BAM_CSOFT_CLIP}, {BAM_CSOFT_CLIP, BAM_CMATCH}};
  int condense_cigar_len[5][3] = {{11, 139}, {35, 115}, {121, 29}, {103, 47}, {23, 127}};
  int overlapping_primer_sizes[] = {0, 2, 2, 0, 0, 0, 0, 2, 2, 1};
  int ctr = 0;
  std::vector<primer> overlapping_primers;
  primer cand_primer;
  while(sam_itr_next(in, iter, aln) >= 0) {
    if((aln->core.flag&BAM_FUNMAP) != 0){
      continue;
    }
    std::cout << bam_get_qname(aln) << std::endl;
    get_overlapping_primers(aln, primers, overlapping_primers);
    if(overlapping_primers.size() != overlapping_primer_sizes[ctr]){
      success = -1;
      std::cout << "Overlappingprimer sizes for " << bam_get_qname(aln)  <<". Expected: " << overlapping_primer_sizes[ctr] << ". Got: " << overlapping_primers.size() << std::endl;
    }
    if(overlapping_primers.size() > 0){
      if(bam_is_rev(aln)){
	cand_primer = get_min_start(overlapping_primers);
	std::cout << primer_ctr << ": " << cand_primer.get_start() << std::endl;
	t = primer_trim(aln, cand_primer.get_start() - 1);
      } else {
	cand_primer = get_max_end(overlapping_primers);
	t = primer_trim(aln, cand_primer.get_end() + 1);
      }
      if(cand_primer.get_indice() != primer_indices[primer_ctr]){
	success = -1;
	std::cout << "Primer indice wrong. Expected: " << primer_indices[primer_ctr] << ". Got: " << cand_primer.get_indice() << std::endl;
      }
      replace_cigar(aln, t.nlength, t.cigar);
      cigar = bam_get_cigar(aln);
      for (int i = 0; i < t.nlength; ++i){
	if(((cigar[i]) & BAM_CIGAR_MASK) != cigar_flag[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar flag didn't match for " << cand_primer.get_indice()  <<  " ! Expected " << cigar_flag[primer_ctr][i]  << " " << "Got " << ((cigar[i]) & BAM_CIGAR_MASK) << std::endl;
	}
	if((((cigar[i]) >> BAM_CIGAR_SHIFT)) != cigar_len[primer_ctr][i]){
	  success = -1;
	  std::cout << "Cigar length didn't match for " << bam_get_qname(aln)  <<  " ! Expected " << cigar_len[primer_ctr][i]  << " " << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
	  std::cout << i << ": " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
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
	  std::cout << "Cigar length didn't match after condense! Expected " << condense_cigar_len[primer_ctr][i]  << " " << "Got " << ((cigar[i]) >> BAM_CIGAR_SHIFT) << std::endl;
	}
      }
      primer_ctr++;
    }
    ctr++;
    std::cout << " ---- " << std::endl;
  }
  // Check if primers found at all
  success = (primer_ctr > 0) ? success : -1;
  return success;
}
