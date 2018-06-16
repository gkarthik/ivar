#include<iostream>
#include<vector>

#include "primer_bed.h"
#include "htslib/sam.h"

#ifndef trim_primer_quality
#define trim_primer_quality

struct cigar_ {
  uint32_t *cigar;
  uint32_t nlength;
  int32_t start_pos;
};

int trim_bam_qual_primer(std::string bam, std::string bed, std::string bam_out, std::string region_);
int32_t get_pos_on_query(uint32_t *cigar, uint32_t ncigar, int32_t pos, int32_t ref_start);
int32_t get_pos_on_reference(uint32_t *cigar, uint32_t ncigar, uint32_t pos, uint32_t ref_start);
void reverse_qual(uint8_t *q, int l);
void reverse_cigar(uint32_t *cigar, int l);
double mean_quality(uint8_t *a, int s, int e);
cigar_ quality_trim(bam1_t* r, int qual_threshold, int sliding_window);
void print_cigar(uint32_t *cigar, int nlength);
cigar_ primer_trim(bam1_t *r, int32_t new_pos);
void replace_cigar(bam1_t *b, int n, uint32_t *cigar);
int16_t get_overlapping_primer_indice(bam1_t* r, std::vector<primer> primers);
cigar_ remove_trailing_query_ref_consumption(uint32_t* cigar, uint32_t n);
cigar_ condense_cigar(uint32_t* cigar, uint32_t n);

#endif
