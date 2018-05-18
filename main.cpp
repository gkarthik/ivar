#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "primer_bed.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

struct cigar_ {
  uint32_t *cigar;
  uint32_t nlength;
  int32_t start_pos;
};

int32_t get_pos_on_query(uint32_t *cigar, uint32_t ncigar, int32_t pos, int32_t ref_start){
  int cig;
  int32_t n;
  int32_t ql = 0, rl = ref_start;
  for (uint32_t i = 0; i < ncigar; ++i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (bam_cigar_type(cig) & 2) { // Only reference consuming
      if (pos <= rl + n) {
	if (bam_cigar_type(cig) & 1) // Only query consuming
	  ql += (pos - rl);	   // n consumed reference, check if it consumes query too.
	return ql;
      }
      rl += n;
    }
    if (bam_cigar_type(cig) & 1) // Only query consuming
      ql += n;
  }
  return ql;
}

int32_t get_pos_on_reference(uint32_t *cigar, uint32_t ncigar, uint32_t pos, uint32_t ref_start){
  int cig;
  int32_t n;
  uint32_t ql = 0, rl = ref_start;
  for (uint32_t i = 0; i < ncigar; ++i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (bam_cigar_type(cig) & 1) { // Only query consuming
      if (pos <= ql + n) {
	if (bam_cigar_type(cig) & 2) // Only reference consuming
	  rl += (pos - ql);	   // n consumed reference, check if it consumes query too.
	return rl;
      }
      ql += n;
    }
    if (bam_cigar_type(cig) & 2) // Only reference consuming
      rl += n;
  }
  return rl;
}

void reverse_qual(uint8_t *q, int l){
  for (int i = 0; i < l/2; ++i){
    q[i]^=q[l-i-1];
    q[l-i-1]^=q[i];
    q[i]^=q[l-i-1];
  }
}

void reverse_cigar(uint32_t *cigar, int l){
  for (int i = 0; i < l/2; ++i){
    cigar[i]^=cigar[l-i-1];
    cigar[l-i-1]^=cigar[i];
    cigar[i]^=cigar[l-i-1];
  }
}

double mean_quality(uint8_t *a, int s, int e){
  double m = 0;
  for (int i = s; i < e; ++i){
    m += (int)a[i];
  }
  m = m/(e-s);
  return m;
}

cigar_ quality_trim(bam1_t* r, int qual_threshold = 20, int sliding_window = 4){
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  uint8_t *qual = bam_get_qual(r);
  int32_t start_pos;
  if(bam_is_rev(r)){
    reverse_qual(qual, r->core.l_qseq);
  }
  double m = 60;
  int del_len, cig, temp;
  uint32_t i = 0, j = 0;
  while(m >= qual_threshold && (i < r->core.l_qseq - sliding_window)){
    m = mean_quality(qual, i, i+sliding_window);
    i++;
  }
  del_len = r->core.l_qseq - i;
  start_pos = get_pos_on_reference(cigar, r->core.n_cigar, del_len, r->core.pos); // For reverse reads need to set core->pos
  if(bam_is_rev(r) && start_pos <= r->core.pos) {
    return {
      cigar,
	r->core.n_cigar,
	r->core.pos
	};
  }
  int32_t n;
  i = 0;
  if(bam_is_rev(r)){
    reverse_cigar(cigar, r->core.n_cigar);
  }
  reverse_cigar(cigar, r->core.n_cigar); // Reverse cigar and trim the beginning of read.
  while(i < r->core.n_cigar){
    if (del_len == 0){
      ncigar[j] = cigar[i];
      i++;
      j++;
      continue;
    }
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if ((bam_cigar_type(cig) & 1)){ // Consumes Query
      if(del_len >= n ){
	ncigar[j] = bam_cigar_gen(n, BAM_CSOFT_CLIP);
      } else if (del_len < n){
	ncigar[j] = bam_cigar_gen(del_len, BAM_CSOFT_CLIP);
      }
      j++;
      temp = n;
      n = std::max(n - del_len, 0);
      del_len = std::max(del_len - temp, 0);
      if(n > 0){
	ncigar[j] = bam_cigar_gen(n, cig);
	j++;
      }
    }
    i++;
  }
  reverse_cigar(ncigar, j);	// Reverse Back
  if(bam_is_rev(r)){
    reverse_cigar(ncigar, j);
  }
  return {
    ncigar,
      j,
      start_pos
  };
}

void print_cigar(uint32_t *cigar, int nlength){
  for (int i = 0; i < nlength; ++i){
    std::cout << ((cigar[i]) & BAM_CIGAR_MASK);
    std::cout << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
}

cigar_ primer_trim(bam1_t *r, int32_t new_pos){
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  uint32_t i = 0, j = 0;
  int del_len, cig, temp;
  if (bam_is_rev(r)){
    del_len = bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r)) - get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos);
    // print_cigar(cigar, r->core.n_cigar);
    reverse_cigar(cigar, r->core.n_cigar);
    // print_cigar(cigar, r->core.n_cigar);
    // std::cout << "Del Length " << del_len << std::endl;
  } else {
    del_len = get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos);
  }
  int32_t n;
  while(i < r->core.n_cigar){
    if (del_len == 0){
      ncigar[j] = cigar[i];
      i++;
      j++;
      continue;
    }
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if ((bam_cigar_type(cig) & 1)){ // Consumes Query
      if(del_len >= n ){
	ncigar[j] = bam_cigar_gen(n, BAM_CSOFT_CLIP);
      } else if (del_len < n){
	ncigar[j] = bam_cigar_gen(del_len, BAM_CSOFT_CLIP);
      }
      j++;
      temp = n;
      n = std::max(n - del_len, 0);
      del_len = std::max(del_len - temp, 0);
      if(n > 0){
	ncigar[j] = bam_cigar_gen(n, cig);
	j++;
      }
    }
    i++;
  }
  if(bam_is_rev(r)){
    reverse_cigar(ncigar, j);
  }
  cigar_ t = {
    ncigar,
    j
  };
  return t;
}

static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
  if (n != b->core.n_cigar) {
    int o = b->core.l_qname + b->core.n_cigar * 4;
    if (b->l_data + (n - b->core.n_cigar) * 4 > b->m_data) {
      b->m_data = b->l_data + (n - b->core.n_cigar) * 4;
      kroundup32(b->m_data);
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->l_data - o);
    memcpy(b->data + b->core.l_qname, cigar, n * 4);
    b->l_data += (n - b->core.n_cigar) * 4;
    b->core.n_cigar = n;
  } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

int16_t get_overlapping_primer_indice(bam1_t* r, std::vector<primer> primers){
  uint32_t query_pos, start_pos, *cigar = bam_get_cigar(r);
  if(bam_is_rev(r)){
    start_pos = bam_endpos(r);
    query_pos = start_pos + (bam_cigar2qlen(r->core.n_cigar, cigar) - get_pos_on_query(cigar, r->core.n_cigar, start_pos, r->core.pos));
  } else {
    start_pos = r->core.pos;
    query_pos = start_pos - get_pos_on_query(cigar, r->core.n_cigar, start_pos, r->core.pos);
  }
  for(std::vector<int>::size_type i = 0; i!=primers.size();i++){
    if(query_pos >= primers[i].get_start() && query_pos <= primers[i].get_end() && start_pos >= primers[i].get_start() && start_pos <= primers[i].get_end()) // Change int to int32_t in primer_bed.cpp
      return i;
  }
  return -1;
}

cigar_ condense_cigar(uint32_t* cigar, uint32_t n){
  int len = 0, cig, next_cig, i = 0;
  while(i< n -1){
    cig = bam_cigar_op(cigar[i]);
    next_cig = bam_cigar_op(cigar[i+1]);
    if(cig == next_cig){
      len = bam_cigar_oplen(cigar[i])+bam_cigar_oplen(cigar[i+1]);
      cigar[i] = bam_cigar_gen(len, bam_cigar_op(cigar[i]));
      for(int j = i+1; j < n - 1; j++){
	cigar[j] = cigar[j+1];
      }
      n--;
    } else {
      i++;
    }
  }
  return {
    cigar,
      n,
      0
      };
}

int main(int argc, char* argv[]) {
  std::cout << "Path " << argv[1] <<std::endl;
  std::string bam = std::string(argv[1]);
  std::string region_;
  std::string bed = std::string(argv[2]);
  std::string bam_out = std::string(argv[3]);
  std::vector<primer> primers = populate_from_file(bed);
  if(argc > 4) {
    region_ = std::string(argv[4]);
  }
  if(!bam.empty()) {
    //open BAM for reading
    samFile *in = hts_open(bam.c_str(), "r");
    BGZF *out = bgzf_open(bam_out.c_str(), "w");
    if(in == NULL) {
      throw std::runtime_error("Unable to open BAM/SAM file.");
    }
    //Load the index
    hts_idx_t *idx = sam_index_load(in, bam.c_str());
    if(idx == NULL) {
      throw std::runtime_error("Unable to open BAM/SAM index."); // Generate index
    }
    //Get the header
    bam_hdr_t *header = sam_hdr_read(in);
    bam_hdr_write(out, header);
    if(header == NULL) {
      sam_close(in);
      throw std::runtime_error("Unable to open BAM header.");
    }
    if (region_.empty()){
      std::cout << "Number of references: " << header->n_targets << std::endl;
      for (int i = 0; i < header->n_targets; ++i){
	std::cout << "Reference Name: " << header->target_name[i] << std::endl;
	std::cout << "Reference Length: " << header->target_len[i] << std::endl;
	if(i==0){
	  region_.assign(header->target_name[i]);
	}
      }
      std::cout << "Using Region: " << region_ << std::endl;
    }
    std::string temp(header->text);
    std::string sortFlag ("SO:coordinate");
    if (temp.find(sortFlag)) {
      std::cout << "Sorted By Coordinate" << std::endl; // Sort by coordinate
    } else {
      std::cout << "Not sorted" << std::endl;
    }
    //Initialize iterator
    hts_itr_t *iter = NULL;
    //Move the iterator to the region we are interested in
    iter  = sam_itr_querys(idx, header, region_.c_str());
    if(header == NULL || iter == NULL) {
      sam_close(in);
      throw std::runtime_error("Unable to iterate to region within BAM.");
    }
    //Initiate the alignment record
    bam1_t *aln = bam_init1();
    int ctr = 0;
    cigar_ t;
    int16_t *p = (int16_t*)malloc(sizeof(int16_t));
    while(sam_itr_next(in, iter, aln) >= 0) {
      // std::cout << "Query Length: " << bam_cigar2qlen(aln->core.n_cigar, cigar) << std::endl;
      // std::cout << "Read Length: " << bam_cigar2rlen(aln->core.n_cigar, cigar) << std::endl;
      // std::cout << "Query Start: " << get_query_start(aln) << std::endl;
      // std::cout << "Query End: " << get_query_end(aln) << std::endl;
      *p = get_overlapping_primer_indice(aln, primers);
      if(*p == -1)
	continue;
      if(bam_is_rev(aln)){
	t = primer_trim(aln, primers[*p].get_start() - 1);
      } else {
	t = primer_trim(aln, primers[*p].get_end() + 1);
	aln->core.pos = primers[*p].get_end() + 1;
      }
      // std::cout << *p << std::endl;
      replace_cigar(aln, t.nlength, t.cigar);
      // Quality Trimming
      t = quality_trim(aln);
      if(bam_is_rev(aln))
	aln->core.pos = t.start_pos;
      // std::cout << "Old: " << t.nlength << std::endl;
      // print_cigar(t.cigar, t.nlength);
      // std::cout << std::endl;
      t = condense_cigar(t.cigar, t.nlength);
      // std::cout << "New: " << t.nlength << std::endl;
      // print_cigar(t.cigar, t.nlength);
      // std::cout << std::endl;
      replace_cigar(aln, t.nlength, t.cigar);
      // std::cout << "Name: " << name  << std::endl;
      // std::cout << "Seq: " << seq << "\tQual: " << qual;
      // std::cout << std::endl;
      // if (strcmp(bam_get_qname(aln), "M01244:143:000000000-BHWC7:1:1104:13925:8758") == 0){
      // 	std::cout << "Start Pos: " << aln->core.pos << std::endl;
      // 	uint8_t *q = bam_get_qual(aln);
      // 	for (int i = 0; i < aln->core.l_qseq -4; ++i){
      // 	  std::cout << mean_quality(q, i, i+4) << " ";
      // 	}
      // 	std::cout << std::endl;
      // }
      int min_length = 30;
      if(bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)) >= min_length){
	// bam1_t *b, const char tag[2], char type, int len, const uint8_t *data
	bam_aux_append(aln, "xa", 'i', 4, (uint8_t*) p);
	bam_write1(out, aln);
      }
      ctr++;
      if(ctr % 100000 == 0){
	std::cout << ctr << std::endl;
      }
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
    bgzf_close(out);
  }
  return 0;
}

