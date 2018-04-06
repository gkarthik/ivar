#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "samtools-1.7/samtools.h"

#include "primer_bed.h"

struct cigar_ {
  uint32_t *cigar;
  int nlength;
};

std::vector<primer> populate_from_file(std::string path){
  std::ifstream  data(path);
  std::string line;
  std::vector<primer> primers;
  while(std::getline(data,line)){
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	p.set_region(cell);
	break;
      case 1:
	p.set_start(std::stoul(cell));
	break;
      case 2:
	p.set_end(std::stoul(cell));
	break;
      case 3:
	p.set_name(cell);
	break;
      case 4:
	p.set_score(std::stoi(cell));
	break;
      case 5:
	p.set_strand(cell[0]);
      }
      ctr++;
    }
    primers.push_back(p);
  }
  return primers;
}

int32_t get_query_start(bam1_t *r){
  int cig;
  int32_t n;
  uint32_t *cigar = bam_get_cigar(r);
  int32_t qs = 0;
  for (int i = 0; i < r->core.n_cigar; ++i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (cig == BAM_CMATCH){
      return qs;
    }
    if (cig != BAM_CDEL && cig != BAM_CHARD_CLIP){	// In deletion don't increment query_start
      qs += n;
    }
  }
  return qs;
}

int32_t get_query_end(bam1_t *r){
  int cig;
  int32_t n;
  uint32_t *cigar = bam_get_cigar(r);
  int32_t qs = 0;
  for (int i = r->core.n_cigar - 1; i >= 0; --i){
    cig  = bam_cigar_op(cigar[i]);
    n = bam_cigar_oplen(cigar[i]);
    if (cig == BAM_CMATCH){
      break;
    }
    qs += n;
  }
  return bam_cigar2rlen(r->core.n_cigar, cigar) - qs;
}

void reverse_cigar(uint32_t *cigar, int l){
  for (int i = 0; i < l/2; ++i){
    cigar[i]^=cigar[l-i-1];
    cigar[l-i-1]^=cigar[i];
    cigar[i]^=cigar[l-i-1];
  }
}

void print_cigar(uint32_t *cigar, int nlength){
  for (int i = 0; i < nlength; ++i){
    std::cout << ((cigar[i]) & BAM_CIGAR_MASK) << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
}

cigar_ primer_trim(bam1_t *r, int32_t new_pos){
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  int i = 0, j = 0, del_len, cig, _cig, temp;
  if (bam_is_rev(r)){
    del_len = bam_endpos(r) + (bam_cigar2rlen(r->core.n_cigar, bam_get_cigar(r)) - get_query_end(r)) - new_pos;
    // print_cigar(cigar, r->core.n_cigar);
    reverse_cigar(cigar, r->core.n_cigar);
    // print_cigar(cigar, r->core.n_cigar);
    // std::cout << "Del Length " << del_len << std::endl;
  } else {
    del_len = new_pos - (r->core.pos - get_query_start(r));
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
    if (bam_cigar_type(cig) & 2){ // Consumes Reference
      _cig = cig;
      if(cig == BAM_CMATCH){
	_cig = BAM_CSOFT_CLIP;
      }
      if(del_len >= n ){
	ncigar[j] = bam_cigar_gen(n, cig);
      } else if (del_len < n){
	ncigar[j] = bam_cigar_gen(del_len, cig);
      }
      j++;
      temp = n;
      n = std::max(n - del_len, 0);
      del_len = std::max(del_len - temp, 0);
    } else if ((bam_cigar_type(cig) & 1)) { // Consumes only Query
      del_len = std::max(del_len - n, 0);
    }
    if(n > 0){
      ncigar[j] = bam_cigar_gen(n, cig);
      j++;
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

int get_overlapping_primer_indice(bam1_t* r, std::vector<primer> primers){
  int32_t pos = r->core.pos - get_query_start(r);
  if(bam_is_rev(r)){
    pos = bam_endpos(r) + (bam_cigar2rlen(r->core.n_cigar, bam_get_cigar(r)) - get_query_end(r));
  }
  for(std::vector<int>::size_type i = 0; i!=primers.size();i++){
    if(pos >= primers[i].get_start() && pos <= primers[i].get_end()){ // Change int to int32_t in primer_bed.cpp
      return i;
    }
  }
  return -1;
}

int main(int argc, char* argv[]) {
  std::cout << "Path " << argv[1] <<std::endl;
  std::string bam = std::string(argv[1]);
  std::string region_;
  std::string bed = std::string(argv[2]);
  std::vector<primer> primers = populate_from_file(bed);
  if(argc > 3) {
    region_ = std::string(argv[2]);
  }
  if(!bam.empty()) {
    //open BAM for reading
    samFile *in = sam_open(bam.c_str(), "r");
    samFile *out = sam_open("temp.bam", "w");
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
    sam_hdr_write(out, header);
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
    std::string t(header->text);
    std::string sortFlag ("SO:coordinate");
    if (t.find(sortFlag)) {
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
    while(sam_itr_next(in, iter, aln) >= 0) {
      std::string seq, qual;
      uint8_t *quali = bam_get_qual(aln);
      uint8_t *seqi = bam_get_seq(aln);
      uint32_t *cigar = bam_get_cigar(aln);
      char *name = bam_get_qname(aln);
      for (int i = 0; i < aln->core.l_qseq; i++) {
	seq += seq_nt16_str[bam_seqi(seqi, i)];
	qual += 33 + quali[i];
      }
      // std::cout << "Query Length: " << bam_cigar2qlen(aln->core.n_cigar, cigar) << std::endl;
      // std::cout << "Read Length: " << bam_cigar2rlen(aln->core.n_cigar, cigar) << std::endl;
      // std::cout << "Query Start: " << get_query_start(aln) << std::endl;
      // std::cout << "Query End: " << get_query_end(aln) << std::endl;
      int p = get_overlapping_primer_indice(aln, primers);
      if(p != -1){
	cigar_ t;
	if(bam_is_rev(aln)){
	  t = primer_trim(aln, primers[p].get_start() - 1);
	} else {
	  t = primer_trim(aln, primers[p].get_end() + 1);
	}
	replace_cigar(aln, t.nlength, t.cigar);
	// std::cout << std::endl;
      } else {
	// std::cout << "On Primer: FALSE" << std::endl;
      }
      // std::cout << std::endl;
      // std::cout << "Name: " << name  << std::endl;
      // std::cout << "Seq: " << seq << "\tQual: " << qual;
      // std::cout << std::endl;
      sam_write1(out, header, aln);
      ctr++;
      if(ctr % 1000 == 0){
	std::cout << ctr << std::endl;
      }
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);
  }
  return 0;
}
