#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "samtools-1.7/samtools.h"

#include "primer_bed.h"

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
    cig  = (cigar[i]) & BAM_CIGAR_MASK;
    n = (cigar[i]) >> BAM_CIGAR_SHIFT;
    if (cig == BAM_CMATCH){
      return qs;
    }
    qs += n;
  }
  return qs;
}


int32_t get_query_end(bam1_t *r){
  int cig;
  int32_t n;
  uint32_t *cigar = bam_get_cigar(r);
  int32_t qs = 0;
  for (int i = r->core.n_cigar - 1; i >= 0; --i){
    cig  = (cigar[i]) & BAM_CIGAR_MASK;
    n = (cigar[i]) >> BAM_CIGAR_SHIFT;
    if (cig == BAM_CMATCH){
      return qs;
    }
    qs += n;
  }
  return bam_cigar2qlen(r->core.n_cigar, cigar) - qs;
}

int get_overlapping_primer_indice(bam1_t* r, std::vector<primer> primers){
  int32_t pos = r->core.pos - get_query_start(r);
  if(bam_is_rev(r)){
    pos = bam_endpos(r) + (bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r)) - get_query_end(r)) ;
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
      std::cout << "Read Chr: " << header->target_name[aln->core.tid];
      std::cout << "\tPos: " << aln->core.pos << std::endl;
      std::string seq, qual;
      uint8_t *quali = bam_get_qual(aln);
      uint8_t *seqi = bam_get_seq(aln);
      uint32_t *cigar = bam_get_cigar(aln);
      char *name = bam_get_qname(aln);
      for (int i = 0; i < aln->core.l_qseq; i++) {
	seq += seq_nt16_str[bam_seqi(seqi, i)];
	qual += 33 + quali[i];
      }
      std::cout << "Cigar: ";
      int cig, ncig;
      for (int i = 0; i < aln->core.n_cigar; ++i){
        cig  = (cigar[i]) & BAM_CIGAR_MASK;
        ncig = (cigar[i]) >> BAM_CIGAR_SHIFT;
	std::cout << cig << "-" << ncig << " ";
      }
      std::cout << "Query Length: " << bam_cigar2qlen(aln->core.n_cigar, cigar) << std::endl;
      std::cout << "Read Length: " << bam_cigar2rlen(aln->core.n_cigar, cigar) << std::endl;
      std::cout << "Query Start: " << get_query_start(aln) << std::endl;
      std::cout << "Query End: " << get_query_end(aln) << std::endl;
      int p = get_overlapping_primer_indice(aln, primers);
      if(p != -1){
	std::cout << "On Primer: " << primers[p].get_name()  << std::endl;
      } else {
	std::cout << "On Primer: FALSE" << std::endl;
      }
      std::cout << std::endl;
      std::cout << "Name: " << name  << std::endl;
      std::cout << "Seq: " << seq << "\tQual: " << qual;
      std::cout << std::endl;
      if (ctr > 10){
	break;
      }
      ctr++;
    }
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
  }
  return 0;
}
