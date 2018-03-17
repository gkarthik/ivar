#include <iostream>
#include <stdexcept>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "samtools-1.7/samtools.h"

int main(int argc, char* argv[]) {
  printf("Path %s\n",argv[1]);
  std::string bam = std::string(argv[1]);
  std::string region_ = "";
  if(argc > 2) {
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
      throw std::runtime_error("Unable to open BAM/SAM index.");
    }
    //Get the header
    bam_hdr_t *header = sam_hdr_read(in);
    if(header == NULL) {
      sam_close(in);
      throw std::runtime_error("Unable to open BAM header.");
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
      std::cout << "\tPos: " << aln->core.pos;
      std::string seq, qual;
      uint8_t *quali = bam_get_qual(aln);
      uint8_t *seqi = bam_get_seq(aln);
      uint32_t *cigar = bam_get_cigar(aln);
      char *name = bam_get_qname(aln);
      for (int i = 0; i < aln->core.l_qseq; i++) {
	seq += seq_nt16_str[bam_seqi(seqi, i)];
	qual += 33 + quali[i];
      }
      std::cout << "Name: " << cigar << "\t Cigar" << cigar << std::endl;
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
