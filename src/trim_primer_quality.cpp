#include "trim_primer_quality.h"

#define round_int(x,total) ((int) (0.5 + ((float)x / float(total)) * 10000))/(float)100

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

// Number of bases from 3' end for reverse reads
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
    m += (double)a[i];
  }
  m = m/(e-s);
  return m;
}

cigar_ quality_trim(bam1_t* r, uint8_t qual_threshold, uint8_t sliding_window){
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
  cigar_ t;
  while(m >= qual_threshold && (i < (uint32_t)r->core.l_qseq)){
    m = mean_quality(qual, i, i+sliding_window);
    i++;
    if(i > (uint32_t)r->core.l_qseq - sliding_window)
      sliding_window--;
  }
  // Reverse qual back.
  if(bam_is_rev(r)){
    reverse_qual(qual, r->core.l_qseq);
  }
  del_len = r->core.l_qseq - i;
  start_pos = get_pos_on_reference(cigar, r->core.n_cigar, del_len, r->core.pos); // For reverse reads need to set core->pos.
  if(bam_is_rev(r) && start_pos <= r->core.pos) {
    t.cigar = cigar;
    t.nlength = r->core.n_cigar;
    t.start_pos = r->core.pos;
    return t;
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
  t.cigar = ncigar;
  t.nlength = j;
  t.start_pos = start_pos;
  return t;
}

void print_cigar(uint32_t *cigar, int nlength){
  for (int i = 0; i < nlength; ++i){
    std::cout << ((cigar[i]) & BAM_CIGAR_MASK);
    std::cout << "-" << ((cigar[i]) >> BAM_CIGAR_SHIFT) << " ";
  }
  std::cout << std::endl;
}

cigar_ primer_trim(bam1_t *r, int32_t new_pos){
  uint32_t *ncigar = (uint32_t*) malloc(sizeof(uint32_t) * (r->core.n_cigar + 1)), // Maximum edit is one more element with soft mask
    *cigar = bam_get_cigar(r);
  uint32_t i = 0, j = 0;
  int del_len, cig, temp;
  if (bam_is_rev(r)){
    del_len = bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r)) - get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos) - 1;
    reverse_cigar(cigar, r->core.n_cigar);
  } else {
    del_len = get_pos_on_query(cigar, r->core.n_cigar, new_pos, r->core.pos);
  }
  int32_t n;
  while(i < r->core.n_cigar){
    if (del_len == 0){ // del_len ==0 or consumes only reference
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
    } else { // Only consumes reference: deletion or consumes nothing
      ncigar[j] = cigar[i];
      j++;
    }
    i++;
  }
  if(bam_is_rev(r)){
    reverse_cigar(ncigar, j);
  }
  cigar_ t = {
    ncigar,
    j,
    0
  };
  return t;
}

void replace_cigar(bam1_t *b, uint32_t n, uint32_t *cigar){
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

void get_overlapping_primers(bam1_t* r, std::vector<primer> primers, std::vector<primer> &overlapped_primers){
  overlapped_primers.clear();
  uint32_t start_pos;
  if(bam_is_rev(r)){
    start_pos = bam_endpos(r)-1;
  } else {
    start_pos = r->core.pos;
  }
  for(std::vector<primer>::iterator it = primers.begin(); it != primers.end(); ++it) {
    if(start_pos >= it->get_start() && start_pos <= it->get_end())
      overlapped_primers.push_back(*it);
  }
}

cigar_ remove_trailing_query_ref_consumption(uint32_t* cigar, int32_t n){
  int i = 0, len = 0, cig, start_pos = 0;
  cigar_ t;
  while(i < n){
    cig = bam_cigar_op(cigar[i]);
    len = bam_cigar_oplen(cigar[i]);
    if((bam_cigar_type(cig) & 2) && (bam_cigar_type(cig) & 1))
      break;
    if(bam_cigar_type(cig) & 1){	// Consumes only query - insertion
      cigar[i] = bam_cigar_gen(len, BAM_CSOFT_CLIP);
    } else {	// Consumes only reference - deletion or consumes nothing like padding
      for (int j = i; j < n-1; ++j){
	cigar[j] = cigar[j+1];
      }
      n--;
      i--;
      start_pos += len;
    }
    i++;
  }
  reverse_cigar(cigar, n);
  i = 0, len = 0;
  while(i < n){			// TODO: This repeats loop from above. Restructure!
    cig = bam_cigar_op(cigar[i]);
    len = bam_cigar_oplen(cigar[i]);
    if((bam_cigar_type(cig) & 2) && (bam_cigar_type(cig) & 1))
      break;
    if(bam_cigar_type(cig) & 1){	// Consumes only query - insertion
      cigar[i] = bam_cigar_gen(len, BAM_CSOFT_CLIP);
    } else {	// Consumes only reference - deletion or consumes nothing like padding
      for (int j = i; j < n-1; ++j){
	cigar[j] = cigar[j+1];
      }
      n--;
      i--;
      start_pos += len;
    }
    i++;
  }
  reverse_cigar(cigar, n);
  t.cigar = cigar;
  t.nlength = n;
  t.start_pos = start_pos;	// Here start_pos is number of softclipped bases that don't consume query
  return t;
}

cigar_ condense_cigar(uint32_t* cigar, uint32_t n){
  uint32_t i = 0, len = 0, cig, next_cig, start_pos = 0;
  cigar_ t = remove_trailing_query_ref_consumption(cigar, n);
  cigar = t.cigar;
  n = t.nlength;
  start_pos = t.start_pos;
  while(i< n -1){
    cig = bam_cigar_op(cigar[i]);
    next_cig = bam_cigar_op(cigar[i+1]);
    if(cig == next_cig){
      len = bam_cigar_oplen(cigar[i])+bam_cigar_oplen(cigar[i+1]);
      cigar[i] = bam_cigar_gen(len, bam_cigar_op(cigar[i]));
      for(uint32_t j = i+1; j < n - 1; j++){
	cigar[j] = cigar[j+1];
      }
      n--;
    } else {
      i++;
    }
  }
  t.cigar = cigar;
  t.nlength = n;
  t.start_pos = start_pos;
  return t;
}

void add_pg_line_to_header(bam_hdr_t** hdr, char *cmd){
  size_t len = strlen((*hdr)->text) + strlen(cmd)+1;
  char * new_text = (char *)malloc(len);
  memcpy(new_text, (*hdr)->text, strlen((*hdr)->text));
  new_text[strlen((*hdr)->text)] = '\0';
  strcat(new_text, cmd);
  free((*hdr)->text);
  (*hdr)->text = new_text;
  new_text = NULL;
  (*hdr)->l_text = len-1;
}

int trim_bam_qual_primer(std::string bam, std::string bed, std::string bam_out, std::string region_, uint8_t min_qual, uint8_t sliding_window, std::string cmd, bool write_no_primer_reads, int min_length = 30){
  std::vector<primer> primers;
  if(!bed.empty()){
    primers = populate_from_file(bed);
  }
  // if(primers.size() == 0){
  //   return 0;
  // }
  if(bam.empty()){
    std::cout << "Bam file in empty." << std::endl;
    return -1;
  }
  bam_out += ".bam";
  samFile *in = hts_open(bam.c_str(), "r");
  BGZF *out = bgzf_open(bam_out.c_str(), "w");
  if(in == NULL) {
    std::cout << ("Unable to open BAM/SAM file.") << std::endl;
    return -1;
  }
  //Load the index
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  if(idx == NULL) {
    if(sam_index_build2(bam.c_str(), 0, 0)< 0){
      std::cout << ("Unable to open BAM/SAM index.") << std::endl;
      return -1;
    } else {
      idx = sam_index_load(in, bam.c_str());
    }
  }
  //Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  if(header == NULL) {
    sam_close(in);
    std::cout << "Unable to open BAM/SAM header." << std::endl;
  }
  add_pg_line_to_header(&header, const_cast<char *>(cmd.c_str()));
  if(bam_hdr_write(out, header) < 0){
    std::cout << "Unable to write BAM header to path." << std::endl;
    sam_close(in);
    return -1;
  }
  // Get relevant region
  int region_id;
  uint64_t unmapped, mapped, log_skip;
  std::cout << std::endl << "Number of references in file: " << header->n_targets << std::endl;
  for (int i = 0; i < header->n_targets; ++i){
    std::cout << header->target_name[i] << std::endl;
    if(region_.compare(std::string(header->target_name[i])) == 0){
      region_id = i;
    }
    if(i==0){			// Reading only first reference
      region_.assign(header->target_name[i]);
      region_id = i;
    }
  }
  std::cout << "Using Region: " << region_ << std::endl << std::endl;
  // Get index stats
  hts_idx_get_stat(idx, region_id, &mapped, &unmapped);
  std::cout << "Found " << mapped << " mapped reads" << std::endl;
  std::cout << "Found " << unmapped << " unmapped reads" << std::endl;
  std::string hdr_text(header->text);
  if (hdr_text.find(std::string("SO:coordinate")) != std::string::npos) {
    std::cout << "Sorted By Coordinate" << std::endl; // Sort by coordinate
  } else if(hdr_text.find(std::string("SO:queryname")) != std::string::npos) {
    std::cout << "Sorted By Query Name" << std::endl; // Sort by name
  } else {
    std::cout << "Not sorted" << std::endl;
  }
  std::cout << "-------" << std::endl;
  log_skip = (mapped + unmapped > 10) ? (mapped + unmapped)/10 : 2;
  //Initialize iterator
  hts_itr_t *iter = NULL;
  //Move the iterator to the region we are interested in
  iter  = sam_itr_querys(idx, header, region_.c_str());
  if(header == NULL || iter == NULL) {
    sam_close(in);
    std::cout << "Unable to iterate to region within BAM/SAM." << std::endl;
    return -1;
  }
  //Initiate the alignment record
  bam1_t *aln = bam_init1();
  int ctr = 0;
  cigar_ t;
  uint32_t primer_trim_count = 0, no_primer_counter = 0, low_quality = 0;
  bool unmapped_flag = false;
  uint32_t unmapped_counter = 0;
  primer cand_primer;
  std::vector<primer> overlapping_primers;
  std::vector<primer>::iterator cit;
  while(sam_itr_next(in, iter, aln) >= 0) {
    unmapped_flag = false;
    if((aln->core.flag&BAM_FUNMAP) == 0){
      get_overlapping_primers(aln, primers, overlapping_primers);
      if(overlapping_primers.size() > 0){
	primer_trim_count++;
	if(bam_is_rev(aln)){
	  cand_primer = get_min_start(overlapping_primers);
	  t = primer_trim(aln, cand_primer.get_start() - 1);
	} else {
	  cand_primer = get_max_end(overlapping_primers);
	  t = primer_trim(aln, cand_primer.get_end() + 1);
	  aln->core.pos = cand_primer.get_end() + 1;
	}
	replace_cigar(aln, t.nlength, t.cigar);
	// Add count to primer
	cit = std::find(primers.begin(), primers.end(), cand_primer);
	if(cit != primers.end())
	  cit->add_read_count(1);
      }
      t = quality_trim(aln, min_qual, sliding_window);	// Quality Trimming
      if(bam_is_rev(aln))
	aln->core.pos = t.start_pos;
      t = condense_cigar(t.cigar, t.nlength);
      // aln->core.pos += t.start_pos;
      replace_cigar(aln, t.nlength, t.cigar);
    } else {
      unmapped_flag = true;
      unmapped_counter++;
      continue;
    }
    if(bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)) >= min_length){
      if(overlapping_primers.size() > 0){	// Write to BAM only if primer found.
	int16_t cand_ind = cand_primer.get_indice();
	bam_aux_append(aln, "XA", 's', sizeof(cand_ind), (uint8_t*) &cand_ind);
	if(bam_write1(out, aln) < 0){
	  std::cout << "Not able to write to BAM" << std::endl;
	  hts_itr_destroy(iter);
	  hts_idx_destroy(idx);
	  bam_destroy1(aln);
	  bam_hdr_destroy(header);
	  sam_close(in);
	  bgzf_close(out);
	  return -1;
	}
      } else {
	if((primers.size() == 0 || write_no_primer_reads) && !unmapped_flag){ // Write mapped reads to BAM if -e flag given
	  if(bam_write1(out, aln) < 0){
	    std::cout << "Not able to write to BAM" << std::endl;
	    hts_itr_destroy(iter);
	    hts_idx_destroy(idx);
	    bam_destroy1(aln);
	    bam_hdr_destroy(header);
	    sam_close(in);
	    bgzf_close(out);
	    return -1;
	  }
	}
	no_primer_counter++;
      }
    } else {
      low_quality++;
    }
    ctr++;
    if(ctr % log_skip == 0){
      std::cout << "Processed " << (ctr/log_skip) * 10 << "% reads ... " << std::endl;
    }
  }
  std::cout << std::endl << "-------" << std::endl;
  std::cout << "Results: " << std::endl;
  std::cout << "Primer Name" << "\t" << "Read Count" << std::endl;
  for(cit = primers.begin(); cit != primers.end(); ++cit) {
    std::cout << cit->get_name() << "\t" << cit->get_read_count() << std::endl;
  }
  std::cout << std::endl << "Trimmed primers from " << round_int(primer_trim_count, mapped) << "% (" << primer_trim_count <<  ") of reads." << std::endl;
  std::cout << round_int( low_quality, mapped) << "% (" << low_quality << ") of reads were quality trimmed below the minimum length of " << min_length << " bp and were not writen to file." << std::endl;
  if(write_no_primer_reads){
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter << ") of reads started outside of primer regions. Since the -e flag was given, these reads were written to file." << std::endl;
  } else if (primers.size() == 0) {
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter << ") of reads started outside of primer regions. Since there were no primers found in BED file, these reads were written to file." << std::endl;
  } else {
    std::cout << round_int(no_primer_counter, mapped) << "% ("  << no_primer_counter << ") of reads that started outside of primer regions were not written to file." << std::endl;
  }
  if(unmapped_counter > 0){
    std::cout << unmapped_counter << " unmapped reads were not written to file." << std::endl;
  }
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);
  bgzf_close(out);
  return 0;
}
