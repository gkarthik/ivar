#include "call_variants.h"

const char gap='N';
const float sig_level = 0.01;

std::vector<allele>::iterator get_ref_allele(std::vector<allele> &ad, char ref){
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc[0] == ref)
      return it;
  }
  // print_allele_depths(ad);
  return ad.end();
}

double* get_frequency_depth(allele a, uint32_t pos_depth, uint32_t total_depth){ // pos_depth is ungapped depth after passing quality. total_depth is depth without quality filtering for indels.
  double *val = new double[2];
  if(a.nuc[0] == '+'){	// For insertions use depth discarding quality
    val[0] = a.depth/(double)total_depth;
    val[1] = total_depth;
    return val;
  }
  val[0] = a.depth/(double)pos_depth;
  val[1] = pos_depth;
  return val;
}

// char get_ref_base(char* ref_seq, int64_t pos_offset, int64_t edit_pos, std::string edit_seq){
//   if(edit_pos == -1 || pos_offset < edit_pos){
//     return *(ref_seq + pos_offset); // Before edit_pos or no edit
//   } else if(pos_offset > edit_pos - 1 + edit_seq.size()){ // Curent position is after edit so just change positions
//     return *(ref_seq + pos_offset + edit_seq.size());
//   }
//   // Current position is within the edited region: (pos_offset >= edit_pos - 1 && pos_offset <= edit_pos -1 + edit_seq.size())
//   return edit_seq[pos_offset - (edit_pos - 1)];
// }

// int64_t calculate_codon_start_position(uint64_t feature_start, uint64_t current_pos, int phase, int ref_len, int64_t edit_pos, std::string edit_seq){
//   int64_t start_pos;
//   if(edit_pos == -1 || current_pos < edit_pos)
//     start_pos = (feature_start - 1) + (((current_pos - (feature_start + phase)))/3)*3;
//   else if(current_pos >= edit_pos)
//     start_pos = (feature_start - 1) + (((current_pos - (feature_start + phase)))/3)*3; // Add edit_seq.size() here
//   return start_pos;
// }

// int write_aa(std::ofstream &fout, int64_t start_pos, uint64_t pos, char *ref_seq, char alt, std::string orf_id, int ref_len, int64_t edit_pos, std::string edit_seq){
//   if(start_pos < 0 || start_pos + 2 > ref_len + edit_seq.size()){ // Cases where the current_pos is not part of a codon presnt within reference sequence length plus edit_seq which is 0 by default
//     fout << EMPTY_AA_FIELDS << std::endl;
//     fout << std::endl;
//   }
//   int tmp;
//   char *aa_codon = new char[3], *ref_codon = new char[3];
//   fout << orf_id << "\t";
//   for (tmp = 0 ; tmp < 3; ++tmp) {
//     ref_codon[tmp] = get_ref_base(ref_seq, start_pos + tmp, edit_pos, edit_seq);
//     fout << ref_codon[tmp];
//   }
//   fout << "\t";
//   fout << codon2aa(ref_codon[0], ref_codon[1], ref_codon[2]) << "\t";
//   for (tmp = 0 ; tmp < 3; ++tmp) {
//     if(pos - 1 == start_pos + tmp){
//       fout << alt;
//       aa_codon[tmp] = alt;
//     } else {
//       fout << get_ref_base(ref_seq, start_pos + tmp, edit_pos, edit_seq);
//       aa_codon[tmp] = get_ref_base(ref_seq, start_pos + tmp, edit_pos, edit_seq);
//     }
//   }
//   fout << "\t";
//   fout << codon2aa(aa_codon[0], aa_codon[1], aa_codon[2]);
//   fout << std::endl;
//   return 0;
// }

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual, double min_threshold, uint8_t min_depth, std::string ref_path, std::string gff_path){
  std::string line, cell, bases, qualities, region;
  ref_antd ref_antd(ref_path, gff_path);
  char *ref_codon = new char[3], *alt_codon = new char[3];
  std::ostringstream out_str;
  std::ofstream fout((out_file+".tsv").c_str());
  fout << "REGION"
    "\tPOS"
    "\tREF"
    "\tALT"
    "\tREF_DP"
    "\tREF_RV"
    "\tREF_QUAL"
    "\tALT_DP"
    "\tALT_RV"
    "\tALT_QUAL"
    "\tALT_FREQ"
    "\tTOTAL_DP"
    "\tPVAL"
    "\tPASS"
    "\tGFF_FEATURE"
    "\tREF_CODON"
    "\tREF_AA"
    "\tALT_CODON"
    "\tALT_AA"
       << std::endl;
  int ctr = 0, pos = 0, tmp;
  int64_t start_pos = 0;
  uint32_t mdepth = 0, pdepth = 0; // mpdepth for mpileup depth and pdeth for ungapped depth at position
  double pval_left, pval_right, pval_twotailed, *freq_depth, err;
  std::stringstream line_stream;
  char ref;
  std::vector<allele> ad;
  std::vector<allele>::iterator ref_it;
  while (std::getline(cin, line)){
    line_stream << line;
    ctr = 0;
    while(std::getline(line_stream,cell,'\t')){
      switch(ctr){
      case 0:
	region = cell;
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	// Read from ref if ref_seq is set, else read from mpileup
	ref = ref_antd.get_base(pos, region);
	ref = (ref == 0) ? cell[0] : ref;
	break;
      case 3:
	mdepth = stoi(cell);
	break;
      case 4:
	bases = cell;
	break;
      case 5:
	qualities = cell;
	break;
      case 6:
	break;
      }
      ctr++;
    }
    if(mdepth < min_depth) {	// Check for minimum depth
      line_stream.clear();
      continue;
    }
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    if(ad.size() == 0){
      line_stream.clear();
      continue;
    }
    ref_it = get_ref_allele(ad, ref);
    if (ref_it == ad.end()) {	// If ref not present in reads.
      allele a;
      a.nuc = ref;
      a.depth = 0;
      a.reverse = 0;
      a.mean_qual = 0;
      ad.push_back(a);
      ref_it = ad.end() - 1;
    }
    // Get ungapped coverage
    pdepth = 0;
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if(it->nuc[0]=='*' || it->nuc[0] == '+' || it->nuc[0] == '-')
	continue;
      pdepth += it->depth;
    }
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if((*it == *ref_it) || it->nuc[0]=='*')
	continue;
      freq_depth = get_frequency_depth(*it, pdepth, mdepth);
      if(freq_depth[0] < min_threshold)
	continue;
      out_str << region << "\t";
      out_str << pos << "\t";
      out_str << ref << "\t";
      out_str << it->nuc << "\t";
      out_str << ref_it->depth << "\t";
      out_str << ref_it->reverse << "\t";
      out_str << (uint16_t)ref_it->mean_qual << "\t";
      out_str << it->depth << "\t";
      out_str << it->reverse << "\t";
      out_str << (uint16_t) it->mean_qual << "\t";
      out_str << freq_depth[0] << "\t";
      out_str << freq_depth[1] << "\t";
      /*
	    | Var   | Ref      |
	Exp | Error | Err free |
	Obs | AD    | RD       |
       */
      err = pow(10, ( -1 * (it->mean_qual)/10));
      kt_fisher_exact((err * mdepth), (1-err) * mdepth, it->depth, ref_it->depth, &pval_left, &pval_right, &pval_twotailed);
      out_str << pval_left << "\t";
      if(pval_left <= sig_level){
	out_str << "TRUE" << "\t";
      } else {
	out_str << "FALSE" << "\t";
      }
      if(it->nuc[0] != '+'){
	ref_antd.codon_aa_stream(out_str, fout, pos, it->nuc[0]);
      } else {
	fout << "\t\t\t\t\t";
      }
    }
    line_stream.clear();
  }
  fout.close();
  return 0;
}
