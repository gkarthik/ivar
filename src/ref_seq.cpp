#include "ref_seq.h"

char ref_antd::get_base(int64_t pos, std::string region){
  if(seq == nullptr)
    return UNKNOWN_BASE;
  return *(seq + (pos - 1));
}

char* ref_antd::get_codon(int64_t pos, std::string region, gff3_feature feature){
  
}

char* ref_antd::get_codon(int64_t pos, std::string region, gff3_feature feature, char alt){
  
}

int ref_antd::add_gff(std::string path){
  // Read GFF file
  if(!path.empty())
    gff.read_file(path);
  return 0;
}

int ref_antd::add_seq(std::string path){
  // Read reference file
  if(!path.empty())
    fai = fai_load(path.c_str());
  if(!fai && !path.empty()){
    std::cout << "Reference file does not exist at " << path << std::endl;
    return -1;
  }
  return 0;
}

ref_antd::ref_antd(std::string ref_path, std::string gff_path){
  this->seq = nullptr;
  this->add_seq(ref_path);
  this->add_gff(gff_path);
}

int ref_antd::codon_aa_stream(std::string region, std::ostringstream &line_stream, std::ofstream &fout, int64_t pos, char alt){
  std::vector<gff3_feature> features = gff.query_features(pos, "CDS");
  std::vector<gff3_feature>::iterator it;
  char *ref_codon = new char[3], *alt_codon = new char[3];
  for(it = features.begin(); it != features.end(); it++){
    fout << line_stream.str() << "\t";
    fout << it->get_attribute("ID") << "\t";
    ref_codon = this->get_codon(pos, region, *it);
    fout << ref_codon << "\t";
    fout << codon2aa(ref_codon[0], ref_codon[1], ref_codon[2]) << "\t";
    alt_codon = this->get_codon(pos, region, *it, alt);
    fout << alt_codon << "\t";
    fout << codon2aa(alt_codon[0], alt_codon[1], alt_codon[2]) << "\t";
    fout << std::endl;
  }
  line_stream.str("");
  line_stream.clear();
  return 0;
}
