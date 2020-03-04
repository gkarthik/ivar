#include "ref_seq.h"

char ref_antd::get_base(int64_t pos, std::string region){
  
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
  this->add_seq(ref_path);
  this->add_gff(gff_path);
}

std::ostringstream ref_antd::codon_aa_stream(std::ostringstream &line_stream, std::ofstream &fout, int64_t pos, char alt){
  
}
