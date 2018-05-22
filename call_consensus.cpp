#include "htslib/vcf.h"

#include<iostream>
#include<fstream>

int main(int argc, char* argv[]) {
  std::cout << "Path " << argv[1] <<std::endl;
  std::string bcf = std::string(argv[1]);
  std::string region_;
  std::string out_path = std::string(argv[2]);
  if(argc > 3) {
    region_ = std::string(argv[3]);
  }
  if(!bcf.empty()) {
    //open BCF for reading
    htsFile *in = bcf_open(bcf.c_str(), "r");
    std::ofstream out(out_path);
    out << ">" << out_path << "\n";
    if(in == NULL) {
      throw std::runtime_error("Unable to open BCF/VCF file.");
    }
    int fmt_id;
    //Get the header
    bcf_hdr_t *header = bcf_hdr_read(in);
    if(header == NULL) {
      bcf_close(in);
      throw std::runtime_error("Unable to open BCF header.");
    }
    std::cout << "AD: " << bcf_hdr_id2int(header, BCF_DT_ID, "AD") << std::endl;
    for (int i = 0; i < 3; ++i){
      std::cout << "Key: " << header->id[i]->key << std::endl;
      std::cout << "Value: " <<  header->id[i]->val << std::endl;
    }
    fmt_id = bcf_hdr_id2int(header, BCF_DT_ID, "AD");
    //Initialize iterator
    bcf1_t *v = bcf_init();
    int ctr = 0;
    while(bcf_read(in, header, v) ==0){
      // std::cout << "Position: " << v->pos << std::endl;
      // std::cout << "Number of Alleles: " << v->n_allele << std::endl;
      if(v->n_allele == 1){
	bcf_unpack(v, BCF_UN_STR);
	char *ALT = v->d.als;
	out << ALT[0];
      } else {
	bcf_unpack(v, BCF_UN_FMT);
	bcf_fmt_t *fmt = v->d.fmt;
	bcf_info_t *info = v->d.info;
	int offset = 0;
	for (int i = 0; i < v->n_fmt; ++i){
	  // std::cout << "ID: " << fmt[i].id << std::endl;
	  // std::cout << " Type: " << fmt[i].type << std::endl;
	  // std::cout <<" Size: " << fmt[i].size << std::endl;
	  // std::cout << "String: " << header->id[BCF_DT_ID][fmt[i].id].key << std::endl;
	  if(strcmp(header->id[BCF_DT_ID][fmt[i].id].key, "AD") == 0){
	    int16_t *p= (int16_t *) (fmt[i].p);
	    int maj_ind = -1;
	    int16_t maj_val = 0;
	    for (int j = 0; j < fmt[i].n; ++j){
	      // std::cout << p[j] << ":";
	      if(p[j] > maj_val){
		maj_ind = j;
		maj_val = p[j];
	      }
	    }
	    // std::cout << "Indice: " << maj_ind << std::endl;
	    bcf_unpack(v, BCF_UN_STR);
	    // std::cout << v->d.als[maj_ind] << std::endl;
	    out << v->d.als[maj_ind];
	    // std::cout << std::endl << std::endl;
	  }
	}
      }
      ctr++;
      if(ctr % 1000 == 0)
	std::cout << ctr << std::endl;
    }
    bcf_hdr_destroy(header);
    bcf_destroy(v);
    bcf_close(in);
  }
  return 0;
}
