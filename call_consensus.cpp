#include "htslib/vcf.h"

#include<iostream>
#include<fstream>

struct ad {
  int8_t *v8;
  int16_t *v16;
  int32_t *v32;
  float *vf;
  int type;
  int8_t maj_v8;
  int16_t maj_v16;
  int32_t maj_v32;
  float maj_vf;
};

int get_majority(bcf_fmt_t fmt, ad val){
  int maj_ind = -1;
  val.maj_v8 = 1;
  val.maj_v16 = 1;
  val.maj_v32 = 1;
  val.maj_vf = 1;
  for (int j = 0; j < fmt.n; ++j){
    // std::cout << p[j] << ":";
    switch(val.type){
    case 1:
      if(val.v8[j] >= val.maj_v8){
	maj_ind = j;
	val.maj_v8 = val.v8[j];
      }
      break;
    case 2:
      if(val.v16[j] >= val.maj_v16){
	maj_ind = j;
	val.maj_v16 = val.v16[j];
      }
      break;
    case 3:
      if(val.v32[j] >= val.maj_v32){
	maj_ind = j;
	val.maj_v32 = val.v32[j];
      }
      break;
    case 4:
      if(val.vf[j] >= val.maj_vf){
	maj_ind = j;
	val.maj_vf = val.vf[j];
      }
      break;
    }
  }
  return maj_ind;
}

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
    int prev_pos = 0;
    int s;
    while(bcf_read(in, header, v) ==0){
      // std::cout << "Position: " << v->pos << std::endl;
      // std::cout << "Number of Alleles: " << v->n_allele << std::endl;
      for (int i = prev_pos; i < v->pos - 1; ++i){
	out << "?";
      }
      bcf_unpack(v, BCF_UN_FMT);
      bcf_fmt_t *fmt = v->d.fmt;
      for (int i = 0; i < v->n_fmt; ++i){
	// std::cout << "ID: " << fmt[i].id << std::endl;
	// std::cout << " Type: " << fmt[i].type << std::endl;
	// std::cout <<" Size: " << fmt[i].size << std::endl;
	// std::cout << "String: " << header->id[BCF_DT_ID][fmt[i].id].key << std::endl;
	if(strcmp(header->id[BCF_DT_ID][fmt[i].id].key, "AD") == 0){
	  int maj_ind = -1;
	  ad tmp;
	  tmp.type = fmt[i].type;
	  switch(fmt[i].type){
	  case BCF_BT_INT8: tmp.v8 = (int8_t*) fmt[i].p; break;
	  case BCF_BT_INT16: tmp.v16 = (int16_t*) fmt[i].p; break;
	  case BCF_BT_INT32: tmp.v32 = (int32_t*) fmt[i].p; break;
	  case BCF_BT_FLOAT: tmp.vf = (float*) fmt[i].p; break;
	  default: hts_log_error("Unexpected type %d", fmt[i].type); exit(1);
	  }
	  maj_ind = get_majority(fmt[i], tmp);
	  if(maj_ind == -1){
	    out << "?";
	  } else {
	    bcf_unpack(v, BCF_UN_STR);
	    int d = strlen(v->d.allele[0]) - strlen(v->d.allele[maj_ind]); // Ref - Allele
	    std::cout << "Length " << d << std::endl;
	    if(d == 0){
	      if(strlen(v->d.allele[maj_ind]) == 1){
		out << v->d.allele[maj_ind];
	      }		       // Else there is no indel as it matches reference
	    } else if(d > 0) {	// Deletion
	      s = strlen(v->d.allele[0]);
	      for (int j = d; j >0; --j){
		out << v->d.allele[0][s-d];
	      }
	    } else if(d < 0) {	// Insertion
	      s = strlen(v->d.allele[maj_ind]);
	      d = -1 * d;
	      for (int j = d; j > 0; --j){
		out << v->d.allele[maj_ind][s-d];
	      }
	    }
	  }
	  // std::cout << "Indice: " << maj_ind << std::endl;
	  // std::cout << v->d.als[maj_ind] << std::endl;

	  // std::cout << std::endl << std::endl;
	}
      }
      prev_pos = v->pos;
      ctr++;
      if(ctr % 1000 == 0)
	std::cout << ctr << std::endl;
    }
    bcf_destroy(v);
    bcf_hdr_destroy(header);
    bcf_close(in);
    out.close();
  }
  return 0;
}
