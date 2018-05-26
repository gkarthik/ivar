#include "htslib/vcf.h"

#include<iostream>
#include<fstream>

struct maj {
  int *ind;
  int n;
};

#define GET_MAX_INDICES(type_t, pu, no, m) {		\
    type_t v = 1;					\
    type_t *p = (type_t*) (pu);			\
    for (int k = 0; k < no; ++k){			\
      if(p[k] > v){					\
	v = p[k];					\
	m.n = 0;					\
	m.ind[m.n] = k;				\
      } else if(p[k] == v) {				\
	m.n++;						\
	m.ind[m.n] = k;				\
      }						\
    }							\
  }

// From bcftools.h - https://github.com/samtools/bcftools/blob/b0376dff1ed70603c9490802f37883b9009215d2/bcftools.h#L48
static inline char gt2iupac(char a, char b)
{
  static const char iupac[4][4] = { {'A','M','R','W'},{'M','C','S','Y'},{'R','S','G','K'},{'W','Y','K','T'} };
  if ( a>='a' ) a -= 'a' - 'A';
  if ( b>='a' ) b -= 'a' - 'A';
  if ( a=='A' ) a = 0;
  else if ( a=='C' ) a = 1;
  else if ( a=='G' ) a = 2;
  else if ( a=='T' ) a = 3;
  else return 'N';
  if ( b=='A' ) b = 0;
  else if ( b=='C' ) b = 1;
  else if ( b=='G' ) b = 2;
  else if ( b=='T' ) b = 3;
  else return 'N';
  return iupac[(int)a][(int)b];
}

char get_maj_nuc_code(char **allele, maj m){
  char nuc = allele[m.ind[0]][0];
  for (int k = 1; k < m.n ; ++k){
    nuc = gt2iupac(nuc, allele[m.ind[k + 1]][0]);
  }
  return nuc;
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
    int ad_id;
    //Get the header
    bcf_hdr_t *header = bcf_hdr_read(in);
    if(header == NULL) {
      bcf_close(in);
      throw std::runtime_error("Unable to open BCF header.");
    }
    ad_id = bcf_hdr_id2int(header, BCF_DT_ID, "AD");
    std::cout << ad_id << " " << std::endl;
    //Initialize iterator
    bcf1_t *v = bcf_init();
    int ctr = 0;
    int prev_pos = 0;
    int s;
    maj m;
    while(bcf_read(in, header, v) ==0){
      // std::cout << "Position: " << v->pos << std::endl;
      // std::cout << "Number of Alleles: " << v->n_allele << std::endl;
      for (int i = prev_pos; i < v->pos - 1; ++i){
	out << "N";
      }
      bcf_unpack(v, BCF_UN_FMT);
      bcf_unpack(v, BCF_UN_STR);
      bcf_fmt_t *fmt = v->d.fmt;
      for (int i = 0; i < v->n_fmt; ++i){
	if(strcmp(header->id[BCF_DT_ID][fmt[i].id].key, "AD") == 0){
	  m.ind = (int *) malloc(fmt[i].n * sizeof(int));
	  m.n = -1;
	  switch(fmt[i].type){
	  case BCF_BT_INT8: GET_MAX_INDICES(int8_t, fmt[i].p, fmt[i].n, m); break;
	  case BCF_BT_INT16: GET_MAX_INDICES(int16_t, fmt[i].p, fmt[i].n, m); break;
	  case BCF_BT_INT32: GET_MAX_INDICES(int32_t, fmt[i].p, fmt[i].n, m); break;
	  case BCF_BT_FLOAT: GET_MAX_INDICES(float, fmt[i].p, fmt[i].n, m); break;
	  default: hts_log_error("Unexpected type %d", fmt[i].type); exit(1);
	  }
	  // maj_ind = get_majority(fmt[i], tmp);
	  if(m.n == -1){
	    out << "N";
	  } else {
	    if(bcf_get_variant_types(v) == VCF_REF){
	      out << v->d.allele[0];
	    } else {
	      if(bcf_get_variant_types(v) == VCF_SNP || bcf_get_variant_types(v) == VCF_MNP){
		char nuc = get_maj_nuc_code(v->d.allele, m);
		out << nuc;
		std::cout << nuc << std::endl;
	      } else {
		// TODO
		int d, ind, max_len = 0;
		for (int k = 0; i <= m.n; ++k){
		  if(strlen(v->d.allele[m.ind[k]]) > max_len)
		    max_len = m.ind[k];
		}
		for (int l = 0; l < max_len; ++l){
		  for (int k = 0; k < m.n; ++k){
		    d = strlen(v->d.allele[0]) - strlen(v->d.allele[m.ind[k]]); // Ref - Allele
		    if(d > 0) { // Deletion
		      s = strlen(v->d.allele[0]);
		      ind = 0;//Reference Allele
		    } else if(d < 0) { // Insertion
		      s = strlen(v->d.allele[m.ind[k]]);
		      d = -1 * d;
		      ind = m.ind[k];
		    }
		    for (int j = d; j >0; --j){
		      out << v->d.allele[0][s-d];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      prev_pos = v->pos;
      ctr++;
      if(ctr % 1000 == 0){
	std::cout << ctr << std::endl;
	break;
      }
    }
    bcf_destroy(v);
    bcf_hdr_destroy(header);
    bcf_close(in);
    out.close();
  }
  return 0;
}
