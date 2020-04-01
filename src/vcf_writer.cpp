#include "vcf_writer.h"

int vcf_writer::init_header(){
  this->hdr = bcf_hdr_init("w");
  int res;
  res = bcf_hdr_append(this->hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
  res = bcf_hdr_append(this->hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total read depth for each allele on forward strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total read depth for each allele on reverse strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)\">");
  res = bcf_hdr_append(this->hdr, ("##contig=<ID="+this->region+">").c_str());
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Total read depth for each allele on forward strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Total read depth for each allele on reverse strand (based on minimum quality threshold)\">");
  res = bcf_hdr_append(this->hdr, "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)\">");
  res = bcf_hdr_append(this->hdr, "##FILTER=<ID=fobs,Description=\"Fisher's exact test to check if frequency of the iSNV is significantly higher than the mean error rate at that position\">");
  bcf_hdr_add_sample(this->hdr, this->sample_name.c_str());
  if(res != 0)
    return -1;
  return 0;
}

vcf_writer::~vcf_writer(){
  bcf_hdr_destroy(this->hdr);
  bcf_close(this->file);
}

vcf_writer::vcf_writer(char _mode, std::string fname, std::string region, std::string sample_name, std::string ref_path){
  this->ref = new ref_antd(ref_path);
  this->region = region;
  this->sample_name = sample_name;
  std::string mode = "w";
  mode += _mode;
  this->file = bcf_open(fname.c_str(), mode.c_str());
  this->init_header();
  if(this->hdr == NULL || this->file == NULL){
    std::cout << "Unable to write BCF/VCF file" << std::endl;
  }
  if(bcf_hdr_write(this->file, this->hdr) < 0){
    std::cout << "Unable to write BCF/VCF file" << std::endl;
  }
}

int vcf_writer::write_record(uint32_t pos, allele aalt, allele aref, uint32_t depth){
  bcf1_t *rec = bcf_init();
  rec->rid = bcf_hdr_name2id(this->hdr, this->region.c_str());
  rec->pos  = pos - 1;		// converts to pos + 1 on write
  rec->qual = aref.mean_qual;
  bcf_update_id(this->hdr, rec, ".");
  std::string allele_str;
  if(aalt.nuc[0]=='-'){
    int del_len = 0;
    while(del_len < aalt.nuc.length() - 1){
      allele_str += this->ref->get_base(pos + del_len, this->region);
      del_len++;
    }
  } else if (aalt.nuc[0] == '+') {
    allele_str = aref.nuc + "," + aref.nuc + aalt.nuc.substr(1);
  } else {
    allele_str = aref.nuc +"," + aalt.nuc;
  }
  int res;
  bcf_update_alleles_str(this->hdr, rec, allele_str.c_str());
  // INFO
  int32_t *tmp = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
  int32_t tmpi = 1;
  // NS
  bcf_update_info_int32(this->hdr, rec, "DP", &(depth), bcf_hdr_nsamples(this->hdr));
  bcf_update_info_int32(this->hdr, rec, "NS", &tmpi, 1);
  tmpi = bcf_hdr_id2int(this->hdr, BCF_DT_ID, "PASS");
  bcf_update_filter(this->hdr, rec, &tmpi, 1);
  bcf_update_info_int32(this->hdr, rec, "AD", &(aref.depth), bcf_hdr_nsamples(this->hdr));
  tmpi = aref.depth - aref.reverse;
  bcf_update_info_int32(this->hdr, rec, "ADF", &tmpi, bcf_hdr_nsamples(this->hdr));
  bcf_update_info_int32(this->hdr, rec, "ADR", &(aref.reverse), bcf_hdr_nsamples(this->hdr));
  float tmpf = (float)aref.depth/(float)depth;
  bcf_update_info_float(this->hdr, rec, "AF", &tmpf, bcf_hdr_nsamples(this->hdr));
  // Set genotype to ref since it crossed threshold.
  tmp[0] = bcf_gt_phased(1);
  bcf_update_genotypes(this->hdr, rec, tmp, 1);
  // FORMAT
  tmpi = aalt.depth - aalt.reverse;
  bcf_update_format_int32(this->hdr, rec, "AD", &(aalt.depth), bcf_hdr_nsamples(this->hdr));
  bcf_update_format_int32(this->hdr, rec, "GQ", &(aalt.mean_qual), bcf_hdr_nsamples(this->hdr));
  bcf_update_format_int32(this->hdr, rec, "ADF", &tmpi, bcf_hdr_nsamples(this->hdr));
  bcf_update_format_int32(this->hdr, rec, "ADR", &(aalt.reverse), bcf_hdr_nsamples(this->hdr));
  tmpf = (float)aalt.depth/(float)depth;
  bcf_update_format_float(this->hdr, rec, "AF", &tmpf, bcf_hdr_nsamples(this->hdr));
  res = bcf_write1(this->file, this->hdr, rec);
  if(res != 0)
    std::cout << "Unable to write to VCF/BCF file!" << std::endl;
  delete tmp;
 return res;
}

int main(int argc, char *argv[])
{
  vcf_writer *vw = new vcf_writer(0, "./test.vcf", "test", "test_sample", "../data/db/test_ref.fa");
  allele r = {
  depth: 5,
  reverse: 2,
  nuc: "C",
  mean_qual:30
  };
  allele a = {
  depth: 10,
  reverse: 6,
  nuc: "-TA",
  mean_qual:25
  };
  int res = vw->write_record(5, a, r, 15);
  delete vw;
  return 0;
}
