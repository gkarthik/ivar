#include<iostream>
#include<fstream>
#include<string>
#include "../src/allele_functions.h"
#include "../src/call_variants.h"

int call_var_check_outfile(std::string prefix, uint8_t min_qual, uint8_t min_depth, double min_threshold, std::string out1, std::string out2){
  std::string path = "../data/test.indel.mpileup";
  std::ifstream mplp(path);
  call_variants_from_plup(mplp, prefix, min_qual, min_threshold, min_depth);
  std::ifstream outFile(prefix+".tsv");
  std::string l;
  getline(outFile, l);		// Ignore first line
  int comp = 0;
  getline(outFile, l);
  comp += l.compare(out1);
  std::cout << out1 << std::endl;
  std::cout << l << std::endl;
  if(!out2.empty()){
    getline(outFile, l);
    comp += l.compare(out2);
  }
  return comp;
}

int main() {
  int num_success = 0;
  // Quality threshold 20. Frequency threshold: 0.03. Total_DP = 3. Indel passes filters with total_depth 4. Has two lines.
  std::string t_20_02_1 = "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0.2\tFALSE";
  std::string t_20_02_2 = "test\t210\tA\t+GT\t1\t1\t41\t1\t0\t20\t0.25\t4\t0.4\tFALSE";
  // Quality threshold 3-. Frequency threshold: 0.03. Total_DP = 3. Freq = 0.666667 No Indel
  std::string t_20_03 = "test\t210\tA\tT\t1\t1\t41\t2\t1\t40\t0.666667\t3\t0.2\tFALSE";
  // Quality threshold 25. Frequency threshold: 0.03. Total_DP = 2. Freq = 0.5 No Indel.
  std::string t_25_03 = "test\t210\tA\tT\t1\t1\t41\t1\t1\t58\t0.5\t2\t0.4\tFALSE";
  // Minimum depth threshold. Should be empty
  std::string t_25_03_20 = "";
  num_success = call_var_check_outfile("../data/test.indel", 20, 0, 0.02, t_20_02_1, t_20_02_2);
  std::cout << num_success << std::endl;
  num_success += call_var_check_outfile("../data/test.indel", 20, 0, 0.03, t_20_03, "");
  std::cout << num_success << std::endl;
  num_success += call_var_check_outfile("../data/test.indel", 25, 0, 0.03, t_25_03, "");
  std::cout << num_success << std::endl;
  num_success += call_var_check_outfile("../data/test.indel", 25, 20, 0.03, t_25_03_20, "");
  std::cout << num_success << std::endl;
  if(num_success == 0)
    return 0;
  return -1;
}
