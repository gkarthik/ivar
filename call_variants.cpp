#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>
#include<htslib/kfunc.h>

#include "allele_functions.h"

const char gap='N';

int call_variants_from_plup(std::istream &cin, std::string out_file, uint8_t min_qual){
  std::string line, cell;
  std::cout << "Min Qual: " << (uint16_t)min_qual << std::endl;
  std::ofstream fout(out_file+".tsv");
  fout << "POS\tREF\tALT\tAD\tRAD\tDP\tQUAL"<<std::endl;
  int ctr = 0, pos = 0, mdepth = 0;
  std::stringstream lineStream;
  char ref;
  std::string bases;
  std::string qualities;
  std::vector<allele> ad;
  while (std::getline(cin, line)){
    lineStream << line;
    ctr = 0;
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	break;
      case 1:
	pos = stoi(cell);
	break;
      case 2:
	ref = cell[0];
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
    ad = update_allele_depth(ref, bases, qualities, min_qual);
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if((it->nuc[0] == ref && it->nuc.length() == 1) || it->nuc[0]=='*')
	continue;
      fout << pos << "\t";
      fout << ref << "\t";
      fout << it->nuc << "\t";
      fout << it->depth << "\t";
      fout << it->reverse << "\t";
      fout << mdepth << "\t";
      fout << (uint16_t)it->mean_qual << std::endl;
    }
    lineStream.clear();
  }
  fout.close();
  return 0;
}

// int main_old(int argc, char* argv[]) {
//   std::string line, cell;
//   std::string out_file = argv[1];
//   uint8_t min_qual = 20;
//   if(argc > 2)
//     min_qual = atoi(argv[2]);
//   std::cout << "Min Qual: " << (uint16_t)min_qual << std::endl;
//   std::ofstream fout(out_file+".tsv");
//   fout << "POS\tREF\tALT\tAD\tRAD\tDP\tQUAL"<<std::endl;
//   int ctr = 0, pos = 0, mdepth = 0;
//   std::stringstream lineStream;
//   char ref;
//   std::string bases;
//   std::string qualities;
//   std::vector<allele> ad;
//   while (std::getline(std::cin, line)){
//     lineStream << line;
//     ctr = 0;
//     while(std::getline(lineStream,cell,'\t')){
//       switch(ctr){
//       case 0:
// 	break;
//       case 1:
// 	pos = stoi(cell);
// 	break;
//       case 2:
// 	ref = cell[0];
// 	break;
//       case 3:
// 	mdepth = stoi(cell);
// 	break;
//       case 4:
// 	bases = cell;
// 	break;
//       case 5:
// 	qualities = cell;
// 	break;
//       case 6:
// 	break;
//       }
//       ctr++;
//     }
//     ad = update_allele_depth(ref, bases, qualities, min_qual);
//     for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
//       if((it->nuc[0] == ref && it->nuc.length() == 1) || it->nuc[0]=='*')
// 	continue;
//       fout << pos << "\t";
//       fout << ref << "\t";
//       fout << it->nuc << "\t";
//       fout << it->depth << "\t";
//       fout << it->reverse << "\t";
//       fout << mdepth << "\t";
//       fout << (uint16_t)it->mean_qual << std::endl;
//     }
//     lineStream.clear();
//   }
//   fout.close();
//   return 0;
// }
