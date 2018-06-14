#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>

#include "allele_functions.h"

// struct allele{
//   std::string nuc;
//   uint32_t depth;
//   uint32_t reverse;
//   bool operator < (const allele& a) const{
//     return (nuc.compare(a.nuc) > 0) ? true : false;
//   }
// };

const char gap='N';

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

void format_alleles(std::vector<allele> &ad, char ref){
  std::vector<allele>::iterator it = ad.begin();
  uint32_t rdepth = 0;
  while(it < ad.end()){
    if(it->nuc[0] == '-'){	// Remove deletions
      it = ad.erase(it);
    } else if (it->nuc[0] == '+'){
      it->nuc[0] = ref;
      rdepth += it->depth;
    }
    ++it;
  }
  int ref_ind = find_ref_in_allele(ad, ref);
  std::cout << ref_ind << std::endl;
  if(ref_ind!=-1)
    ad.at(ref_ind).depth -= rdepth;
}

std::string get_consensus_allele(std::vector<allele> ad, char ref){
  std::cout << "Consensus ";
  print_allele_depths(ad);
  format_alleles(ad, ref);
  print_allele_depths(ad);
  if(ad.size()==0)
    return "";
  if(ad.size() == 1)
    return (ad.at(0).nuc.compare("*") == 0) ? "" : ad.at(0).nuc;
  std::string cnuc = "";
  char n;
  int max_l = 0, mdepth = 0, tdepth = 0;
  uint32_t gaps = 0;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end()-1; ++it) {
    if(it->nuc.length() > max_l){
      max_l = it->nuc.length();
    }
  }
  for (int i = 0; i < max_l; ++i){
    n = '*';
    mdepth = 0;
    tdepth = 0;
    std::vector<allele>::iterator it = ad.begin();
    gaps = 0;
    while(it!=ad.end()){
      if(it->nuc[i]=='+'){
	it++;
	continue;
      }
      if(!(i < it->nuc.length())){
	gaps += it->depth;
	it++;
	continue;
      }
      tdepth = it->depth;
      while(it!=ad.end() -1 && it->nuc[i] == (it+1)->nuc[i] && (i+1) < (it+1)->nuc.length()){ // Third condition for cases with more than 1 base in insertion.
	tdepth += (it+1)->depth;
	it++;
      }
      if(tdepth > mdepth){
	n = it->nuc[i];
	mdepth = tdepth;
      } else if(tdepth == mdepth){
	n = gt2iupac(n, it->nuc[i]);
      }
      it++;
    }
    if(n!='*' && mdepth >= gaps) // TODO: Check what to do when equal.
      cnuc += n;
  }
  std::cout << " " << cnuc << std::endl;
  return cnuc;
}

// int check_allele_exists(std::string n, std::vector<allele> ad){
//   for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
//     if(it->nuc.compare(n) == 0){
//       return it - ad.begin();
//     }
//   }
//   return -1;
// }

// std::vector<allele> update_allele_depth(char ref,std::string bases, std::string qualities){
//   std::vector<allele> ad;
//   std::string indel;
//   int i = 0, n =0, j = 0;
//   while (i < bases.length()){
//     if(bases[i] == '^'){
//       i += 2;			// Skip mapping quality as well (i+1) - 33
//       continue;
//     }
//     if(bases[i] == '$'){
//       i++;
//       continue;
//     }
//     std::string b;
//     allele tmp;
//     bool forward= true;
//     switch(bases[i]){
//     case '.':
//       b = ref;
//       break;
//     case ',':
//       b = ref;
//       forward = false;
//       break;
//     case '*':
//       b = bases[i];
//       break;
//     default:
//       int asc_val = bases[i];
//       if(asc_val >= 65 && asc_val <= 90){
// 	b = bases[i];
//       } else if(asc_val>=97 && asc_val<=122){
// 	b = bases[i] - 32;
// 	forward = false;
//       } else {
// 	i++;
// 	continue;
//       }
//     }
//     if(bases[i+1]=='+' || bases[i+1]=='-'){		// Deletions are ignored since subsequent bases take care of bases
//       j = i+2;
//       while(isdigit(bases[j])){
// 	j++;
//       }
//       j = j - (i+2);
//       n = stoi(bases.substr(i+2, j));
//       if(bases[i+1]== '+'){
// 	indel = bases.substr(i+2+j, n);
// 	transform(indel.begin(), indel.end(), indel.begin(),::toupper);
// 	b += indel;
//       }
//       i += n + 2;
//     }
//     int ind = check_allele_exists(b, ad);
//     if (ind==-1){
//       tmp.nuc = b;
//       tmp.depth = 1;
//       if(!forward)
// 	tmp.reverse = 1;
//       else
// 	tmp.reverse = 0;
//       ad.push_back(tmp);
//     } else {
//       ad.at(ind).depth += 1;
//       if(!forward)
// 	ad.at(ind).reverse += 1;
//     }
//     i++;
//   }
//   std::sort(ad.begin(), ad.end());
//   return ad;
// }

int main(int argc, char* argv[]) {
  std::string line, cell;
  std::string out_file = argv[1];
  uint8_t min_qual = 20;
  uint32_t min_depth = 10;
  if(argc > 1){
    min_qual = atoi(argv[2]);
  }
  if(argc > 2){
    min_depth = atoi(argv[3]);
  }
  std::ofstream fout(out_file+".fa");
  fout << ">Consensus"<<std::endl;
  int ctr = 0, pos = 0, mdepth = 0;
  std::stringstream lineStream;
  char ref;
  std::string bases;
  std::string qualities;
  std::vector<allele> ad;
  while (std::getline(std::cin, line)){
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
    // print_allele_depths(ad);
    fout << get_consensus_allele(ad, ref);
    // std::cout << std::endl << "Consensus: " << get_consensus_allele(ad) << std::endl;
    lineStream.clear();
  }
  fout.close();
  return 0;
}
