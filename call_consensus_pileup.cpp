#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>

struct allele{
  std::string nuc;
  uint32_t depth;
  uint32_t reverse;
};

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

std::string get_consensus_allele(std::vector<allele> ad){
  if(ad.size()==0)
    return "";
  if(ad.size() == 1)
    return ad.at(0).nuc;
  std::string cnuc = "";
  char n;
  int max_l = 0;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end()-1; ++it) {
    if(it->nuc.length() > max_l){
      max_l = it->nuc.length();
    }
  }
  for (int i = 0; i < max_l; ++i){
    n = '*';
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end()-1; ++it) {
      if(i < it->nuc.length() && i < (it+1)->nuc.length()){
	n = gt2iupac(it->nuc[i], (it+1)->nuc[i]);
      } else if(i < it->nuc.length()){
	n = it->nuc[i];
      } else if(i < (it+1)->nuc.length()){
	n = (it+1)->nuc[i];
      }
    }
    if(n!='*')
      cnuc += n;
  }
  return cnuc;
}

std::vector<allele> get_major_alleles(std::vector<allele> ad){
  std::vector<allele> maj_ad;
  uint32_t max = 0;
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc.compare("*") == 0)
      continue;
    if(it->depth > max){
      maj_ad.clear();
      maj_ad.push_back(ad.at(it-ad.begin()));
      max = it->depth;
    } else if(it->depth == max){
      maj_ad.push_back(ad.at(it-ad.begin()));
    }
  }
  return maj_ad;
}

void print_allele_depths(std::vector<allele> ad){
  std::cout << "Print AD" << " ";
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    std::cout << it->nuc << " ";
    std::cout << it->depth << " ";
    std::cout << it->reverse;
  }
  std::cout << std::endl;
}

int check_allele_exists(std::string n, std::vector<allele> ad){
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    if(it->nuc.compare(n) == 0){
      return it - ad.begin();
    }
  }
  return -1;
}

std::vector<allele> update_allele_depth(char ref,std::string bases, std::string qualities){
  std::vector<allele> ad;
  int i = 0, n =0, j = 0;
  while (i < bases.length()){
    std::string b;
    allele tmp;
    bool forward= true;
    switch(bases[i]){
    case '.':
      b = ref;
      break;
    case ',':
      b = ref;
      forward = false;
      break;
    case '*':
      b = bases[i];
      break;
    default:
      int asc_val = bases[i];
      if(asc_val >= 65 && asc_val <= 90){
	b = bases[i];
      } else if(asc_val>=97 && asc_val<=122){
	b = bases[i] - 32;
	forward = false;
      } else {
	i++;
	continue;
      }
    }
    if(bases[i+1]=='+' || bases[i+1]=='-'){		// Deletions are ignored since subsequent bases take care of bases
      j = i+2;
      while(isdigit(bases[j])){
	j++;
      }
      j = j - (i+2);
      n = stoi(bases.substr(i+2, j));
      if(bases[i+1]== '+')
	b += bases.substr(i+3, n);
      i += n + 2;
    }
    int ind = check_allele_exists(b, ad);
    if (ind==-1){
      tmp.nuc = b;
      tmp.depth = 1;
      if(!forward)
	tmp.reverse = 1;
      else
	tmp.reverse = 0;
      ad.push_back(tmp);
    } else {
      ad.at(ind).depth += 1;
      if(!forward)
	ad.at(ind).reverse += 1;
    }
    i++;
  }
  return ad;
}

int main(int argc, char* argv[]) {
  std::string line, cell;
  std::string out_file = argv[1];
  std::ofstream fout(out_file+".fa");
  fout << ">Consensus"<<std::endl;
  std::cout << "Lines: " << std::endl;
  int ctr = 0;
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
	break;
      case 2:
	ref = cell[0];
	break;
      case 3:
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
    ad = update_allele_depth(ref, bases, qualities);
    // print_allele_depths(ad);
    std::vector<allele> mad = get_major_alleles(ad);
    fout << get_consensus_allele(mad);
    // std::cout << "Consensus: " << get_consensus_allele(mad) << std::endl;
    lineStream.clear();
  }
  fout.close();
  return 0;
}
