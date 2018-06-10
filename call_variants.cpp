#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<string>
#include<regex>
#include<htslib/kfunc.h>

struct allele{
  std::string nuc;
  uint32_t depth;
  uint32_t reverse;
  uint8_t mean_qual;
  bool operator < (const allele& a) const{
    return (nuc.compare(a.nuc) > 0) ? true : false;
  }
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
  std::string indel;
  int i = 0, n =0, j = 0, q_ind = 0;
  uint8_t q;
  while (i < bases.length()){
    if(bases[i] == '^'){
      i += 2;			// Skip mapping quality as well (i+1) - 33
      continue;
    }
    if(bases[i] == '$'){
      i++;
      continue;
    }
    std::string b;
    allele tmp;
    bool forward= true;
    q = qualities[q_ind] - 33;
    switch(bases[i]){
    case '.':
      b = ref;
      break;
    case ',':
      b = ref;
      forward = false;
      break;
    // case '*':
    //   b = bases[i];
    //   break;
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
      indel = bases.substr(i+2+j, n);
      transform(indel.begin(), indel.end(), indel.begin(),::toupper);
      b += indel;
      i += n + 2;
    }
    int ind = check_allele_exists(b, ad);
    if (ind==-1){
      tmp.nuc = b;
      tmp.depth = 1;
      tmp.mean_qual = q;
      if(!forward)
	tmp.reverse = 1;
      else
	tmp.reverse = 0;
      ad.push_back(tmp);
    } else {
      ad.at(ind).mean_qual = ((ad.at(ind).mean_qual * ad.at(ind).depth) + q)/(ad.at(ind).depth + 1);
      ad.at(ind).depth += 1;
      if(!forward)
	ad.at(ind).reverse += 1;
    }
    i++;
    q_ind++;
  }
  std::sort(ad.begin(), ad.end());
  return ad;
}

int main(int argc, char* argv[]) {
  std::string line, cell;
  std::string out_file = argv[1];
  std::ofstream fout(out_file+".tsv");
  fout << "Position\tRef\tAlt\tAllele Depth\tReverse Reads Depth\tTotal Depth\tMean Quality"<<std::endl;
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
    ad = update_allele_depth(ref, bases, qualities);
    for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
      if(it->nuc[0] == ref && it->nuc.length() == 1)
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
