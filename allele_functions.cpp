#include<string>
#include<vector>
#include<algorithm>
#include<iostream>

#include "allele_functions.h"

void print_allele_depths(std::vector<allele> ad){
  std::cout << "AD Size: " << ad.size() << " ";
  for(std::vector<allele>::iterator it = ad.begin(); it != ad.end(); ++it) {
    std::cout << it->nuc << " ";
    std::cout << it->depth << " ";
    std::cout << it->reverse << " ";
    char o = (it->prev_base!=0) ? it->prev_base : '0';
    std::cout << o  << "\t";
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

int find_ref_in_allele(std::vector<allele> ad, char ref){
  std::string ref_s(1, ref);
  std::vector<allele>::iterator it = ad.begin();
  while(it < ad.end()){
    if(it->nuc.compare(ref_s) == 0)
      return (it - ad.begin());
    it++;
  }
  return -1;
}

std::vector<allele> update_allele_depth(char ref,std::string bases, std::string qualities, uint8_t min_qual){
  std::vector<allele> ad;
  std::string indel;
  int i = 0, n =0, j = 0, q_ind = 0;
  uint8_t q;
  char prev_base;
  while (i < bases.length()){
    if(bases[i] == '^'){
      i += 2;			// Skip mapping quality as well (i+1) - 33
      continue;
    }
    if(bases[i] == '$'){
      i++;
      continue;
    }
    q = qualities[q_ind] - 33;
    if(q < min_qual){
      i++;
      q_ind++;
      continue;
    }
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
    case '+': case '-':
      j = i+1;
      while(isdigit(bases[j])){
	j++;
      }
      j = j - (i+1);
      n = stoi(bases.substr(i+1, j));
      indel = bases.substr(i+1+j, n);
      transform(indel.begin(), indel.end(), indel.begin(),::toupper);
      b = bases[i] + indel;	// + for Insertion and - for Deletion
      i += n + 1;
      if(indel[0]>=97 && indel[0] <= 122)
	forward=false;
      break;
    default:
      int asc_val = bases[i];
      if(asc_val >= 65 && asc_val <= 90){
	b = bases[i];
      } else if(asc_val>=97 && asc_val<=122){
	b = bases[i] - 32;
	forward = false;
      }
    }
    int ind = check_allele_exists(b, ad);
    if(q >= min_qual){
      if (ind==-1){
	tmp.prev_base = 0;
	if(b[0]=='+' || b[0]=='-')
	  tmp.prev_base = prev_base;
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
      prev_base = b[0];
    }
    i++;
    if(b[0] !='+' && b[0]!='-')
      q_ind++;
  }
  std::sort(ad.begin(), ad.end());
  return ad;
}
