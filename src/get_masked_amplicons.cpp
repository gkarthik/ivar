#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "primer_bed.h"

int get_primer_indice(std::vector<primer> p, unsigned int pos){
  for(std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if(it->get_start() <= pos && it->get_end() >= pos){
      return it - p.begin();
    }
  }
  return -1;
}

int get_primers_with_mismatches(std::string bed, std::string vpath, std::string out){
  std::vector<primer> primers = populate_from_file(bed);
  std::vector<unsigned int> indices;
  std::string line, cell;
  std::ifstream fin(vpath.c_str());
  out += ".txt";
  std::ofstream fout(out.c_str());
  unsigned int ctr, pos;
  int ind;
  std::stringstream  line_stream;
  while (std::getline(fin, line)){
    line_stream << line;
    ctr = 0;
    pos = 0;
    while(std::getline(line_stream,cell,'\t')){
      switch(ctr){
      case 1:
	if(cell != "POS")
	  pos = atoi(cell.c_str());
	break;
      default:
	break;
      }
      ctr++;
    }
    if(pos == 0){
      line_stream.clear();
      continue;
    }
    ind = get_primer_indice(primers, pos);
    if(std::find(indices.begin(), indices.end(), ind) == indices.end() && ind != -1){
      indices.push_back(ind);
    }
    line_stream.clear();
  }
  for(std::vector<unsigned int>::iterator it = indices.begin(); it != indices.end(); ++it) {
    fout << *it;
    std::cout << *it;
    if(it != indices.end() - 1){
      fout << " ";
      std::cout << " ";
    }
  }
  std::cout << std::endl;
  return 0;
}
