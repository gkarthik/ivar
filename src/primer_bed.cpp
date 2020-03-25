#include "primer_bed.h"

std::string primer::get_name(){
  return name;
}

std::string primer::get_region(){
  return region;
}

int primer::get_score(){
  return score;
}

unsigned int primer::get_start(){
  return start;
}

unsigned int primer::get_end(){
  return end;
}

char primer::get_strand(){
  return strand;
}

int primer::get_length(){
  return end - start + 1;
}

int primer::get_pair_indice(){
  return pair_indice;
}

void primer::set_start(unsigned int s){
  start = s;
}

void primer::set_end(unsigned int e){
  end = e;
}

void primer::set_strand(char s){
  strand = s;
}

void primer::set_region(std::string r){
  region = r;
}

void primer::set_name(std::string n){
  name = n;
}

void primer::set_score(int s){
  score = s;
}

void primer::set_pair_indice(int i){
  pair_indice = i;
}

void print_bed_format(){
  std::cout << "iVar uses the standard 6 column BED format as defined here - https://genome.ucsc.edu/FAQ/FAQformat.html#format1." << std::endl;
  std::cout << "It requires the following columns delimited by a tab: chrom, chromStart, chromEnd, name, score, strand" << std::endl;
}

std::vector<primer> populate_from_file(std::string path){
  std::ifstream  data(path.c_str());
  std::string line;
  std::vector<primer> primers;
  while(std::getline(data,line)){ // Remove extra lineStream
    std::stringstream lineStream(line);
    std::string cell;
    int ctr = 0;
    primer p;
    while(std::getline(lineStream,cell,'\t')){
      switch(ctr){
      case 0:
	p.set_region(cell);
	break;
      case 1:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_start(std::stoul(cell));
	} else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
	break;
      case 2:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_end(std::stoul(cell)-1); // Bed format - End is not 0 based
	} else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
	break;
      case 3:
	p.set_name(cell);
	break;
      case 4:
	if(std::all_of(cell.begin(), cell.end(), ::isdigit)) {
	  p.set_score(stoi(cell));
	} else {
	  print_bed_format();
	  primers.clear();
	  return primers;
	}
	break;
      case 5:
	p.set_strand(cell[0]);
      }
      ctr++;
    }
    p.set_pair_indice(-1);
    primers.push_back(p);
  }
  if(primers.size() == 0){
    print_bed_format();
    std::cout << "Found 0 primers in BED file" << std::endl;
  }
  std::cout << "Found " << primers.size() << " primers in BED file" << std::endl;
  return primers;
}

int get_primer_indice(std::vector<primer> p, unsigned int pos){
  for(std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if(it->get_start() <= pos && it->get_end() >= pos){
      return it - p.begin();
    }
  }
  return -1;
}

int get_primer_indice(std::vector<primer> p, std::string name){
  for(std::vector<primer>::iterator it = p.begin(); it != p.end(); ++it) {
    if(it->get_name().compare(name) == 0){
      return it - p.begin();
    }
  }
  return -1;
}

int populate_pair_indices(std::vector<primer> &primers, std::string path){
  std::ifstream fin(path.c_str());
  std::string line, cell, p1,p2;
  std::stringstream line_stream;
  std::vector<primer>::iterator it;
  int indice = -1;
  while (std::getline(fin, line)){
    line_stream << line;
    std::getline(line_stream, cell, '\t');
    p1 = cell;
    line_stream.clear();
    std::getline(line_stream, cell, '\t');
    p2 = cell;
    line_stream.clear();
    if(!p1.empty() && !p2.empty()){
      for(it = primers.begin(); it != primers.end(); ++it) {
	if (it->get_name() == p1) {
	  indice = get_primer_indice(primers, p2);
	  if (indice != -1)
	    it->set_pair_indice(indice);
	  else
	    std::cout << "Primer pair for " << p1 << " not found in BED file." <<std::endl;
	} else if (it->get_name() == p2){
	  indice = get_primer_indice(primers, p1);
	  if(indice != -1)
	    it->set_pair_indice(indice);
	  else
	    std::cout << "Primer pair for " << p2 << " not found in BED file." << std::endl;
	}
      }
    }
  }
  return 0;
}
