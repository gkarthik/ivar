#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

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
	p.set_start(std::stoul(cell));
	break;
      case 2:
	p.set_end(std::stoul(cell)-1); // Bed format - End is not 0 based
	break;
      case 3:
	p.set_name(cell);
	break;
      case 4:
	p.set_score(atoi(cell.c_str()));
	break;
      case 5:
	p.set_strand(cell[0]);
      }
      ctr++;
    }
    primers.push_back(p);
  }
  return primers;
}
