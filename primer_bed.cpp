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
