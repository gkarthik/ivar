#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <map>
#include <iterator>

#include "parse_gff.h"

gff3_feature::gff3_feature(std::string line){
  int ctr = 0;
  std::stringstream line_stream;
  std::string cell;
  line_stream << line;
  while(std::getline(line_stream,cell,'\t')){
    switch(ctr){
    case 0:			// Name
      this->seqid = cell;
      break;
    case 1:
      this->source = cell;
      break;
    case 2:
      this->type = cell;
      break;
    case 3:
      this->start = atoi(cell.c_str());
      break;
    case 4:
      this->end = atoi(cell.c_str());
      break;
    case 5:
      this->score = atof(cell.c_str());
      break;
    case 6:
      this->strand = cell[0];
      break;
    case 7:
      this->phase = cell[0];
      break;
    case 8:
      this->set_attributes(cell);
      break;
    }
    ctr++;
  }
  line_stream.clear();
}

int gff3_feature::print(){
  std::cout << seqid << "\t"
	    << source << "\t"
    	    << type << "\t"
    	    << start << "\t"
	    << end << "\t"
	    << score << "\t"
	    << strand << "\t"
	    << phase << "\t";
  std::map<std::string, std::string>::iterator it;
  for (it = attributes.begin(); it != attributes.end(); it++){
    std::cout << it->first << ": " << it->second << "; ";
  }
  std::cout << std::endl;
  return 0;
}

std::string gff3_feature::get_attr(std::string key){
  std::string val;
  if(attributes.find(key) != attributes.end()){
    val = attributes[key];
  }
  return val;
}

int gff3_feature::set_attributes(std::string attr){
  std::string key, val;
  std::regex exp("[^;]+");
  std::regex_iterator<std::string::iterator> it(attr.begin(), attr.end(), exp);
  std::regex_iterator<std::string::iterator> rend;
  std::string delimiter = "=";
  while (it!=rend) {
    key = it->str().substr(0, it->str().find(delimiter));
    val = it->str().substr(it->str().find(delimiter)+1, it->str().length());
    if(!key.empty() && !val.empty()){
      this->attributes[key] = val;
    }
    ++it;
  }
  return 0;
}

int gff3::print(){
  std::vector<gff3_feature>::iterator it;
  for(it = features.begin(); it != features.end(); it++){
    it->print();
  }
  return 0;
}

std::vector<gff3_feature> gff3::get_features(){
  return features;
}

gff3::gff3(){}

gff3::gff3(std::string path){
  this->read_file(path);
}

int gff3::read_file(std::string path){
  std::ifstream fin = std::ifstream(path);
  if(!fin){
    std::cout << "GFF file does not exist at " << path << std::endl;
    return -1;
  }
  std::string line;
  while (std::getline(fin, line)){
    if(line[0] == '#' && line[1] == '#') // Avoid comments in GFF file
      continue;
    features.push_back(gff3_feature(line));
  }
  return 0;
}
