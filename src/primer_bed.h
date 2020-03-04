#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#ifndef primer_bed
#define primer_bed

class primer {
 private:
  std::string region;
  unsigned int start;
  unsigned int end;
  std::string name;
  int score;
  char strand;
  int pair_indice;

 public:
  std::string get_name();
  std::string get_region();
  int get_score();
  unsigned int get_start();
  unsigned int get_end();
  char get_strand();
  int get_length();
  int get_pair_indice();
  void set_start(unsigned int s);
  void set_end(unsigned int e);
  void set_strand(char s);
  void set_region(std::string r);
  void set_name(std::string n);
  void set_score(int s);
  void set_pair_indice(int i);
};

std::vector<primer> populate_from_file(std::string path);
int get_primer_indice(std::vector<primer> p, unsigned int pos);
int get_primer_indice(std::vector<primer> p, std::string name);
int populate_pair_indices(std::vector<primer> &primers, std::string path);

#endif
