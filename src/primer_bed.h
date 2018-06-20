#include<string>
#include<vector>

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

 public:
  std::string get_name();
  std::string get_region();
  int get_score();
  unsigned int get_start();
  unsigned int get_end();
  char get_strand();
  int get_length();
  void set_start(unsigned int s);
  void set_end(unsigned int e);
  void set_strand(char s);
  void set_region(std::string r);
  void set_name(std::string n);
  void set_score(int s);
};

std::vector<primer> populate_from_file(std::string path);

#endif
