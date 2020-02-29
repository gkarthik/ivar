#include <iostream>
#include <fstream>
#include <sstream>

int read_gff(std::string file_path){
  std::string line, cell;
  std::stringstream line_stream;
  std::ifstream fin = std::ifstream(file_path);
  int ctr = 0;
  while (std::getline(fin, line)){
    line_stream << line;
    ctr = 0;
    while(std::getline(line_stream,cell,'\t')){
      switch(ctr){
      case 0:

	break;
      case 1:

	break;
      case 2:

	break;
      case 3:

	break;
      case 4:
	break;
      case 5:
	break;
      case 6:
	break;
      }
      ctr++;
    }
    line_stream.clear();
  }
  return 0;
}
