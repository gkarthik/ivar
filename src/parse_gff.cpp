#include <iostream>
#include <fstream>
#include <sstream>

/* 
GFF3 file format:
| Col        | Desc                                                          |
| seqid      | ID of ref                                                     |
| source     | algorithm/operating procedure                                 |
| type       | type of feature                                               |
| start      | 1-based start                                                 |
| end        | start <= end                                                  |
| score      | score of feature                                              |
| strand     | +/-                                                           |
| phase      | Number of bases to be removed with reference to reading frame |
| attributes | list of attributes in the format tag=value;tag-value..        |
*/

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
      case 7:
	break;
      case 8:
	break;
      default:
	break;
      }
      ctr++;
    }
    line_stream.clear();
  }
  return 0;
}
