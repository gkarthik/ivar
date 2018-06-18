#include<iostream>

/* Substitution Matrix

   A  T  G  C
A  1 -1 -1 -1
T -1  1 -1 -1
G -1 -1  1 -1
C -1 -1 -1  1

*/

const int substitution[4][4] = {
  {3,-3,-3,3},
  {-3,3,-3,-3},
  {-3,-3,3,-3},
  {-3,-3,-3,3}
};

const int gap_open = 2;
const int gap_extension = -1;
const int max_read_size = 500;
const int max_adapter_size = 30;

int get_sub_score(char a, char b){
  int i,j;
  switch(a){
  case 'A':
    i = 0;
    break;
  case 'T':
    i = 1;
    break;
  case 'G':
    i = 2;
    break;
  case 'C':
    i = 3;
    break;
  }
  switch(b){
  case 'A':
    j = 0;
    break;
  case 'T':
    j = 1;
    break;
  case 'G':
    j = 2;
    break;
  case 'C':
    j = 3;
    break;
  }
  return substitution[i][j];
}

int get_gap_penalty(int k){
  return k * gap_open; // Linear Gap Penalty
}

void print_matrix(int h[][max_adapter_size], int r, int c, std::string read, std::string adap){
  for(int i =0; i< r;i++){
    if(i == 0){
      std::cout << "    ";
      for(int k =0; k < c-1; k++) std::cout << adap[k] << " ";
      std::cout << "\n";
    }
    for(int j = 0; j < c;j++){
      if(j == 0 && i > 0){
	if(i -1 < read.length())
	  std::cout << read[i-1] << " ";
      }
      if(i == 0 && j == 0)
	std::cout << " " << " ";
      std::cout << h[i][j] << " ";
    }
    std::cout << "\n";
  }
}

// Return Value 1 - Diag, 2->Right, 3->Down, 0-> 0
int*  get_score_cell(int h[][max_adapter_size], int i, int j, std::string read, std::string adap){
  int s[] = {0,0,0};
  int *tmp = new int[2];
  tmp[0] = 0;
  tmp[1] = 0;
  s[0] = h[i-1][j-1] + get_sub_score(read[i-1], adap[j-1]);
  s[1] = h[i-1][j] - get_gap_penalty(1);
  s[2] = h[i][j-1] - get_gap_penalty(1);
  int max = (s[0]>=s[1]) ? s[0] : s[1];
  max = (max>=s[2]) ? max : s[2];
  max = (max > 0) ? max : 0;
  tmp[0] = max;
  if(max == s[0]){
    tmp[1] = 1;
  } else if(max == s[1]){
    tmp[1] = 3;
  } else if(max == s[2]){
    tmp[1] = 2;
  }
  return tmp;
}

void print_alignment(char a[2][max_read_size], int n){
  for(int j =0;j<2;j++){
    for(int i = n -1; i >= 0;i--){
      std::cout << a[j][i] << " ";
    }
    std::cout << std::endl;
  }
}

int align_seqs(std::string read, std::string adap){
  int h[max_read_size][max_adapter_size];
  int t[max_read_size][max_adapter_size];
  int m = read.length() + 1,
    n=adap.length()+ 1,
    max_i = 0,
    max_j = 0,
    max_v = -1;
  int *tmp = new int[2];
  // Initialize first column to 0
  for(int i = 0; i < n; i++) h[0][i] = 0;
  // Initialize first row to 0
  for(int i = 0; i < m; i++) h[i][0] = 0;
  // print_matrix(h, m, n, read, adap);
  for(int i = 1; i < m; i++){
    for(int j = 1; j < n; j++){
      tmp = get_score_cell(h, i, j, read, adap);
      h[i][j] = *tmp;
      t[i][j] = *(tmp+1);
      if(*tmp > max_v){
	max_i = i;
	max_j = j;
	max_v = *tmp;
      }
    }
  }
  // print_matrix(h, m, n, read, adap);
  // print_matrix(t, m, n, read, adap);
  // Traceback
  int _t, _l, _d, _align_n = 0, end_i = max_i, end_j = max_j;
  char _align[2][max_read_size];
  // 1 - Diag, 2->Right, 3->Down
  while(max_v != 0){
    _d = h[max_i-1][max_j-1];
    _l = h[max_i][max_j-1];
    _t = h[max_i-1][max_j];
    switch(t[max_i][max_j]){
    case 1:
      max_v = _d;
      max_i = max_i -1;
      max_j = max_j - 1;
      _align[0][_align_n] = read[max_i];
      _align[1][_align_n] = adap[max_j];
      _align_n++;
      break;
    case 2:
      max_v = _l;
      max_j = max_j - 1;
      _align[0][_align_n] = '-';
      _align[1][_align_n] = read[max_i];
      _align_n++;
      break;
    case 3:
      max_v = _t;
      max_i = max_i - 1;
      _align[0][_align_n] = adap[max_i];
      _align[1][_align_n] = '-';
      _align_n++;
      break;
    case 0:
      max_v = 0;
      break;
    }
  }
  print_alignment(_align, _align_n);
  return 0;
}

int main(){
  std::string read = "CAATGGTTTTGCTTTGGCCTGGTTGGCAATACGAGCGATGGCTGTTCCACGCACTGACAACATCACCTTGGCAATCCTAGCTGCTCTGACACCACTGGCCCGAGGCACACTGCTTGTAGCGTGGAGAGCAGGCCTTGCTACTTGTGGGGGGTTCATGCTCCTCTCTCTGAAGGGGAAAGGTAGTGTGAAGAAGAACCTACCATTTGTCATGGCCTTGGGACTAACCGCTGTGAGGCTGGTTGACCCCATC";
  std::string adap = "GCACCTTG";
  align_seqs(read, adap);
  return 0;
}
