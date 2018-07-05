#include "iostream"
#include "algorithm"

const std::string alphabet = "XATGCRYSWKMBDHVN-";
const int MAX_CHAR = alphabet.length();
const int MAX_SIZE = 300;

class suffix_node{
public:
  int begin, end, nchildren;
  suffix_node **children, *parent, *suffix_link;

  suffix_node(int b, int e, suffix_node *p, suffix_node *l){
    this->children = new suffix_node*[MAX_CHAR];
    for(int i = 0; i < MAX_CHAR;i++)
      this->children[i] = 0;
    this->begin = b;
    this->end = e;
    this->nchildren = 0;
    this->parent = p;
    this->suffix_link =l;
  }

  bool contains_indice(int ind){
    return (ind >= this->begin && ind <= this->end);
  }

  int get_length(){
    if(this->end == -1)
      return 0;			// For root
    return this->end - this->begin + 1;
  }

  int get_depth(){
    int ctr = 0;
    suffix_node* n = this;
    while(n->begin!=-1){
      n = n->parent;
      ctr++;
    }
    return ctr;
  }

  std::string get_path(std::string s){
    return s.substr(this->begin, this->end - this->begin + 1);
  }

  void extend_path(int e){
    this->end = e;
  }

  suffix_node* add_leaf(int ext_ind, int* str_ind, int suf_length){
    int ext = str_ind[ext_ind];
    suffix_node* new_int_node;
    this->children[ext] = new suffix_node(ext_ind,ext_ind, this, 0);
    new_int_node = new suffix_node(this->begin + suf_length, this->end, this, 0);
    this->children[str_ind[this->begin + suf_length]] = new_int_node;
    this->end = this->begin+suf_length - 1;
    return new_int_node;
  }

  void print(std::string s){
    for(int i = 0; i < this->get_depth(); i++){
      std::cout << " ";
    }
    std::string t = (this->begin == -1) ? "R" : " "+s.substr(this->begin, this->end - this->begin + 1);
    std::cout << t  << std::endl;
    for(int i = 0; i<alphabet.length();i++){
      if(this->children[i]!=0)
	this->children[i]->print(s);
    }
  }

  suffix_node* walk_path(int beg, int* str_ind, int &s){
    suffix_node* node = this;
    if(node->get_length()<s){
      std::cout << "Node Path: " << alphabet[str_ind[beg+node->get_length()]] << std::endl;
      while(s >= 0 && node->children[str_ind[beg+node->get_length()]] != 0){
	node = node->children[str_ind[beg+node->get_length()]];
	s -= node->parent->get_length();
	beg += node->parent->get_length();
	std::cout << alphabet[str_ind[beg+node->get_length()]] << std::endl;
      }
    }
    return node;
  }
};

suffix_node* build_suffix_tree(std::string s){
  suffix_node *root = new suffix_node(-1,-1,0,0);
  root->suffix_link = root;
  suffix_node *cur_node = root, *new_internal_node = 0, *new_node = 0;
  int str_ind[MAX_SIZE], n_str_ind = 0, ext, suffix_length = 0, ctr = 0;
  for(int i = 0; i < s.length();i++){
    str_ind[i] = alphabet.find(s[i]);
    n_str_ind++;
  }
  root->children[str_ind[0]] = new suffix_node(0,0,root,root);
  for(int i = 0; i < n_str_ind - 1;i++){
    ext = str_ind[i+1];		// Index in alphabet
    for(int j = 0; j<=i+1;j++){
      suffix_length = i - j + 1;
      // Go to suffix link by traversing at most 1 edge
      if(cur_node->suffix_link == 0) // Since root has a suffix link to itself
	cur_node = cur_node->parent->suffix_link;
      else
	cur_node = cur_node->suffix_link;
      if(suffix_length > cur_node->get_length()){
	while(suffix_length > 0 && !cur_node->contains_indice(i)){
	  if(cur_node->children[str_ind[cur_node->end + 1]] == 0)
	    break;
	  cur_node = cur_node->children[str_ind[cur_node->end + 1]];
	  suffix_length -= cur_node->parent->get_length();
	}
      }
      // if(cur_node->get_length() == suffix_length){ // Rule 1
      // 	cur_node->extend_path(i+1);
      // 	if(new_node != 0)
      // 	  new_node->suffix_link = cur_node;
      // 	new_node = new_internal_node = 0;
      // } else if (cur_node->get_length() > suffix_length){ // Rule 2
      // 	if(alphabet[cur_node->begin + (cur_node->end - suffix_length)] != ext){
      // 	  new_internal_node = cur_node->add_leaf(i+1, str_ind, suffix_length);
      // 	  if(new_node != 0)
      // 	    new_node->suffix_link = new_internal_node;
      // 	  new_node = new_internal_node;
      // 	}
      // } else if(cur_node->get_length() < suffix_length){	// Create new Node
      // 	new_internal_node = new suffix_node(i+1, i+1, cur_node, 0);
      // 	cur_node->children[ext] = new_internal_node;
      // 	if(new_node != 0)
      // 	  new_node->suffix_link = new_internal_node;
      //   new_node = new_internal_node;
      // }
    }
  }
  return root;
}

int main(){
  std::string s = "XABC";
  // suffix_node* root = build_suffix_tree(s);
  // std::cout << "Tree: " << std::endl;
  // root->print(s);
  int str_ind[10];
  for(int i = 0;i<s.length();i++){
    str_ind[i] = alphabet.find(s[i]);
  }
  suffix_node* root = new suffix_node(-1,-1,0,0);
  root->children[str_ind[0]] = new suffix_node(0,2,root,root);
  root->children[str_ind[0]]->add_leaf(3,str_ind,1);
  root->print(s);
  int size = 3;
  suffix_node* node = root->walk_path(0, str_ind, size);
  std::cout << "Size: " << size << " " << node->get_path(s) << std::endl;
  return 0;
}
