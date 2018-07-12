#include "iostream"
#include "algorithm"

const std::string alphabet = "ATGCRYSWKMBDHVN-X#@";
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

  bool is_leaf_node(){
    return (nchildren == 0);
  }

  bool contains_child(int ext){
    return !(this->children[ext] == 0);
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
    if(this->begin == -1)
      return "R";
    return s.substr(this->begin, this->end - this->begin + 1);
  }

  void extend_path(int e){
    this->end = e;
  }

  suffix_node* add_child(int ext, int b, int e, suffix_node* l){
    suffix_node *n = new suffix_node(b, e, this, l);
    this->children[ext] = n;
    this->nchildren++;
    return n;
  }

  void add_child(suffix_node* c, int ext){
    this->children[ext] = c;
    c->parent = this;
    this->nchildren++;
  }

  suffix_node* get_child(int ext){
    return this->children[ext];
  }

  suffix_node* add_internal_node(int ext_ind, int* str_ind, int beg){
    suffix_node *new_int_node;
    int ext = str_ind[ext_ind];
    // Add new Internal Node
    new_int_node = this->parent->add_child(str_ind[this->begin], this->begin, this->begin + beg - 1, 0);
    //Add new node as child
    new_int_node->add_child(ext, ext_ind, ext_ind, 0);
    //Add older node as child
    this->begin = this->begin+beg;
    new_int_node->add_child(this, str_ind[this->begin]);
    // this->add_child(str_ind[this->begin+beg], this->begin + beg, this->end, this, 0);
    // this->end = this->begin + beg - 1;
    return new_int_node;
  }

  void print(std::string s){
    for(int i = 0; i < this->get_depth(); i++){
      std::cout << " ";
    }
    std::string t = (this->begin == -1) ? "R" : " "+s.substr(this->begin, this->end - this->begin + 1);
    std::cout << t;
    if(this->suffix_link != 0)
      std::cout << " --- ( " << this->suffix_link->parent->get_path(s) << " " << this->suffix_link->get_path(s) << ")";
    if(this->begin!=-1)
      std::cout << " - " << this->parent->get_path(s);
    std::cout << std::endl;
    for(int i = 0; i<alphabet.length();i++){
      if(this->children[i]!=0)
	this->children[i]->print(s);
    }
  }

  suffix_node* walk_path(int &beg, int* str_ind, int &s){
    suffix_node *node = this;
    while(s >= 0){
      if(node->get_length() <= s){
	if(node->children[str_ind[beg+node->get_length()]] == 0)
	  break;
	node = node->children[str_ind[beg+node->get_length()]];
	s -= node->parent->get_length();
	beg += node->parent->get_length();
      } else {
	return node;
      }
      // std::cout << alphabet[str_ind[beg+node->get_length()]] << std::endl;
    }
    return node;
  }
};

suffix_node* build_suffix_tree(std::string s){
  suffix_node *root = new suffix_node(-1,-1,0,0);
  root->suffix_link = root;
  root->parent = root;
  suffix_node *cur_node = root, *new_internal_node = 0, *new_node = 0;
  int str_ind[MAX_SIZE], n_str_ind = 0, ext, suffix_length = 0, beg = 0, end_length = 0;
  for(int i = 0; i < s.length();i++){
    str_ind[i] = alphabet.find(s[i]);
    n_str_ind++;
  }
  root->add_child(str_ind[0], 0, 0, root);
  for(int i = 0; i < n_str_ind - 1;i++){
    ext = str_ind[i+1];		// Index in alphabet
    new_node = 0;
    for(int j = 0; j<=i+1;j++){
      beg = j;
      std::cout << std::endl;
      std::cout << "To be inserted: " << s.substr(beg, suffix_length) << std::endl;
      // Go to suffix link by traversing at most 1 edge
      if(cur_node->suffix_link == 0){ // Since root has a suffix link to itself
	if(cur_node->parent->get_length() == 1)
	  cur_node = root;	// Single character nodes assumed to have suffix link to root.
	else
	  cur_node = cur_node->parent->suffix_link;
      } else {
	cur_node = cur_node->suffix_link;
      }
      beg += cur_node->get_depth();
      suffix_length = (i + 1 - beg) + 1;
      cur_node = cur_node->walk_path(beg, str_ind, suffix_length);
      end_length = suffix_length;
      if(cur_node->get_length() < end_length){
	if(cur_node->is_leaf_node()){	       // If its a leaf node Rule 1
	  std::cout << "Rule 1" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->parent->get_path(s) << " - " << cur_node->get_path(s) << std::endl;
	  cur_node->extend_path(i+1); // Speedup Trick 2
	  // if(new_node != 0)
	  //   new_node->suffix_link = cur_node;
	  // new_node = 0;
	} else {		// Extension Rule 2
	  std::cout << "Extension Rule 2" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->parent->get_path(s) << " - " << cur_node->get_path(s) << std::endl;
	  cur_node->add_child(ext, i+1, i+1, 0);
	  if(new_node != 0)
	    new_node->suffix_link = cur_node;
	  new_node = cur_node;
	}
      } else if (cur_node->get_length() >= end_length){
	if(str_ind[cur_node->begin + suffix_length -1] == ext){ // Rule 3
	  std::cout << "Rule 3" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->parent->get_path(s) << " - " << cur_node->get_path(s) << std::endl;
	  break;
	  // if(new_node != 0)
	  //   new_node->suffix_link = cur_node;
	  // new_node = 0;
	} else{			// Rule 2
	  std::cout << "Rule 2" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->parent->get_path(s) << " - " << cur_node->get_path(s) << std::endl;
	  cur_node = cur_node->add_internal_node(i+1, str_ind, suffix_length - 1);
	  if(new_node != 0)
	    new_node->suffix_link = cur_node;
	  new_node = cur_node; // New Internal Created
	  // std::cout << "New Created Node: " << std::endl;
	  // cur_node->print(s);
	  // std::cout << std::endl;
	}
      }
      std::cout << "Tree:" << std::endl;
      root->print(s);
      std::cout << std::endl;
    }
  }
  return root;
}

int main(){
  std::cout << alphabet.length() << " " << MAX_CHAR << std::endl;
  std::string s = "XABXA@BABXBA#";
  suffix_node *root = build_suffix_tree(s);
  std::cout << "Tree: " << std::endl;
  root->print(s);
  // int str_ind[10];
  // for(int i = 0;i<s.length();i++){
  //   str_ind[i] = alphabet.find(s[i]);
  // }
  // suffix_node* root = new suffix_node(-1,-1,0,0);
  // root->children[str_ind[1]] = new suffix_node(1,4,root,root);
  // root->print(s);
  // root->children[str_ind[1]]->add_leaf(1,str_ind,1);
  // root->print(s);
  // int size = 1;
  // int beg = 4;
  // suffix_node* node = root->walk_path(beg, str_ind, size);
  // std::cout << "Size: " << size << " " << node->get_path(s) << std::endl;
  // size = 2;
  // beg = 0;
  // node = root->walk_path(beg, str_ind, size);
  // std::cout << "Size: " << size << " " << node->get_path(s) << std::endl;
  // size = 1;
  // beg = 2;
  // node = root->walk_path(beg, str_ind, size);
  // std::cout << "Size: " << size << " " << node->get_path(s) << std::endl;
  // beg = 1;
  // size = 1;
  // node = root->walk_path(beg, str_ind, size);
  // std::cout << "Size: " << size << " " << node->get_path(s) << std::endl;
  return 0;
}
