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

  bool is_leaf_node(){
    return (nchildren == 0);
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
    if(this->begin == -1)
      return "R";
    return s.substr(this->begin, this->end - this->begin + 1);
  }

  void extend_path(int e){
    this->end = e;
  }

  suffix_node* add_child(int ext, int b, int e, suffix_node* p, suffix_node* l){
    suffix_node *n = new suffix_node(b, e, p, l);
    this->children[ext] = n;
    this->nchildren++;
    return n;
  }

  suffix_node* add_leaf(int ext_ind, int* str_ind, int beg){
    suffix_node *new_int_node;
    int ext = str_ind[ext_ind];
    this->add_child(ext, ext_ind, ext_ind, this, 0);
    this->add_child(str_ind[this->begin+beg], this->begin + beg, this->end, this, 0);
    this->end = this->begin + beg - 1;
    return this;
  }

  void print(std::string s){
    for(int i = 0; i < this->get_depth(); i++){
      std::cout << " ";
    }
    std::string t = (this->begin == -1) ? "R" : " "+s.substr(this->begin, this->end - this->begin + 1);
    std::cout << t;
    if(this->suffix_link != 0)
      std::cout << " --- " << this->suffix_link->get_path(s);
    std::cout << std::endl;
    for(int i = 0; i<alphabet.length();i++){
      if(this->children[i]!=0)
	this->children[i]->print(s);
    }
  }

  suffix_node* walk_path(int &beg, int* str_ind, int &s){
    suffix_node *node = this;
    while(s >= 0){
      std::cout << "Walk Start: " << node->begin << " Length: " << node->get_length() << " Suffix:" << s << std::endl;
      if(node->get_length() < s){
	std::cout << "Walk (Length < s): " << node->begin << std::endl;
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
  suffix_node *cur_node = root, *new_internal_node = 0, *new_node = 0;
  int str_ind[MAX_SIZE], n_str_ind = 0, ext, suffix_length = 0, beg = 0, end_length = 0;
  for(int i = 0; i < s.length();i++){
    str_ind[i] = alphabet.find(s[i]);
    n_str_ind++;
  }
  root->add_child(str_ind[0], 0, 0, root, root);
  for(int i = 0; i < n_str_ind - 1;i++){
    ext = str_ind[i+1];		// Index in alphabet
    new_node = 0;
    for(int j = 0; j<=i+1;j++){
      suffix_length = (i + 1 - j) + 1;
      beg = j;
      // Go to suffix link by traversing at most 1 edge
      std::cout << std::endl;
      std::cout << "Previous Current Node: " << cur_node->get_path(s) << std::endl;
      if(cur_node->suffix_link == 0){ // Since root has a suffix link to itself
	if(cur_node->parent->get_length() == 1)
	  cur_node = root;	// Single character nodes assumed to have suffix link to root.
	else
	  cur_node = cur_node->parent->suffix_link;
      } else {
	cur_node = cur_node->suffix_link;
      }
      std::cout << "To be inserted: " << s.substr(beg, suffix_length) << std::endl;
      std::cout << "Current Node: " << cur_node->get_path(s) << std::endl;
      cur_node = cur_node->walk_path(beg, str_ind, suffix_length);
      // if(cur_node->begin!=-1)
      // 	end_length = suffix_length + (beg - cur_node->begin);
      // else
      // 	end_length = suffix_length;
      end_length = suffix_length;
      if(cur_node->get_length() < end_length){
	if(cur_node->is_leaf_node()){	       // If its a leaf node Rule 1
	  std::cout << "Rule 1" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->get_path(s) << std::endl;
	  cur_node->extend_path(i+1); // Speedup Trick 2
	  // if(new_node != 0)
	  //   new_node->suffix_link = cur_node;
	  // new_node = 0;
	} else {		// Extension Rule 2
	  std::cout << "Extension Rule 2" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->get_path(s) << std::endl;
	  cur_node->add_child(ext, i+1, i+1, cur_node, 0);
	  if(new_node != 0)
	    new_node->suffix_link = cur_node;
	  new_node = cur_node;
	}
      } else if (cur_node->get_length() >= end_length){
	if(str_ind[cur_node->begin + suffix_length -1] == ext){ // Rule 3
	  std::cout << "Rule 3" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->get_path(s) << std::endl;
	  // if(new_node != 0)
	  //   new_node->suffix_link = cur_node;
	  // new_node = 0;
	} else{			// Rule 2
	  std::cout << "Rule 2" << std::endl;
	  std::cout << "Beg: " << beg << std::endl;
	  std::cout << "suffix_length: " << suffix_length << std::endl;
	  std::cout << "node path: " << cur_node->get_path(s) << std::endl;
	  cur_node->add_leaf(i+1, str_ind, suffix_length - 1);
	  if(new_node != 0)
	    new_node->suffix_link = cur_node;
	  new_node = cur_node;
	}
      }
      std::cout << "Tree: " << std::endl;
      root->print(s);
      // if(cur_node->get_length() == suffix_length){ // Rule 1
      // 	std::cout << "Rule 1" << std::endl;
      // 	cur_node->extend_path(i+1);
      // 	if(new_node != 0)
      // 	  new_node->suffix_link = cur_node;
      // 	new_node = new_internal_node = 0;
      // } else if (cur_node->get_length() > suffix_length){ // Rule 2
      // 	std::cout << "Rule 2" << std::endl;
      // 	if(str_ind[cur_node->begin + suffix_length-1] != ext){
      // 	  new_internal_node = cur_node->add_leaf(i+1, str_ind, suffix_length);
      // 	  if(new_node != 0)
      // 	    new_node->suffix_link = new_internal_node;
      // 	  new_node = new_internal_node;
      // 	} else {		// Rule 3
      // 	  std::cout << "Rule 3" << std::endl;
      // 	}
      // } else if(cur_node->get_length() < suffix_length - 1){	// Create new Node
      // 	std::cout << "Extension Rule 2" << std::endl;
      // 	cur_node->children[ext] = new suffix_node(i+1, i+1, cur_node, 0);
      // 	if(new_node != 0)
      // 	  new_node->suffix_link = cur_node;
      // 	new_node = cur_node;
      // 	std::cout << "Tree: " << std::endl;
      // 	root->print(s);
      // }
    }
  }
  return root;
}

int main(){
  std::string s = "GATAGACA-";
  suffix_node* root = build_suffix_tree(s);
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
