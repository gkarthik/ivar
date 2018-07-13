#include "iostream"
#include "algorithm"

const std::string alphabet = "ATGCRYSWKMBDHVN-X#@";
const int MAX_CHAR = alphabet.length();
const int MAX_SIZE = 300;

class suffix_node{
public:
  int begin, nchildren, *end;
  suffix_node **children, *parent, *suffix_link;

  suffix_node(int b, int *e, suffix_node *p, suffix_node *l){
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
    if(this->end == 0)
      return 0;			// For root
    return *(this->end) - this->begin + 1;
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
    return s.substr(this->begin, *(this->end) - this->begin + 1);
  }

  void extend_path(int *e){
    this->end = e;
  }

  suffix_node* add_child(int ext, int b, int *e, suffix_node* l){
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

  bool contains_depth(int depth){
    return this->get_depth() <= depth && depth <= (this->get_length() + this->get_depth());
  }

  void print(std::string s){
    for(int i = 0; i < this->get_depth(); i++){
      std::cout << " ";
    }
    std::string t = (this->begin == -1) ? "R" : " "+s.substr(this->begin, *(this->end) - this->begin + 1);
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

  bool walk_next(int &beg, int &suffix_length){
    suffix_node *node = this;
    if(suffix_length >= node->get_length()){
      beg += node->get_length();
      suffix_length -= node->get_length();
      return true;
    }
    return false;
  }
};

suffix_node* build_suffix_tree(std::string s){
  suffix_node *root = new suffix_node(-1,0,0,0);
  root->suffix_link = root;
  root->parent = root;
  suffix_node *cur_node = root, *new_node = 0, *new_cur_node = 0;
  int str_ind[MAX_SIZE], n_str_ind = 0, suffix_length = 0, j =0, *leaf_end=  new int, suffix_count = 0, beg = -1, *end;
  for(int i = 0; i < s.length();i++){
    str_ind[i] = alphabet.find(s[i]);
    n_str_ind++;
  }
  *leaf_end = 0;
  // root->add_child(str_ind[0], 0, leaf_end, root);
  for(int i = 0; i < n_str_ind;i++){
    *leaf_end = i;		// Handle Extension Rule 1
    suffix_count++;
    new_node = 0;
    while(suffix_count > 0){
      if(suffix_length == 0)
	beg = i;
      root->print(s);
      std::cout << "Cur Node: " << cur_node->get_path(s) << std::endl;
      if(cur_node->children[str_ind[beg]] == 0){
	std::cout << "First Rule 2" << std::endl;
	cur_node->add_child(str_ind[beg], beg, leaf_end, 0); // Trick 3
	if(new_node != 0){
	  new_node->suffix_link = cur_node;
	  new_node = 0;		// No internal node created in Rule 2 here.
	}
      } else {
	new_cur_node = cur_node->children[str_ind[beg]];
	if(new_cur_node->walk_next(beg, suffix_length)){
	  std::cout << "Continue!" << std::endl;
	  cur_node = new_cur_node;
	  continue;
	}
	if(str_ind[new_cur_node->begin + suffix_length] == str_ind[i]){
	  std::cout << "Rule 3" << std::endl;
	  if(new_node != 0 && cur_node->end != 0){
	    new_node->suffix_link = cur_node;
	    new_node = 0;		// No internal node created in Rule 2 here.
	  }
	  suffix_length++;
	  break;
	}
	std::cout << "Rule 2" << std::endl;
	suffix_node *new_int_node;
	end = new int;
	*end = new_cur_node->begin + suffix_length - 1;
	//New internal node
	new_int_node = cur_node->add_child(str_ind[beg], new_cur_node->begin, end, 0);
	new_int_node->add_child(str_ind[i], i, leaf_end, 0);
	new_cur_node->begin += suffix_length;
	new_int_node->add_child(new_cur_node, str_ind[new_cur_node->begin]);
	if(new_node != 0)
	  new_node->suffix_link = new_int_node;
	new_node = new_int_node;
      }
      suffix_count--;
      if(cur_node->end==0 && suffix_length > 0){
	suffix_length--;
	beg = i - suffix_count + 1;
      } else if(cur_node->end!=0){
	if(cur_node->suffix_link != 0){
	  cur_node = cur_node->suffix_link;
	} else {
	  cur_node = cur_node->parent->suffix_link;
	}
      } else if(cur_node->get_length() == 1){
	cur_node = root;
      }
      std::cout << "Tree:" << std::endl;
      std::cout << std::endl;
    }
  }
  return root;
}

int main(){
  std::cout << alphabet.length() << " " << MAX_CHAR << std::endl;
  std::string s = "XABXA#BABXBA@";
  suffix_node *root = build_suffix_tree(s);
  std::cout << "Tree: " << std::endl;
  root->print(s);
  return 0;
}
