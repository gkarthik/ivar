#include "iostream"
#include "algorithm"

const int MAX_LENGTH = 50;

class suffix_node{
public:
  char a;
  suffix_node **children, *parent, *suffix_link;
  int nchildren;

  suffix_node(char n, suffix_node *p){
    this->children = new suffix_node*[MAX_LENGTH];
    this->parent = p;
    this->a = n;
    this->nchildren = 0;
  }

  int get_depth(){
    suffix_node *cur = this;
    int depth = 1;
    while(cur->parent){
      depth += 1;
      cur = cur->parent;
    }
    return depth;
  }

  std::string get_path(){
    suffix_node *cur = this;
    std::string s = "";
    while(cur->parent){
      s += cur->a;
      cur = cur->parent;
    }
    std::reverse(s.begin(), s.end());
    return s;
  }

  void create_path(std::string s){
    suffix_node *cur = this;
    for(unsigned int i =0;i<s.length();i++){
      cur = cur->add_child(s[i]);
    }
  }

  void print_tree(){
    for(int i = 0;i < this->get_depth();i++){
      std::cout << " ";
    }
    char n = (this->a == 0) ? '1' : this->a;
    std::cout << n << std::endl;
    for(int i = 0;i<this->nchildren;i++){
      for(int j = 0;j < this->children[i]->get_depth();j++){
	std::cout << " ";
      }
      this->children[i]->print_tree();
    }
  }

  suffix_node* add_child(char n){
    suffix_node *node = new suffix_node(n, this);
    this->children[nchildren] = node;
    this->nchildren++;
    return node;
  }

  suffix_node* find_end_node(std::string s){
    for(int i =0; i < this->nchildren; i++){
      if(s[0] == this->children[i]->a){
	return this->children[i]->find_end_node(s.substr(1, s.length() - 1));
      }
    }
    return this;
  }
};

class suffix_tree{
public:
  suffix_node *root;

  suffix_tree(std::string s){
    this->root = new suffix_node(0, 0);
    this->build_suffix_tree(s);
    this->print_suffix_tree();
  }

  void build_suffix_tree(std::string s){
    suffix_node *tmp;
    // this->root = new suffix_node(s[0], 0);
    std::string suf;
    for(unsigned int i = 0; i<s.length();i++){
      for(unsigned int j = 0;j<=i;j++){
	suf = s.substr(j, (i-j+1));
	tmp = this->root->find_end_node(suf);
	if(tmp->get_path().compare(suf) != 0)
	  tmp->add_child(s[i]);
      }
    }
  }

  void print_suffix_tree(){
    std::cout << "Tree:" << std::endl;
    this->root->print_tree();
  }
};

int main(int argc, char* argv[]){
  std::string s = "axabxb$";
  suffix_tree *suf_tree = new suffix_tree(s);
  return 0;
}
