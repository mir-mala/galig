//=================================
// include guard
#ifndef _MEMSLIST_HPP_
#define _MEMSLIST_HPP_

//=================================
// forward declared dependencies
//class foo;

//=================================
// included dependencies
#include <string>
#include <vector>
#include <forward_list>

struct Mem {
    int t;
    int p;
    int l;

    Mem(int t, int p, int l) {
	this->t = t;
	this->p = p;
	this->l = l;
    }

    std::string toStr() {
	return "(" + std::to_string(t) + "," +
	    std::to_string(p) + "," +
	    std::to_string(l) + ")";
    }
};

class MemsList {
private:
    int length;
    std::vector<std::forward_list<Mem> > mems;
public:
    MemsList(const int& l);
    void addMem(const int& t, const int& p, const int& l);
    std::forward_list<Mem> getMems(const int& i);
    int getLength();
};

#endif