#include "Tensor.h"
#include "Dimension.h"

#include <map>
#include <set>
#include <iostream>

template <class T>
unsigned int Tensor<T>::size() {
    unsigned int s = 1;
    for(int i=0; i<geometry.size(); ++i) {
        s *= geometry[i].size;
    }
    return s;
}

template <class T>
T & Tensor<T>::at(Tensor<T>::Coordinates & coord) {
    if(coord.size() != geometry.size())
        std::cout << "ERROR: Coordinate size does not match geometry!" << std::endl;
    unsigned int total = 0;
    unsigned int product = 1;
    for(int i=0; i<geometry.size(); ++i) {
        if(coord[i] >= geometry[i].size)
            std::cout << "ERROR: Out of bounds coordinate!" << std::endl;
        total += product*coord[i];
        product *= geometry[i].size;
    }
    return elements[total];
}

template <class T>
Tensor<T>::Tensor(Geometry & geo) {
    geometry = geo;
    elements = new T[this->size()];
}

template <class T>
Tensor<T> Tensor<T>::contract(std::vector<Tensor<T> > & tensors) {
    typedef Dimension::ID ID;
    std::map<ID, unsigned int *> original_indices;
    
    std::set<ID> unique_indices_set;
    std::set<ID> contracted_indices_set;

    for(int t_i=0; t_i<tensors.size(); ++t_i) {
        for(int dim_i=0; dim_i<tensors[t_i].geometry.size(); ++dim_i) {
            Dimension & dim = tensors[t_i].geometry[dim_i];
            if(original_indices.find(dim.id) != original_indices.end()) {
                unique_indices_set.erase(dim.id);
                contracted_indices_set.insert(dim.id);
            }
            else{
                original_indices[dim.id] = new unsigned int;
                unique_indices_set.insert(dim.id);
            }
        }
    }

}
