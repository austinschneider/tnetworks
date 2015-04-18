#include "Tensor.h"

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
