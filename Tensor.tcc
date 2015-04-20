#ifndef TENSOR_TCC_
#define TENSOR_TCC_

#include "Tensor.h"
#include "Utilities.h"
#include "Dimension.h"
#include "ContractionIterator.h"

#include <map>
#include <set>
#include <cstring>
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
T & Tensor<T>::at(Coordinates & coord) {
  unsigned int I = coordinate_transoform_nd_to_1d(geometry, coord);
  Coordinates coord_result = coordinate_transoform_1d_to_nd(geometry, I);
  if(coord_result.size() != coord.size())
    std::cout << "ERROR: Coordinate size mismatch!" << std::endl;
  for(int i=0; i<coord_result.size(); ++i) {
    if(coord_result[i] != coord[i])
    std::cout << "ERROR: Coordinate mismatch!" << std::endl;
  }
  return elements[coordinate_transoform_nd_to_1d(geometry, coord)];
}

template <class T>
T & Tensor<T>::at(ContractionIterator & it) {
  unsigned int total = 0;
  unsigned int product = 1;
  for(int i=0; i<geometry.size(); ++i) {
    std::map<Dimension::ID, unsigned int>::iterator m_it = it.indices_map.find(geometry[i].id);
    if(m_it == it.indices_map.end())
      std::cout << "ERROR: ContractionIterator does not contain index!" << std::endl;
    else {
      if(m_it->second >= geometry[i].size)
        std::cout << "ERROR: Out of bounds coordinate!" << std::endl;
      total += product*m_it->second;
      product *= geometry[i].size;
    }
  }
  return elements[total];
}

template <class T>
void Tensor<T>::zero() {
  std::memset(elements, 0, sizeof(elements));
}

template <class T>
Tensor<T>::Tensor(Geometry & geo) {
  geometry = geo;
  elements = new T[this->size()];
}

template <class T>
Tensor<T> Tensor<T>::contract(std::vector<Tensor<T> > & tensors) {
  typedef Dimension::ID ID;
  std::map<ID, unsigned int> original_indices;

  std::set<ID> unique_indices_set;
  std::set<ID> contracted_indices_set;

  for(int t_i=0; t_i<tensors.size(); ++t_i) {
    for(int dim_i=0; dim_i<tensors[t_i].geometry.size(); ++dim_i) {
      Dimension & dim = tensors[t_i].geometry[dim_i];
      if(original_indices.find(dim.id) != original_indices.end()) {
        if(dim.size != original_indices[dim.id])
          std::cout << "ERROR: Contracted index sizes do not match!" << std::endl;
        unique_indices_set.erase(dim.id);
        contracted_indices_set.insert(dim.id);
      }
      else{
        original_indices[dim.id] = dim.size;
        unique_indices_set.insert(dim.id);
      }
    }
  }

  ContractionIterator c_it;
  c_it.n_unique_indices = unique_indices_set.size();
  for(std::set<ID>::iterator it = unique_indices_set.begin(); it != unique_indices_set.end(); ++it) {
    c_it.indices_map[*it] = 0;
    Dimension d(original_indices[*it], *it);
    c_it.indices.push_back(d);
  }

  Tensor<T> result(c_it.indices);

  for(std::set<ID>::iterator it = contracted_indices_set.begin(); it != contracted_indices_set.end(); ++it) {
    c_it.indices_map[*it] = 0;
    c_it.indices.push_back(Dimension(original_indices[*it], *it));
  }

  T * total = NULL;

  for(; c_it.is_done == false; ++c_it) {
    if(c_it.last_is_unique) {
      total = & result.at(c_it);
    }

    T product = 1;
    for(int t_i=0; t_i<tensors.size(); ++t_i) {
      product *= tensors[t_i].at(c_it);
    }
    (*total) += product;
  }
  return result;
}

#endif
