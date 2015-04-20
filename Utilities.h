#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <iostream>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

#include "Tensor.h"

unsigned int coordinate_transoform_nd_to_1d(Geometry &, Coordinates &);

Coordinates coordinate_transoform_1d_to_nd(Geometry &, unsigned int);

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> unfold(Tensor<T> & t, unsigned int dim) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  unsigned int dim_size = t.geometry[dim].size;
  unsigned int m_dim_size = t.size() / dim_size;

  Geometry g;
  for(int i=0; i<t.geometry.size(); ++i) {
    if(i != dim)
      g.push_back(t.geometry[i]);
  }

  M result(dim_size, m_dim_size);

  Coordinates coords_nd[m_dim_size];
  for(int j=0; j<m_dim_size; ++j) {
    coords_nd[j] = coordinate_transoform_1d_to_nd(g, j);
    coords_nd[j].insert(coords_nd[j].begin()+dim, j);
  }

  for(int i=0; i<dim_size; ++i) {
    for(int j=0; j<m_dim_size; ++j) {
      coords_nd[j][dim] = i;
      result(i,j) = t.at(coords_nd[j]);
    }
  }

  return result;
};

template <class T>
void HOSVD(Tensor<T> & t) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  REDSVD::RedSVD<T> svds[t.geometry.size()];
  for(int i=0; i<t.geometry.size(); ++i) {
    M A = unfold(t, i);
    svds[i] = REDSVD::RedSVD<T>(A, A.cols());
  }
};

#endif
