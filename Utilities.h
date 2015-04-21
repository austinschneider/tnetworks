#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <set>
#include <vector>
#include <iostream>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

#include "UID.h"
#include "Tensor.h"
#include "IndexGroup.h"

unsigned int coordinate_transoform_nd_to_1d(Geometry &, Coordinates &);

Coordinates coordinate_transoform_1d_to_nd(Geometry &, unsigned int);

unsigned int coordinate_transoform_nd_to_1d(Geometry &, IndexGroup &, Coordinates &);

Coordinates coordinate_transoform_1d_to_nd(Geometry &, IndexGroup &, unsigned int);

void coordinate_transform(Geometry &, Geometry &, std::vector<IndexGroup> &, Coordinates &, Coordinates &);


// Unfinished
template <class T>
void reduce_rank(Tensor<T> & t, std::vector<IndexGroup> & reductions) {
  std::set<unsigned int> unchanges_indices;

  for(int i=0; i<t.geometry.size(); ++i) {
    unchanges_indices.insert(i);
  }
  for(int i=0; i<reductions.size(); ++i) {
    for(int j=0; j<reductions[i].size; ++j) {
      unchanges_indices.erase(i);
    }
  }

  std::vector<IndexGroup> igs;

  for(std::set<unsigned int>::iterator it = unchanges_indices.begin(); it != unchanges_indices.end(); ++it) {
    IndexGroup ig(1);
    ig.is_forward = true;
    ig.index = *it;
    (&ig.indices) = igs.size();
    igs.push_back(ig);
  }

  for(int i=0; i<reductions.size(); ++i) {
    IndexGroup & ig = reductions[i];
    ig.is_forward = false;

  }

  Tensor<T> result;
}

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
Tensor<T> fold(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix, Geometry & geometry, unsigned int dim) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  unsigned int dim_size = matrix.rows();
  unsigned int m_dim_size = matrix.cols();

  Tensor<T> t(geometry);

  Geometry g;
  for(int i=0; i<t.geometry.size(); ++i) {
    if(i != dim)
      g.push_back(t.geometry[i]);
  }

  Coordinates coords_nd[m_dim_size];
  for(int j=0; j<m_dim_size; ++j) {
    coords_nd[j] = coordinate_transoform_1d_to_nd(g, j);
    coords_nd[j].insert(coords_nd[j].begin()+dim, j);
  }

  for(int i=0; i<dim_size; ++i) {
    for(int j=0; j<m_dim_size; ++j) {
      coords_nd[j][dim] = i;
      t.at(coords_nd[j]) = matrix(i, j);
    }
  }

  return t;
}



template <class T>
void HOSVD(Tensor<T> & t, Tensor<T> & S, Tensor<T> * U) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  M transforms[t.geometry.size()];
  for(int i=0; i<t.geometry.size(); ++i) {
    M A = unfold(t, i);
    transforms[i] = REDSVD::RedSVD<T>(A, A.cols()).matrixU();
  }

  unsigned int axis = 0;
  M t_unfolded = unfold(t, axis);
  M B(1,1);
  B(0,0) = 1;
  for(int i=0; i<t.geometry.size(); ++i) {
    if(i != axis)
      B = kroneckerProduct(B,transforms[i]).eval();
  }

  M S_unfolded = transforms[axis].adjoint() * t_unfolded * B;

  S = fold(S_unfolded, t.geometry, axis);

  Geometry matrix_geo;
  matrix_geo.push_back(Dimension(0, 0));
  matrix_geo.push_back(Dimension(0, 0));
  for(int i=0; i<t.geometry.size(); ++i) {
    UID::ID id = UID::get_id();
    S.geometry[i].id = id;
    matrix_geo[0].size = t.geometry[i].size;
    matrix_geo[1].size = t.geometry[i].size;
    matrix_geo[0].id = t.geometry[i].id;
    matrix_geo[1].id = id;
    U[i] = fold(transforms[i], matrix_geo, 0);
    std::cout << transforms[i] << std::endl << std::endl;
  }


};

#endif
