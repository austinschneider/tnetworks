#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <set>
#ifndef MY_CMATH_
#define MY_CMATH_
#include <cmath>
#endif
#include <vector>
#include <iostream>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/KroneckerProduct>

//#include "redsvd/redsvd.hpp"
//#include "redsvd/util.hpp"
#include <eigen3/Eigen/SVD>

#include "UID.h"
#include "Tensor.h"
#include "IndexGroup.h"

inline void gdb() {
  return;
}

unsigned int coordinate_transoform_nd_to_1d(Geometry &, Coordinates &);

Coordinates coordinate_transoform_1d_to_nd(Geometry &, unsigned int);

unsigned int coordinate_transoform_nd_to_1d(Geometry &, std::map<unsigned int, unsigned int> &, IndexGroup &, Coordinates &);

Coordinates coordinate_transoform_1d_to_nd(Geometry &, std::map<unsigned int, unsigned int> &, IndexGroup &, unsigned int);

void coordinate_transform(Geometry &, Geometry &, std::vector<IndexGroup> &, Coordinates &, Coordinates &);


template <class T>
void reduce_rank(Tensor<T> & t, std::vector<IndexGroup> & reductions);

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> unfold(Tensor<T> & t, unsigned int dim);

template <class T>
Tensor<T> fold(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix, Geometry & geometry, unsigned int dim);

template <class T>
void get_HOSVD_U(Tensor<T> & t, Tensor<T> & U, unsigned int index_pos);

template <class T>
void HOSVD(Tensor<T> & t, Tensor<T> & S, Tensor<T> * U);

template <class T>
T calculate_cut_norm(Tensor<T> & t, unsigned int index_pos, unsigned int D_cut);

template <class T>
void truncate(Tensor<T> & t, unsigned int D_cut);

template <class T>
void coordinate_transform(Tensor<T> & init_t, Tensor<T> & final_t, std::vector<IndexGroup> & init_ig, Coordinates & init_coords, Coordinates & final_coords);

#include "Utilities.tcc"

#endif
