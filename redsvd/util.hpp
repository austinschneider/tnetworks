/*
 *  Copyright (c) 2011 Daisuke Okanohara
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#ifndef REDSVD_UTIL_HPP__
#define REDSVD_UTIL_HPP__

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <vector>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace REDSVD {

template <class T>
class Util{
public:
  typedef Eigen::SparseMatrix<T, Eigen::RowMajor> SMatrixXf;
  typedef std::vector<std::pair<int, T> > fv_t;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;

  static const T SVD_EPS = 0.0001f;

  static void convertFV2Mat(const std::vector<fv_t>& fvs, SMatrixXf& A);
  static void sampleGaussianMat(Matrix& x);
  static void processGramSchmidt(Matrix& mat);
  static double getSec();

private:
  static void sampleTwoGaussian(T& f1, T& f2);
};

}

#include "util.tpp"

#endif // REDSVD_UTIL_HPP_
