/*
 *  Copyright (c) 2010 Daisuke Okanohara
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

#ifndef REDSVD_HPP__
#define REDSVD_HPP__

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <vector>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "util.hpp"

namespace REDSVD {

template <class T>
class RedSVD {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
  RedSVD(){}

  template <class Mat>
  RedSVD(Mat& A){
    int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
    run(A, r);
  }

  template <class Mat>
  RedSVD(Mat& A, const int rank){
    run(A, rank);
  }

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
    r = (r < A.rows()) ? r : A.rows();

    // Gaussian Random Matrix for A^T
    Matrix O(A.rows(), r);
    Util<T>::sampleGaussianMat(O);

    // Compute Sample Matrix of A^T
    Matrix Y = A.transpose() * O;

    // Orthonormalize Y
    Util<T>::processGramSchmidt(Y);

    // Range(B) = Range(A^T)
    Matrix B = A * Y;

    // Gaussian Random Matrix
    Matrix P(B.cols(), r);
    Util<T>::sampleGaussianMat(P);

    // Compute Sample Matrix of B
    Matrix Z = B * P;

    // Orthonormalize Z
    Util<T>::processGramSchmidt(Z);

    // Range(C) = Range(B)
    Matrix C = Z.transpose() * B;

    Eigen::JacobiSVD<Matrix> svdOfC(C, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // C = USV^T
    // A = Z * U * S * V^T * Y^T()
    matU_ = Z * svdOfC.matrixU();
    matS_ = svdOfC.singularValues();
    matV_ = Y * svdOfC.matrixV();
  }

  const Matrix& matrixU() const {
    return matU_;
  }

  const Vector& singularValues() const {
    return matS_;
  }

  const Matrix& matrixV() const {
    return matV_;
  }

private:
  Matrix matU_;
  Vector matS_;
  Matrix matV_;
};

template <class T>
class RedSymEigen {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
  RedSymEigen(){}

  template <class Mat>
  RedSymEigen(Mat& A, const int rank){
    run(A, rank);
  }

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
    r = (r < A.rows()) ? r : A.rows();

    // Gaussian Random Matrix
    Matrix O(A.rows(), r);
    Util<T>::sampleGaussianMat(O);

    // Compute Sample Matrix of A
    Matrix Y = A.transpose() * O;

    // Orthonormalize Y
    Util<T>::processGramSchmidt(Y);

    Matrix B = Y.transpose() * A * Y;
    Eigen::SelfAdjointEigenSolver<Matrix> eigenOfB(B);

    eigenValues_ = eigenOfB.eigenvalues();
    eigenVectors_ = Y * eigenOfB.eigenvectors();
  }

  const Matrix& eigenVectors() const {
    return eigenVectors_;
  }

  const Vector& eigenValues() const {
    return eigenValues_;
  }

private:
  Vector eigenValues_;
  Matrix eigenVectors_;
};

template <class T>
class RedPCA {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
public:
  RedPCA(){}

  template <class Mat>
  RedPCA(const Mat& A, const int rank) {
    run(A, rank);
  }

  template <class Mat>
  void run(const Mat& A, const int rank) {
    RedSVD<T> redsvd;
    redsvd.run(A, rank);
    const Vector& S = redsvd.singularValues();
    principalComponents_ = redsvd.matrixV();
    scores_              = redsvd.matrixU() * S.asDiagonal();
  }

  const Matrix& principalComponents() const {
    return principalComponents_;
  }

  const Matrix& scores() const {
    return scores_;
  }

 private:
  Matrix principalComponents_;
  Matrix scores_;
};

}

#endif // REDSVD_HPP__
