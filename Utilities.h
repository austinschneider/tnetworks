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
void reduce_rank(Tensor<T> & t, std::vector<IndexGroup> & reductions) {
  std::set<unsigned int> unchanged_indices;

  for(int i=0; i<t.geometry.size(); ++i) {
    unchanged_indices.insert(t.geometry[i].id);
  }
  for(int i=0; i<reductions.size(); ++i) {
    for(int j=0; j<reductions[i].size; ++j) {
      unchanged_indices.erase(reductions[i][j]);
    }
  }

  std::vector<IndexGroup> igs;
  std::map<UID::ID, unsigned int> map_ids_to_igs;

  for(std::set<unsigned int>::iterator it = unchanged_indices.begin(); it != unchanged_indices.end(); ++it) {
    IndexGroup ig(1);
    ig.is_forward = true;
    ig.index = *it;
    (*ig.indices) = *it;
    map_ids_to_igs[ig.index] = igs.size();
    igs.push_back(ig);
  }

  for(int i=0; i<reductions.size(); ++i) {
    IndexGroup ig = reductions[i].copy();
    ig.is_forward = false;
    IndexGroup::Index min;
    for(int r_i=0; r_i<ig.size; ++r_i) {
      if(ig[r_i] < min || r_i == 0) {
        min = ig[r_i];
      }
    }
    ig.index = min;
    map_ids_to_igs[ig.index] = igs.size();
    igs.push_back(ig);
  }

  std::map<unsigned int, unsigned int> & init_map = t.map();

  Geometry final_g;
  for(int i=0; i<t.geometry.size(); ++i) {
    if(map_ids_to_igs.find(t.geometry[i].id) != map_ids_to_igs.end()) {
      unsigned size = 1;
      IndexGroup & ig = igs[map_ids_to_igs[t.geometry[i].id]];
      for(int j=0; j<ig.size; ++j) {
        size *= t.geometry[init_map[ig[j]]].size;
      }
      final_g.push_back(Dimension(size,t.geometry[i].id));
    }
  }

  Tensor<T> result(final_g);

  Coordinates init_coords(t.geometry.size(), 0);
  Coordinates final_coords(result.geometry.size(), 0);

  // Check is for NULL pointers
  /*for(int i=0; i<igs.size(); ++i) {
    if(igs[i].indices == NULL)
      std::cout << "ERROR: igs[" << i << "] has NULL pointer" << std::endl;
    std::cout << "IndexGroup: " << igs[i].index;
    if(igs[i].is_forward)
      std::cout << " ->";
    else
      std::cout << " <-";
    for(int j=0; j<igs[i].size; ++j)
      std::cout << " " << igs[i][j];
    std::cout << std::endl;
  }*/

  for(int it=0; it<t.size(); ++it) {
    //for(int c=0; c<init_coords.size(); ++c) {
    //  std::cout << init_coords[c] << " ";
    //}
    //std::cout << std::endl;
    coordinate_transform(t, result, igs, init_coords, final_coords);

    //for(int c=0; c<final_g.size(); ++c) {
    //  std::cout << final_g[c].size << " ";
    //}
    //std::cout << std::endl;
    //for(int c=0; c<final_coords.size(); ++c) {
    //  std::cout << final_coords[c] << " ";
    //}
    //std::cout << std::endl;
    result.at(final_coords) = t.at(init_coords);
    // Increment init_coords
    for(int i=init_coords.size()-1; i>=0; --i) {
      Dimension & dim = t.geometry[i];
      init_coords[i] += 1;
      if(init_coords[i] >= dim.size) {
        init_coords[i] = 0;
        if(i==0)
          break;
      }
      else {
        break;
      }
    }
  }

  for(int i=0; i<igs.size(); ++i) {
    igs[i].destroy();
  }

  t.destroy();
  t = result;
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> unfold(Tensor<T> & t, unsigned int dim) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  unsigned int dim_size = t.geometry[dim].size;
  unsigned int m_dim_size = t.size() / dim_size;

  Geometry g(t.geometry.size() - 1, Dimension());
  int ii=0;
  for(int i=0; i<t.geometry.size(); ++i) {
    if(i != dim) {
      g[ii] = t.geometry[i];
      ++ii;
    }
  }

  M result(dim_size, m_dim_size);

  Coordinates coords_nd[m_dim_size];
  for(int j=0; j<m_dim_size; ++j) {
    coords_nd[j] = coordinate_transoform_1d_to_nd(g, j);
    coords_nd[j].insert(coords_nd[j].begin()+dim, 0);
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
void get_HOSVD_U(Tensor<T> & t, Tensor<T> & U, unsigned int index_pos) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  M A = unfold(t, index_pos);
  Geometry matrix_geo;
  matrix_geo[0].size = t.geometry[index_pos].size;
  matrix_geo[1].size = t.geometry[index_pos].size;
  matrix_geo[0].id = 0;
  matrix_geo[1].id = 0;
  M m_U = Eigen::JacobiSVD<M>(A, Eigen::ComputeFullU).matrixU();
  U = fold(m_U, matrix_geo, 0);
}

template <class T>
void HOSVD(Tensor<T> & t, Tensor<T> & S, Tensor<T> * U) {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M;
  M transforms[t.geometry.size()];
  //std::cout << "Unfolded:" << std::endl;
  for(int i=0; i<t.geometry.size(); ++i) {
    M A = unfold(t, i);
    //std::cout << A << std::endl << std::endl;
    //transforms[i] = REDSVD::RedSVD<T>(A, A.cols()).matrixU().eval();
    transforms[i] = Eigen::JacobiSVD<M>(A, Eigen::ComputeFullU).matrixU();
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

  //std::cout << "S Unfolded:" << std::endl;
  //std::cout << S_unfolded << std::endl << std::endl;

  S = fold(S_unfolded, t.geometry, axis);

  //std::cout << "U:" << std::endl;
  Geometry matrix_geo;
  matrix_geo.push_back(Dimension(0, 0));
  matrix_geo.push_back(Dimension(0, 0));
  for(int i=0; i<t.geometry.size(); ++i) {
    UID::ID id = UID::get_id();
    //std::cout << "New ID: " << id << std::endl;
    S.geometry[i].id = id;
    matrix_geo[0].size = t.geometry[i].size;
    matrix_geo[1].size = t.geometry[i].size;
    matrix_geo[0].id = t.geometry[i].id;
    matrix_geo[1].id = id;
    //std::cout << "matrix_geo[0].id: " << matrix_geo[0].id << std::endl;
    //std::cout << "matrix_geo[1].id: " << matrix_geo[1].id << std::endl;
    U[i] = fold(transforms[i], matrix_geo, 0);
    //std::cout << "U[i].geometry[0].id: " << U[i].geometry[0].id << std::endl;
    //std::cout << "U[i].geometry[1].id: " << U[i].geometry[1].id << std::endl;
    //std::cout << transforms[i] << std::endl << std::endl;
  }

  /*
  std::vector<unsigned int> coord(S.geometry.size(), 0);

  std::cout << "S:" << std::endl;
  for(int i=0; i<S.geometry[0].size; ++i) {
    coord[0] = i;
    for(int j=0; j<S.geometry[1].size; ++j) {
      coord[1] = j;
      for(int k=0; k<S.geometry[2].size; ++k) {
        coord[2] = k;
        std::cout << S.at(coord) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  */
};

template <class T>
T calculate_cut_norm(Tensor<T> & t, unsigned int index_pos, unsigned int D_cut) {
  typedef Dimension::ID ID;
  if(t.geometry[index_pos].size <= D_cut)
    return 0;
  ContractionIterator c_it;
  c_it.n_unique_indices = 1;
  c_it.indices_map[t.geometry[index_pos].id] = D_cut;
  c_it.indices.push_back(t.geometry[index_pos]);

  for(int i=0; i<t.geometry.size(); ++i) {
    if(i != index_pos) {
      c_it.indices_map[t.geometry[i].id] = 0;
      c_it.indices.push_back(t.geometry[i]);
    }
  }

  T total = 0;

  for(; c_it.is_done == false; ++c_it) {
    total += std::pow(std::abs(t.at(c_it)), 2.0);
  }

  return total;
}

template <class T>
void truncate(Tensor<T> & t, unsigned int D_cut) {
  typedef Dimension::ID ID;
  std::map<unsigned int, unsigned int> sizes; // map index size to position in truncation sets
  std::vector<std::set<unsigned int> *> truncation_sets; // sets of positions to truncate
  for(int i=0; i<t.geometry.size(); ++i) {
    if(t.geometry[i].size > D_cut) {
      unsigned int size = t.geometry[i].size;
      if(sizes.find(size) == sizes.end()) {
        sizes[size] = truncation_sets.size();
        truncation_sets.push_back(new std::set<unsigned int>());
      }
      int map_size = sizes.size();
      if(map_size == 0)
        gdb();
      int accessor = sizes[size];
      truncation_sets[accessor]->insert(i);
    }
  }

  if(truncation_sets.size() == 0)
    return;

  Tensor<T> * current = &t;
  Tensor<T> * next = NULL;
  for(int i=0; i<truncation_sets.size(); ++i) {
    std::set<unsigned int> & t_set = *truncation_sets[i];
    Tensor<T> S;
    Tensor<T> U[current->geometry.size()];
    HOSVD(*current, S, U);
    S.destroy();
    unsigned int min_norm_pos = 0;
    T min_norm = -1;

    for(std::set<unsigned int>::iterator it = t_set.begin(); it != t_set.end(); ++it) {
      T norm = calculate_cut_norm(*current, *it, D_cut);
      if(norm < min_norm | min_norm < 0) {
        min_norm_pos = *it;
        min_norm = norm;
      }
    }

    unsigned int orig_size = current->geometry[min_norm_pos].size;

    Geometry trunc_geo;
    trunc_geo.push_back(Dimension(orig_size , 0));
    trunc_geo.push_back(Dimension(D_cut, 0));

    Tensor<T> original_truncation_matrix(trunc_geo);
    std::vector<unsigned int> coords(2, 0);
    for(int j=0; j<orig_size; ++j) {
      coords[0] = j;
      for(int k=0; k<D_cut; ++k) {
        coords[1] = k;
        original_truncation_matrix.at(coords) = U[min_norm_pos].at(coords);
      }
    }

    for(int U_i=0; U_i<current->geometry.size(); ++U_i) {
      U[U_i].destroy();
    }

    Tensor<T> truncation_matrices[t_set.size()];
    std::vector<Tensor<T> > contraction_tensors;
    contraction_tensors.push_back(*current);
    for(std::set<unsigned int>::iterator it = t_set.begin(); it != t_set.end(); ++it) {
      truncation_matrices[i] = original_truncation_matrix;
      ID new_id = UID::get_id();
      truncation_matrices[i].geometry[0].id = new_id;
      truncation_matrices[i].geometry[1].id = current->geometry[i].id;
      current->geometry[i].id = new_id;
      contraction_tensors.push_back(truncation_matrices[i]);
    }

    next = new Tensor<T>();
    *next = Tensor<T>::contract(contraction_tensors);
    current->destroy();
    if(i > 0)
      delete current;
    original_truncation_matrix.destroy();
    current = next;
    next = NULL;
  }

  for(int i=0; i<truncation_sets.size(); ++i) {
    delete truncation_sets[i];
  }

  t = *current;
}

template <class T>
void coordinate_transform(Tensor<T> & init_t, Tensor<T> & final_t, std::vector<IndexGroup> & init_ig, Coordinates & init_coords, Coordinates & final_coords) {
  Geometry & init_g = init_t.geometry;
  Geometry & final_g = final_t.geometry;
  std::map<unsigned int, unsigned int> & init_map = init_t.map();
  std::map<unsigned int, unsigned int> & final_map = final_t.map();

  /*std::cout << "init_map:" << std::endl;
  for(std::map<unsigned int, unsigned int>::iterator it = init_map.begin(); it != init_map.end(); ++it) {
    std::cout << "ID: " << it->first << " = " << it->second << std::endl;
  }

  std::cout << "final_map:" << std::endl;
  for(std::map<unsigned int, unsigned int>::iterator it = final_map.begin(); it != final_map.end(); ++it) {
    std::cout << "ID: " << it->first << " = " << it->second << std::endl;
  }

  // Check is for NULL pointers
  std::cout << "init_ig:" << std::endl;
  for(int i=0; i<init_ig.size(); ++i) {
    if(init_ig[i].indices == NULL)
      std::cout << "ERROR: init_ig[" << i << "] has NULL pointer" << std::endl;
    std::cout << "IndexGroup: " << init_ig[i].index;
    if(init_ig[i].is_forward)
      std::cout << " ->";
    else
      std::cout << " <-";
    for(int j=0; j<init_ig[i].size; ++j)
      std::cout << " " << init_ig[i][j];
    std::cout << std::endl;
  }*/

  final_coords.clear();
  final_coords.resize(final_g.size());
  for(int i=0; i<init_ig.size(); ++i) {
    IndexGroup & ig = init_ig[i];
    if(ig.is_forward) {
      if(ig.is_single) {
        if(init_map.find(ig.index) == init_map.end())
          std::cout << "ERROR: init_map does not contain index id: " << ig.index << std::endl;
        final_coords[final_map[ig[0]]] = init_coords[init_map[ig.index]];
      }
      else {
        if(init_map.find(ig.index) == init_map.end())
          std::cout << "ERROR: init_map does not contain index id: " << ig.index << std::endl;
        Coordinates temp_coords = coordinate_transoform_1d_to_nd(init_g, init_map, ig, init_coords[init_map[ig.index]]);
        for(int j=0; j<temp_coords.size(); ++j) {
          final_coords[final_map[ig[0]]] = temp_coords[j];
        }
      }
    }
    else {
      if(ig.is_single) {
        if(init_map.find(ig.index) == init_map.end())
          std::cout << "ERROR: init_map does not contain index id: " << ig.index << std::endl;
        final_coords[final_map[ig.index]] = init_coords[init_map[ig[0]]];
      }
      else {
        final_coords[final_map[ig.index]] = coordinate_transoform_nd_to_1d(init_g, init_map, ig, init_coords);
      }
    }
  }
}

#endif
