
#include <stdlib.h>

#ifndef MY_CMATH_
#define MY_CMATH_
#include <cmath>
#endif
#include <iostream>
#include <sstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "UID.h"
#include "Tensor.h"
#include "Dimension.h"
#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

int main(int argc, char * argv[]) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  // Define parameters
  typedef double T;
  T temp = 0.2;
  unsigned int n_iterations = 2;
  unsigned int D_cut = 24;

  if(argc > 1) {
    std::stringstream ss(argv[1]);
    ss >> temp;
  }
  if(argc > 2) {
    std::stringstream ss(argv[2]);
    ss >> n_iterations;
  }

  // Define IDs
  UID::ID common = UID::get_id();
  UID::ID up = UID::get_id();
  UID::ID down = UID::get_id();
  UID::ID left = UID::get_id();
  UID::ID right = UID::get_id();

  // Define W matrix geometry
  Geometry matrix_geo(2);
  matrix_geo[0].size = 2;
  matrix_geo[1].size = 2;
  matrix_geo[0].id = common;
  matrix_geo[1].id = 0;

  // Define W matrix
  Tensor<T> W(matrix_geo);
  Coordinates coords(2, 0);

  coords[0] = 0; coords[1] = 0; W.at(coords) = std::sqrt(std::cosh(1.0/temp));
  coords[0] = 0; coords[1] = 1; W.at(coords) = std::sqrt(std::sinh(1.0/temp));
  coords[0] = 1; coords[1] = 0; W.at(coords) = std::sqrt(std::cosh(1.0/temp));
  coords[0] = 1; coords[1] = 1; W.at(coords) = -std::sqrt(std::sinh(1.0/temp));

  /*
  coords[0] = 0; coords[1] = 0; W.at(coords) = -std::sinh(2.0/temp)/(2.0*(1.0+std::cosh(2.0/temp)));
  coords[0] = 0; coords[1] = 1; W.at(coords) = std::sinh(2.0/temp)/(2.0*(1.0-std::cosh(2.0/temp)));
  coords[0] = 1; coords[1] = 0; W.at(coords) = -std::sinh(2.0/temp)/(2.0*(1.0+std::cosh(2.0/temp)));
  coords[0] = 1; coords[1] = 1; W.at(coords) = -std::sinh(2.0/temp)/(2.0*(1.0-std::cosh(2.0/temp)));
  */

  // Define W matrices
  Tensor<T> W_up = W; W_up.geometry[1].id = up;
  Tensor<T> W_down = W; W_down.geometry[1].id = down;
  Tensor<T> W_left = W; W_left.geometry[1].id = left;
  Tensor<T> W_right = W; W_right.geometry[1].id = right;

  // Define contraction vector
  std::vector<Tensor<T> > contraction_tensors(4);
  contraction_tensors[0] = W_up;
  contraction_tensors[1] = W_down;
  contraction_tensors[2] = W_left;
  contraction_tensors[3] = W_right;

  // Construct original tensor for square geometry
  Tensor<T> start = Tensor<T>::contract(contraction_tensors);

  Tensor<T> * current = new Tensor<T>;
  *current = start.copy();
  Tensor<T> * next = NULL;

  UID::ID up_p = UID::get_id();
  UID::ID down_p = UID::get_id();
  UID::ID left_p = UID::get_id();
  UID::ID right_p = UID::get_id();

  IndexGroup up_group(2);
  up_group[0] = up;
  up_group[1] = up_p;

  IndexGroup down_group(2);
  down_group[0] = down;
  down_group[1] = down_p;

  std::vector<IndexGroup> post_h_reduction(2);
  post_h_reduction[0] = up_group;
  post_h_reduction[1] = down_group;

  IndexGroup left_group(2);
  left_group[0] = left;
  left_group[1] = left_p;

  IndexGroup right_group(2);
  right_group[0] = right;
  right_group[1] = right_p;

  std::vector<IndexGroup> post_v_reduction(2);
  post_v_reduction[0] = left_group;
  post_v_reduction[1] = right_group;

  UID::ID index_mapping[4];
  UID::ID start_index_mapping[4];

  start_index_mapping[0] = start.map()[up];
  start_index_mapping[1] = start.map()[down];
  start_index_mapping[2] = start.map()[left];
  start_index_mapping[3] = start.map()[right];

  for(int i=0; i<n_iterations-1; ++i) {
    Tensor<T> temp;
    temp = *current;

    index_mapping[0] = current->map()[up];
    index_mapping[1] = current->map()[down];
    index_mapping[2] = current->map()[left];
    index_mapping[3] = current->map()[right];

    start.geometry[start_index_mapping[0]].id = up;
    start.geometry[start_index_mapping[1]].id = common;
    start.geometry[start_index_mapping[2]].id = left_p;
    start.geometry[start_index_mapping[3]].id = right_p;

    std::cout << "up    " << current->geometry[index_mapping[0]].size << std::endl;
    std::cout << "down  " << current->geometry[index_mapping[1]].size << std::endl;
    std::cout << "left  " << current->geometry[index_mapping[2]].size << std::endl;
    std::cout << "right " << current->geometry[index_mapping[3]].size << std::endl;

    current->geometry[index_mapping[0]].id = common;
    current->geometry[index_mapping[1]].id = down;
    current->geometry[index_mapping[2]].id = left;
    current->geometry[index_mapping[3]].id = right;

    std::vector<Tensor<T> > contraction_tensors(2);
    contraction_tensors[0] = *current;
    contraction_tensors[1] = start;
    *current = Tensor<T>::contract(contraction_tensors);

    temp.destroy();

    //for(int in=0; in<temp.geometry.size(); ++in) {
    //  std::cout << "Index " << in << " id = " << temp.geometry[in].id << std::endl;
    //}

    reduce_rank(*current, post_v_reduction);
    std::cout << "up    " << current->geometry[index_mapping[0]].size << std::endl;
    std::cout << "down  " << current->geometry[index_mapping[1]].size << std::endl;
    std::cout << "left  " << current->geometry[index_mapping[2]].size << std::endl;
    std::cout << "right " << current->geometry[index_mapping[3]].size << std::endl;
    truncate(*current, D_cut);

    //for(int in=0; in<current->geometry.size(); ++in) {
    //  std::cout << "Index " << in << " id = " << current->geometry[in].id << std::endl;
    //}
  }

  start.destroy();
  start = current->copy();

  std::cout << "start.geometry.size() = " << start.geometry.size() << std::endl;

  start_index_mapping[0] = start.map()[up];
  start_index_mapping[1] = start.map()[down];
  start_index_mapping[2] = start.map()[left];
  start_index_mapping[3] = start.map()[right];

  std::cout << "start.geometry.size() = " << start.geometry.size() << std::endl;

  for(int i=0; i<n_iterations-1; ++i) {
    Tensor<T> temp;
    temp = *current;

    index_mapping[0] = current->map()[up];
    index_mapping[1] = current->map()[down];
    index_mapping[2] = current->map()[left];
    index_mapping[3] = current->map()[right];

    start.geometry[start_index_mapping[0]].id = up_p;
    start.geometry[start_index_mapping[1]].id = down_p;
    start.geometry[start_index_mapping[2]].id = common;
    start.geometry[start_index_mapping[3]].id = right;

    std::cout << "up    " << current->geometry[index_mapping[0]].size << std::endl;
    std::cout << "down  " << current->geometry[index_mapping[1]].size << std::endl;
    std::cout << "left  " << current->geometry[index_mapping[2]].size << std::endl;
    std::cout << "right " << current->geometry[index_mapping[3]].size << std::endl;

    current->geometry[index_mapping[0]].id = up;
    current->geometry[index_mapping[1]].id = down;
    current->geometry[index_mapping[2]].id = left;
    current->geometry[index_mapping[3]].id = common;

    std::vector<Tensor<T> > contraction_tensors(2);
    contraction_tensors[0] = *current;
    contraction_tensors[1] = start;
    *current = Tensor<T>::contract(contraction_tensors);

    temp.destroy();

    reduce_rank(*current, post_h_reduction);
    std::cout << "up    " << current->geometry[index_mapping[0]].size << std::endl;
    std::cout << "down  " << current->geometry[index_mapping[1]].size << std::endl;
    std::cout << "left  " << current->geometry[index_mapping[2]].size << std::endl;
    std::cout << "right " << current->geometry[index_mapping[3]].size << std::endl;
    truncate(*current, D_cut);

    /*for(int in=0; in<current->geometry.size(); ++in) {
      std::cout << "Index " << in << " id = " << current->geometry[in].id << std::endl;
    }*/

  }

  index_mapping[0] = current->map()[up];
  index_mapping[1] = current->map()[down];
  index_mapping[2] = current->map()[left];
  index_mapping[3] = current->map()[right];

  current->geometry[index_mapping[0]].id = up;
  current->geometry[index_mapping[1]].id = up;
  current->geometry[index_mapping[2]].id = left;
  current->geometry[index_mapping[3]].id = left;

  contraction_tensors.clear();
  contraction_tensors.push_back(*current);
  Tensor<T> result = Tensor<T>::contract(contraction_tensors);

  //std::cout << "result geometry size: " << result.geometry.size() << std::endl;
  std::cout << "result value: " << result.elements[0] << std::endl;

}












