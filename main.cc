
#include <stdlib.h>

#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "UID.h"
#include "Tensor.h"
#include "Dimension.h"
#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

int main() {

  srand (time(NULL));



  unsigned int i_s = 3;
  unsigned int j_s = 3;
  unsigned int k_s = 3;

  Dimension i(i_s, UID::get_id());
  Dimension j(j_s, UID::get_id());
  Dimension k(k_s, UID::get_id());

  std::vector<unsigned int> coord;
  coord.push_back(0); coord.push_back(0); coord.push_back(0);

  std::vector<Dimension> geo_2;
  geo_2.push_back(i);
  geo_2.push_back(j);
  geo_2.push_back(k);

  Tensor<double> t2(geo_2);
  t2.zero();

  double m2[j_s][i_s][k_s] = {
      {
          {0.9073, 0.7158, -0.3698},
          {0.8924, -0.4898, 2.4288},
          {2.1488, 0.3054, 2.3753}
      },
      {
          {1.7842, 1.6970, 0.0151},
          {1.7753, -1.5077, 4.0337},
          {4.2495, 0.3207, 4.7146},
      },
      {
          {2.1236, -0.0740, 1.4429},
          {-0.6631, 1.9103, -1.7495},
          {1.8260, 2.1335, -0.2716},
      }
  };

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        t2.at(coord) = m2[j][i][k];
      }
    }
  }
 std::cout << "Begin Element Printout:" << std::endl;

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        std::cout << "A[" << i << ", " << j << ", " << k << "] = " << t2.at(coord) << std::endl;
      }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;
  }
  std::cout << std::endl;

  //std::cout << unfold(t2, 0) << std::endl << std::endl;
  //std::cout << unfold(t2, 1) << std::endl << std::endl;
  //std::cout << unfold(t2, 2) << std::endl << std::endl;
  //Tensor<double> repeat = fold(unfold(t2, 0), t2.geometry, 0);
  //std::cout << unfold(repeat, 0) << std::endl;

  Tensor<double> S;
  Tensor<double> U[t2.geometry.size()];

  HOSVD(t2, S, U);

  std::cout << "U2:" << std::endl;
  std::vector<Tensor<double> > t_list;
  t_list.push_back(S);
  for(int i=0; i<t2.geometry.size(); ++i) {
    t_list.push_back(U[i]);
    std::cout << unfold(U[i], 0) << std::endl << std::endl;
  }

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        std::cout << "S[" << i << ", " << j << ", " << k << "] = " << S.at(coord) << std::endl;
      }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;
  }
  std::cout << std::endl;

  Tensor<double> A = Tensor<double>::contract(t_list);

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        std::cout << "A[" << i << ", " << j << ", " << k << "] = " << A.at(coord) << std::endl;
      }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Index id's:" << std::endl;
  for(int i=0; i<A.geometry.size(); ++i) {
    std::cout << A.geometry[i].id << " ";
  }
  std::cout << std::endl << std::endl;

  /*
  unsigned int i_s = 2;
  unsigned int j_s = 2;
  unsigned int k_s = 2;

  Dimension i(i_s, UID::get_id());
  Dimension j(j_s, UID::get_id());
  Dimension k(k_s, UID::get_id());

  std::vector<unsigned int> coord(3, 0);

  double m2[i_s][j_s][k_s] = {
      {
          {1, 2},
          {3, 4}
      },
      {
          {5, 6},
          {7, 8}

      }
  };

  std::vector<Dimension> geo_0;
  geo_0.push_back(i);
  geo_0.push_back(j);
  geo_0.push_back(k);

  Tensor<double> t0(geo_0);
  t0.zero();

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        t0.at(coord) = m2[i][j][k];
      }
    }
  }

  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        std::cout << t0.at(coord) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Index id's:" << std::endl;
  for(int i=0; i<t0.geometry.size(); ++i) {
    std::cout << t0.geometry[i].id << " ";
  }
  std::cout << std::endl << std::endl;

  Tensor<double> S;
  Tensor<double> U[t0.geometry.size()];

  HOSVD(t0, S, U);

  std::cout << "Index id's:" << std::endl;
  for(int i=0; i<S.geometry.size(); ++i) {
    std::cout << S.geometry[i].id << " ";
  }
  std::cout << std::endl << std::endl;

  std::cout << "U2:" << std::endl;
  std::vector<Tensor<double> > t_list;
  t_list.push_back(S);
  for(int i=0; i<t0.geometry.size(); ++i) {
    t_list.push_back(U[i]);
    std::cout << unfold(U[i], 0) << std::endl << std::endl;
    std::cout << "Index id's:" << std::endl;
    for(int j=0; j<U[i].geometry.size(); ++j) {
      std::cout << U[i].geometry[j].id << " ";
    }
    std::cout << std::endl << std::endl;
  }

  Tensor<double> A = Tensor<double>::contract(t_list);

  std::cout << "A:" << std::endl;
  for(int i=0; i<i_s; ++i) {
    coord[0] = i;
    for(int j=0; j<j_s; ++j) {
      coord[1] = j;
      for(int k=0; k<k_s; ++k) {
        coord[2] = k;
        std::cout << A.at(coord) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Index id's:" << std::endl;
  for(int i=0; i<A.geometry.size(); ++i) {
    std::cout << A.geometry[i].id << " ";
  }
  std::cout << std::endl << std::endl;

  return 0;
  */
}
