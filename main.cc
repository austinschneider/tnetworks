
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "Tensor.h"
#include "Dimension.h"
#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

int main() {
    Dimension i(3, 0);
    Dimension j(3, 1);
    Dimension k(3, 2);
    /*
    std::vector<Dimension> geo_0;
    geo_0.push_back(i);
    geo_0.push_back(j);

    std::vector<Dimension> geo_1;
    geo_1.push_back(j);
    geo_1.push_back(k);

    Tensor<double> t0(geo_0);
    t0.zero();
    Tensor<double> t1(geo_1);
    t1.zero();

    double m0[3][3] = {
      {7, 5, 2},
      {3, 6, 12},
      {4, 0, 4}
    };

    double m1[3][3] = {
      {8, 2, 3},
      {9, 19, 6},
      {5, 7, 3}
    };
    */
    std::vector<unsigned int> coord;
    coord.push_back(0); coord.push_back(0);
    /*
    for(int i=0; i<3; ++i) {
      coord[0] = i;
      for(int j=0; j<3; ++j) {
        coord[1] = j;
        t0.at(coord) = m0[i][j];
        t1.at(coord) = m1[i][j];
      }
    }

    std::vector<Tensor<double> > ts;
    ts.push_back(t0);
    ts.push_back(t1);

    Tensor<double> t_res = Tensor<double>::contract(ts);

    for(int i=0; i<3; ++i) {
      coord[0] = i;
      for(int j=0; j<3; ++j) {
        coord[1] = j;
        std::cout << t_res.at(coord) << " ";
      }
      std::cout << std::endl;
    }

    Eigen::MatrixXf A = Eigen::MatrixXf::Random(10, 20);
    REDSVD::RedSVD REDSVDOfA(A, 20);
    */

    std::vector<Dimension> geo_2;
    geo_2.push_back(i);
    geo_2.push_back(j);
    geo_2.push_back(k);

    Tensor<double> t2(geo_2);
    t2.zero();

    coord.push_back(0);
    unsigned int ii = 0;
    for(int i=0; i<3; ++i) {
      coord[0] = i;
      for(int j=0; j<3; ++j) {
        coord[1] = j;
        for(int k=0; k<3; ++k) {
          coord[2] = k;
          t2.at(coord) = ii;
          ++ii;
        }
      }
    }

    std::cout << unfold(t2, 0) << std::endl;
    //HOSVD(t2);

    return 0;
}
