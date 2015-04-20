
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "UID.h"
#include "Tensor.h"
#include "Dimension.h"
#include "redsvd/redsvd.hpp"
#include "redsvd/util.hpp"

int main() {
    Dimension i(3, UID::get_id());
    Dimension j(3, UID::get_id());
    Dimension k(3, UID::get_id());

    std::vector<unsigned int> coord;
    coord.push_back(0); coord.push_back(0);

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
    Tensor<double> repeat = fold(unfold(t2, 0), t2.geometry, 0);
    std::cout << unfold(repeat, 0) << std::endl;

    Tensor<double> S;
    Tensor<double> U[t2.geometry.size()];

    HOSVD(t2, S, U);

    return 0;
}
