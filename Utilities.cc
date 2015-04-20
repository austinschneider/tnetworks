#include "Tensor.h"

unsigned int coordinate_transoform_nd_to_1d(Geometry & geometry, Coordinates & coord) {
  if(coord.size() != geometry.size())
    std::cout << "ERROR: Coordinate size does not match geometry!" << std::endl;
  unsigned int total = 0;
  unsigned int product = 1;
  for(int i=0; i<geometry.size(); ++i) {
    if(coord[i] >= geometry[i].size)
      std::cout << "ERROR: Out of bounds coordinate!" << std::endl;
    total += product*coord[i];
    product *= geometry[i].size;
  }
  return total;
}

Coordinates coordinate_transoform_1d_to_nd(Geometry & geometry, unsigned int I) {
  Coordinates coord;
  unsigned int total = I;
  for(int i=0; i<geometry.size(); ++i) {
    coord.push_back(total % geometry[i].size);
    total /= geometry[i].size;
  }
  return coord;
}
