#include "Tensor.h"
#include "IndexGroup.h"

unsigned int coordinate_transoform_nd_to_1d(Geometry & geometry, Coordinates & coord) {
  if(coord.size() != geometry.size())
    std::cout << "ERROR: Coordinate size does not match geometry!" << std::endl;
  unsigned int total = 0;
  unsigned int product = 1;
  //for(int i=0; i<geometry.size(); ++i) {
  for(int i=geometry.size() - 1; i>=0; --i) {
    if(coord[i] >= geometry[i].size)
      std::cout << "ERROR: Out of bounds coordinate!" << std::endl;
    total += product*coord[i];
    product *= geometry[i].size;
  }
  return total;
}

Coordinates coordinate_transoform_1d_to_nd(Geometry & geometry, unsigned int I) {
  Coordinates coord(geometry.size(), 0);
  unsigned int total = I;
  //for(int i=0; i<geometry.size(); ++i) {
  for(int i=geometry.size() - 1; i>=0; --i) {
    coord[i] = total % geometry[i].size;
    total /= geometry[i].size;
  }
  return coord;
}

unsigned int coordinate_transoform_nd_to_1d(Geometry & geometry, IndexGroup & ig, Coordinates & coord) {
  if(coord.size() != ig.size)
      std::cout << "ERROR: IndexGroup size does not match geometry!" << std::endl;
  if(geometry.size() < ig.size)
    std::cout << "ERROR: IndexGroup size is larger than geometry!" << std::endl;
  unsigned int total = 0;
  unsigned int product = 1;
  //for(int i=0; i<ig.size; ++i) {
  for(int i=ig.size - 1; i>=0; --i) {
    if(coord[ig[i]] >= geometry[ig[i]].size)
      std::cout << "ERROR: Out of bounds coordinate!" << std::endl;
    total += product*coord[ig[i]];
    product *= geometry[ig[i]].size;
  }
  return total;
}

Coordinates coordinate_transoform_1d_to_nd(Geometry & geometry, IndexGroup & ig, unsigned int I) {
  Coordinates coord(ig.size, 0);
  unsigned int total = I;
  //for(int i=0; i<ig.size; ++i) {
  for(int i=ig.size - 1; i>=0; --i) {
    coord[i] = total % geometry[ig[i]].size;
    total /= geometry[ig[i]].size;
  }
  return coord;
}

void coordinate_transform(Geometry & init_g, Geometry & final_g, std::vector<IndexGroup> & init_ig, Coordinates & init_coords, Coordinates & final_coords) {
  final_coords.clear();
  final_coords.resize(final_g.size());
  unsigned int f_iter = 0;
  for(int i=0; i<init_g.size(); ++i) {
    IndexGroup & ig = init_ig[i];
    if(ig.is_forward) {
      if(ig.is_single) {
        final_coords[*(ig.indices)] = init_coords[ig.index];
      }
      else {
        Coordinates temp_coords = coordinate_transoform_1d_to_nd(init_g, ig, init_coords[ig.index]);
        for(int j=0; j<temp_coords.size(); ++j) {
          final_coords[*(ig.indices)] = temp_coords[j];
        }
      }
    }
    else {
      if(ig.is_single) {
        final_coords[ig.index] = init_coords[*(ig.indices)];
      }
      else {
        final_coords[ig.index] = coordinate_transoform_nd_to_1d(init_g, ig, init_coords);
      }
    }
  }
}
