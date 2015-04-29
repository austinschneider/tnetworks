#include "Tensor.h"
#include "IndexGroup.h"

unsigned int coordinate_transoform_nd_to_1d(Geometry & geometry, Coordinates & coord) {
  if(coord.size() != geometry.size()) {
    std::cout << "ERROR: Coordinate size does not match geometry!" << std::endl;
    std::cout << "coord.size(): " << coord.size() << std::endl;
    std::cout << "geometry.size(): " << geometry.size() << std::endl;
  }
  unsigned int total = 0;
  unsigned int product = 1;
  //for(int i=0; i<geometry.size(); ++i) {
  for(int i=geometry.size() - 1; i>=0; --i) {
    if(coord[i] >= geometry[i].size) {
      std::cout << "ERROR: Out of bounds coordinate 12!" << std::endl;
      std::cout << "\ti: " << i << " coord[i]: " << coord[i] << " geometry[i].size: " << geometry[i].size << std::endl;
    }
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

unsigned int coordinate_transoform_nd_to_1d(Geometry & geometry, std::map<unsigned int, unsigned int> & index_map, IndexGroup & ig, Coordinates & coord) {
  /*std::cout << "coordinate_transoform_nd_to_1d index_map:" << std::endl;
  for(std::map<unsigned int, unsigned int>::iterator it = index_map.begin(); it != index_map.end(); ++it) {
    std::cout << "ID: " << it->first << " = " << it->second << std::endl;
  }
  std::cout << "coord: ";
  for(int i=0; i<coord.size(); ++i) {
    std::cout << coord[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "geometry: ";
  for(int i=0; i<geometry.size(); ++i) {
    std::cout << geometry[i].size << " ";
  }
  std::cout << std::endl;*/
  if(geometry.size() < ig.size)
    std::cout << "ERROR: IndexGroup size is larger than geometry!" << std::endl;
  if(index_map.size() < ig.size)
    std::cout << "ERROR: IndexGroup size is larger than index map!" << std::endl;
  unsigned int total = 0;
  unsigned int product = 1;
  //for(int i=0; i<ig.size; ++i) {
  for(int i=ig.size - 1; i>=0; --i) {
    if(index_map.find(ig[i]) == index_map.end()) {
      std::cout << "ERROR: Index map does not contain specified index!" << std::endl;
      std::cout << "\tIndex failure: i = " << i << " ig[i] = " << ig[i] << std::endl;
      std::cout << std::endl;
      //gdb();
    }
    if(coord[index_map[ig[i]]] >= geometry[index_map[ig[i]]].size) {
      std::cout << "ERROR: Out of bounds coordinate 62!" << std::endl;
      std::cout << "\ti: " << i << " ig[i]: " << ig[i] <<  " index_map[ig[i]]: " << index_map[ig[i]] << " coord[]: " << coord[index_map[ig[i]]] << " geometry[].size: " << geometry[index_map[ig[i]]].size << std::endl;
    }
    total += product*coord[index_map[ig[i]]];
    product *= geometry[index_map[ig[i]]].size;
  }
  return total;
}

Coordinates coordinate_transoform_1d_to_nd(Geometry & geometry, std::map<unsigned int, unsigned int> & index_map, IndexGroup & ig, unsigned int I) {
  if(geometry.size() < ig.size)
    std::cout << "ERROR: IndexGroup size is larger than geometry!" << std::endl;
  if(index_map.size() < ig.size)
    std::cout << "ERROR: IndexGroup size is larger than index map!" << std::endl;
  Coordinates coord(ig.size, 0);
  unsigned int total = I;
  //for(int i=0; i<ig.size; ++i) {
  for(int i=ig.size - 1; i>=0; --i) {
    if(index_map.find(ig[i]) == index_map.end())
      std::cout << "ERROR: Index map does not contain specified index!" << std::endl;
    coord[i] = total % geometry[ig[i]].size;
    total /= geometry[ig[i]].size;
  }
  return coord;
}
