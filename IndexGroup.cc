
#include <iostream>

#include "UID.h"
#include "IndexGroup.h"

IndexGroup::Index & IndexGroup::operator[](unsigned int a) {
  return indices[a];
}

IndexGroup::IndexGroup(unsigned int s): size(s), is_single(s <= 1) {
    indices = new Index[size];
}

IndexGroup::IndexGroup(): size(0), is_single(size <= 1) {
  indices = NULL;
}

IndexGroup::IndexGroup(const IndexGroup & other) {
  is_single = other.is_single;
  is_forward = other.is_forward;
  size = other.size;
  index = other.index;
  indices = other.indices;
}

IndexGroup::~IndexGroup() {

}

void IndexGroup::destroy() {
  delete[] indices;
}

IndexGroup IndexGroup::copy() {
  IndexGroup copy = *this;
  copy.indices = new Index[size];
  for(int i=0; i<size; ++i) {
    copy[i] = indices[i];
  }
  return copy;
}
