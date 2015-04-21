
#include <iostream>

#include "UID.h"
#include "IndexGroup.h"

IndexGroup::Index & IndexGroup::operator[](unsigned int a) {
  return indices[a];
}

IndexGroup::IndexGroup(unsigned int s): size(s), is_single(s <= 1) {
  if(is_single)
    indices = new Index;
  else
    indices = new Index[size];
}
