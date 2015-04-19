#include "ContractionIterator.h"

ContractionIterator & ContractionIterator::operator++() {
  for(int i=indices.size()-1; i>=0; --i) {
    Dimension & dim = indices[i];
    indices_map[dim.id] += 1;
    if(indices_map[dim.id] >= dim.size) {
      indices_map[dim.id] = 0;
      if(i==0)
        is_done = true;
    }
    else {
      break;
      last_is_unique = i < n_unique_indices;
    }
  }
  return *this;
};

ContractionIterator::ContractionIterator(): n_unique_indices(0), last_is_unique(true), is_done(false) {};
