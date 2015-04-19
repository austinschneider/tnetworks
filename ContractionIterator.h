#ifndef CONTRACTION_ITERATOR_H_
#define CONTRACTION_ITERATOR_H_

#include <map>
#include <vector>

#include "Dimension.h"

struct ContractionIterator {
  typedef Dimension::ID ID;
  std::map<ID, unsigned int> indices_map;
  std::vector<Dimension> indices;
  unsigned int n_unique_indices;

  bool last_is_unique;
  bool is_done;

  ContractionIterator & operator++();

  ContractionIterator();
};

#endif
