#ifndef INDEX_GROUP_H_
#define INDEX_GROUP_H_

#include "UID.h"

struct IndexGroup {
  typedef unsigned int Index;
  bool is_single;
  bool is_forward;
  unsigned int size;
  Index index;
  Index * indices;

  Index & operator[](unsigned int a);
  IndexGroup(unsigned int size);
  void set_sister(IndexGroup *);
};

#endif
