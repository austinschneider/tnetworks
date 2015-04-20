#ifndef INDEX_GROUP_H_
#define INDEX_GROUP_H_

struct IndexGroup {
  typedef unsigned int Index;
  bool is_single;
  unsigned int size;
  Index * indices;
  IndexGroup * sister;

  Index & operator[](unsigned int a);
  IndexGroup(unsigned int size);
};

#endif
