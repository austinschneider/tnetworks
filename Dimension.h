#ifndef DIMENSION_H_
#define DIMENSION_H_

struct Dimension {
  typedef unsigned long int ID;
  unsigned int size;
  ID id;
  Dimension();
  Dimension(unsigned int, ID);
};

#endif
