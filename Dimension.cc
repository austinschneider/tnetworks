#include "Dimension.h"

Dimension::Dimension() {};

Dimension::Dimension(unsigned int s, ID i): size(s), id(i) {};

Dimension::Dimension(const Dimension & other) {
  size = other.size;
  id = other.id;
}

Dimension::~Dimension() {

}
