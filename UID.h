#ifndef UID_H_
#define UID_H_

class ID_store {
public:
   typedef unsigned long int ID;
private:
   ID id;
public:
  ID_store() {
    id = 0;
  }
  ID get_id() {
    return id++;
  }

   // provide some way to get at letters_
};

class UID {
public:
  typedef ID_store::ID ID;
  static ID_store id; // constructor runs once, single instance
  static ID get_id() {
    return id.get_id();
  }
};
#endif
