### TNetworks                    ###
### Written by: Austin Schneider ###

This library contains generic classes for the representation and manipulation of tensors in a tensor network.
Additionally it contains utilities useful in the implementation of tensor network algorithms.
This library leverages two existing C++ libraries: Eigen, REDSVD.

To use this library Eigen 3 must be installed on the host system.
To compile properly the Eigen library includes should be located in the directory "eigen3" and be visible to the compiler.
The known compatible version is Eigen 3.2.4, and can be found on the Eigen homepage:
http://eigen.tuxfamily.org/
http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz

Components of the REDSVD library come packaged with TNetworks and do not need to be installed separately.
The REDSVD library has been modified to be a template library to reflect the generic nature of Eigen and TNetworks.

To compile:
make
