#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif // H5_NO_STD
#endif


#include "H5Cpp.h"

int hdf5Test();

#endif // HDF5BUILDER_H_INCLUDE
