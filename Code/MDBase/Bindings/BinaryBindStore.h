#ifndef BINARYBINDSTORE_H
#define BINARYBINDSTORE_H

#include "../../DataStructures/vector1.h"
#include "../../DataStructures/vector1.cpp"

struct BinaryBindStore {
    vector1<bool> isbound;

    vector1<int> boundto;
};

#endif