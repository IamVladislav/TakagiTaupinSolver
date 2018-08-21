//
//  Vector_And_Operations.hpp
//  TakagiTaupinSolver
//
//  Created by Владислав Метель on 18.08.2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#ifndef Vector_And_Operations_hpp
#define Vector_And_Operations_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>

template <typename type> class Vector {
public:
    type c[3];
};

template <typename type> std::ostream &operator<<(std::ostream &os, Vector<type> const &m) {
    return os << "(" << m.c[0] << "," << m.c[1] << "," << m.c[2] << ")";
}

template <typename type> const bool operator == (const Vector<type> &first_vector, const Vector<type> &second_vector)
{
    if ((first_vector.c[0] == second_vector.c[0]) && (first_vector.c[1] == second_vector.c[1]) && (first_vector.c[2] == second_vector.c[2]))
        return true;
    return false;
}

template <typename type> const bool operator != (const Vector<type> &first_vector, const Vector<type> &second_vector)
{
    return !(first_vector == second_vector);
}

template <typename type_of_input, typename type_of_return> Vector<type_of_return> Vector_Multiplication(Vector<type_of_input> first_vector, Vector<type_of_input> second_vector){
    Vector<type_of_return> third_vector;
    third_vector.c[0]=first_vector.c[1]*second_vector.c[2]-first_vector.c[2]*second_vector.c[1];
    third_vector.c[1]=first_vector.c[2]*second_vector.c[0]-first_vector.c[0]*second_vector.c[2];
    third_vector.c[2]=first_vector.c[0]*second_vector.c[1]-first_vector.c[1]*second_vector.c[0];
    return third_vector;
}

template <typename type_of_input, typename type_of_return> Vector<type_of_return> Vector_Normalization(Vector<type_of_input> vector){
    Vector<type_of_return> vector_normalize;
    double abs_vector=sqrt(vector.c[0]*vector.c[0]+vector.c[1]*vector.c[1]+vector.c[2]*vector.c[2]);
    vector_normalize.c[0]=vector.c[0]/abs_vector;
    vector_normalize.c[1]=vector.c[1]/abs_vector;
    vector_normalize.c[2]=vector.c[2]/abs_vector;
    return vector_normalize;
}

template <typename type_of_input> Vector<type_of_input> Vector_Inverse(Vector<type_of_input> vector){
    Vector<type_of_input> inverse_vector;
    inverse_vector.c[0]=-vector.c[0];
    inverse_vector.c[1]=-vector.c[1];
    inverse_vector.c[2]=-vector.c[2];
    return inverse_vector;
}

#endif /* Vector_And_Operations_hpp */
