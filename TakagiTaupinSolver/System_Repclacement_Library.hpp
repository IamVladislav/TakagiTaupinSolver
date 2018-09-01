//
//  System_Repclacement.hpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#ifndef System_Repclacement_Library_hpp
#define System_Repclacement_Library_hpp

#include <stdio.h>
#include <vector>
#include "Vector_And_Operations.hpp"

class System_Replacement{
public:
    double xc, yc, zc;
    double xcinv, ycinv, zcinv;
    virtual void Basis(double x, double y, double z) = 0;
    virtual void Basis_Inverse(double x, double y, double z) = 0;
};

class Rotation_Matrix_Z: public System_Replacement{
public:
    double phi;
    Rotation_Matrix_Z(double phi);
    void Basis(double x, double y, double z);
    void Basis_Inverse(double x, double y, double z);
};

class Rotation_Matrix_Z_Vector{
public:
    double phi;
    Vector<double> vector,vector_inv;
    Rotation_Matrix_Z_Vector(double phi);
    void Basis(Vector<double> input_vector);
    void Basis_Inverse(Vector<double> input_vector);
};

class Rotation_Matrix_Y: public System_Replacement{
public:
    double phi;
    Rotation_Matrix_Y(double phi);
    void Basis(double x, double y, double z);
    void Basis_Inverse(double x, double y, double z);
};

class Rotation_Matrix_Y_Vector{
public:
    double phi;
    Vector<double> vector,vector_inv;
    Rotation_Matrix_Y_Vector(double phi);
    void Basis(Vector<double> input_vector);
    void Basis_Inverse(Vector<double> input_vector);
};

class Parallel_Shfit: public System_Replacement{
public:
    double x_shift, y_shift, z_shift;
    Parallel_Shfit(double x_shift,double y_shift,double z_shift);
    void Basis(double x, double y, double z);
    void Basis_Inverse(double x, double y, double z);
};

#endif /* System_Repclacement_Library_hpp */
