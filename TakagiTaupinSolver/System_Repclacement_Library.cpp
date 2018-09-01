//
//  System_Repclacement.cpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include "System_Repclacement_Library.hpp"
#include <cmath>
#include <vector>
#include "Vector_And_Operations.hpp"

Rotation_Matrix_Z::Rotation_Matrix_Z(double phi){
    this->phi=phi;
}
void Rotation_Matrix_Z::Basis(double x, double y, double z){
    xc=x*cos(phi)+y*sin(phi);
    yc=-x*sin(phi)+y*cos(phi);
    zc=z;
};
void Rotation_Matrix_Z::Basis_Inverse(double x, double y, double z){
    xcinv=x*cos(phi)-y*sin(phi);
    ycinv=x*sin(phi)+y*cos(phi);
    zcinv=z;
};

Rotation_Matrix_Z_Vector::Rotation_Matrix_Z_Vector(double phi){
    this->phi=phi;
}
void Rotation_Matrix_Z_Vector::Basis(Vector<double> input_vector){
    vector.c[0]=input_vector.c[0]*cos(phi)+input_vector.c[1]*sin(phi);
    vector.c[1]=-input_vector.c[0]*sin(phi)+input_vector.c[1]*cos(phi);
    vector.c[2]=input_vector.c[2];
}
void Rotation_Matrix_Z_Vector::Basis_Inverse(Vector<double> input_vector){
    vector_inv.c[0]=input_vector.c[0]*cos(phi)-input_vector.c[1]*sin(phi);
    vector_inv.c[1]=input_vector.c[0]*sin(phi)+input_vector.c[1]*cos(phi);
    vector_inv.c[2]=input_vector.c[2];
}

Rotation_Matrix_Y::Rotation_Matrix_Y(double phi){
    this->phi=phi;
}
void Rotation_Matrix_Y::Basis(double x, double y, double z){
    xc=x*cos(phi)+z*sin(phi);
    yc=y;
    zc=-x*sin(phi)+z*cos(phi);
}
void Rotation_Matrix_Y::Basis_Inverse(double x, double y, double z){
    xcinv=x*cos(phi)-z*sin(phi);
    ycinv=y;
    zcinv=x*sin(phi)+z*cos(phi);
}

Rotation_Matrix_Y_Vector::Rotation_Matrix_Y_Vector(double phi){
    this->phi=phi;
}
void Rotation_Matrix_Y_Vector::Basis(Vector<double> input_vector){
    vector.c[0]=input_vector.c[0]*cos(phi)+input_vector.c[2]*sin(phi);
    vector.c[1]=input_vector.c[1];
    vector.c[2]=-input_vector.c[0]*sin(phi)+input_vector.c[2]*cos(phi);
}
void Rotation_Matrix_Y_Vector::Basis_Inverse(Vector<double> input_vector){
    vector_inv.c[0]=input_vector.c[0]*cos(phi)-input_vector.c[2]*sin(phi);
    vector_inv.c[1]=input_vector.c[1];
    vector_inv.c[2]=input_vector.c[0]*sin(phi)+input_vector.c[2]*cos(phi);
}

Parallel_Shfit::Parallel_Shfit(double x_shift,double y_shift,double z_shift){
    this->x_shift=x_shift;
    this->y_shift=y_shift;
    this->z_shift=z_shift;
}
void Parallel_Shfit::Basis(double x, double y, double z){
    xc=x-x_shift;
    yc=y-y_shift;
    zc=z-z_shift;
}
void Parallel_Shfit::Basis_Inverse(double x, double y, double z){
    xcinv=x+x_shift;
    ycinv=y+y_shift;
    zcinv=z+z_shift;
}
