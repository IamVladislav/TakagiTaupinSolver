//
//  System_Repclacement.cpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include "System_Repclacement.hpp"
#include <cmath>
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
    Rotation_Matrix_Z(double phi){
        this->phi=phi;
    }
    void Basis(double x, double y, double z) override {
        xc=x*cos(phi)+y*sin(phi);
        yc=-x*sin(phi)+y*cos(phi);
        zc=z;
    }
    void Basis_Inverse(double x, double y, double z) override {
        xcinv=x*cos(phi)-y*sin(phi);
        ycinv=x*sin(phi)+y*cos(phi);
        zcinv=z;
    }
};

class Rotation_Matrix_Z_Vector{
public:
    double phi;
    Vector<double> vector,vector_inv;
    Rotation_Matrix_Z_Vector(double phi){
        this->phi=phi;
    }
    void Basis(Vector<double> input_vector) {
        vector.c[0]=input_vector.c[0]*cos(phi)+input_vector.c[1]*sin(phi);
        vector.c[1]=-input_vector.c[0]*sin(phi)+input_vector.c[1]*cos(phi);
        vector.c[2]=input_vector.c[2];
    }
    void Basis_Inverse(Vector<double> input_vector) {
        vector_inv.c[0]=input_vector.c[0]*cos(phi)-input_vector.c[1]*sin(phi);
        vector_inv.c[1]=input_vector.c[0]*sin(phi)+input_vector.c[1]*cos(phi);
        vector_inv.c[2]=input_vector.c[2];
    }
};

class Rotation_Matrix_Y: public System_Replacement{
public:
    double phi;
    Rotation_Matrix_Y(double phi){
        this->phi=phi;
    }
    void Basis(double x, double y, double z) override {
        xc=x*cos(phi)+z*sin(phi);
        yc=y;
        zc=-x*sin(phi)+z*cos(phi);
    }
    void Basis_Inverse(double x, double y, double z) override {
        xcinv=x*cos(phi)-z*sin(phi);
        ycinv=y;
        zcinv=x*sin(phi)+z*cos(phi);
    }
};

class Rotation_Matrix_Y_Vector{
public:
    double phi;
    Vector<double> vector,vector_inv;
    Rotation_Matrix_Y_Vector(double phi){
        this->phi=phi;
    }
    void Basis(Vector<double> input_vector) {
        vector.c[0]=input_vector.c[0]*cos(phi)+input_vector.c[2]*sin(phi);
        vector.c[1]=input_vector.c[1];
        vector.c[2]=-input_vector.c[0]*sin(phi)+input_vector.c[2]*cos(phi);
    }
    void Basis_Inverse(Vector<double> input_vector) {
        vector_inv.c[0]=input_vector.c[0]*cos(phi)-input_vector.c[2]*sin(phi);
        vector_inv.c[1]=input_vector.c[1];
        vector_inv.c[2]=input_vector.c[0]*sin(phi)+input_vector.c[2]*cos(phi);
    }
};

class Parallel_Shfit: public System_Replacement{
public:
    double x_shift, y_shift, z_shift;
    Parallel_Shfit(double x_shift,double y_shift,double z_shift){
        this->x_shift=x_shift;
        this->y_shift=y_shift;
        this->z_shift=z_shift;
    }
    void Basis(double x, double y, double z) override {
        xc=x-x_shift;
        yc=y-y_shift;
        zc=z-z_shift;
    }
    void Basis_Inverse(double x, double y, double z) override {
        xcinv=x+x_shift;
        ycinv=y+y_shift;
        zcinv=z+z_shift;
    }
};
