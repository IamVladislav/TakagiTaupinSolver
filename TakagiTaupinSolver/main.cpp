//
//  main.cpp
//  TakagiTaupinSolver
//
//  Created by Владислав Метель on 29.07.2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include "Vector_And_Operations.hpp"

const double PI_2=1.5708;
const double PI=3.14159;

class System_Replacment{
public:
    double xc, yc, zc;
    double xcinv, ycinv, zcinv;
    virtual void Basis(double x, double y, double z) = 0;
    virtual void Basis_Inverse(double x, double y, double z) = 0;
};

class Rotation_Matrix_Z: public System_Replacment{
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

class Rotation_Matrix_Y: public System_Replacment{
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

class Parallel_Shfit: public System_Replacment{
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


class DisplacmentGradient{
public:
    virtual ~DisplacmentGradient(){}
    Vector<double> burgers_vector;
    double depth, nu, kappa;
    double segment_lenght;
    virtual double uxxcalc(double x, double y, double z) = 0;
    virtual double uxycalc(double x, double y, double z) = 0;
    virtual double uxzcalc(double x, double y, double z) = 0;
    virtual double uyxcalc(double x, double y, double z) = 0;
    virtual double uyycalc(double x, double y, double z) = 0;
    virtual double uyzcalc(double x, double y, double z) = 0;
    virtual double uzxcalc(double x, double y, double z) = 0;
    virtual double uzycalc(double x, double y, double z) = 0;
    virtual double uzzcalc(double x, double y, double z) = 0;
};

class AngularDislocation: public DisplacmentGradient{
public:
    AngularDislocation(){
    }
    double uxxcalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(r1 + r2 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r1 + r2 + r5,1.5);
        double r7 = pow(r1 + r2 + r5,2);
        double r8 = pow(r1 + r2 + r5,2.5);
        double r9 = pow(x,3);
        double r10 = sin(kappa);
        double r11 = sqrt(r1 + r2 + r3);
        double r12 = cos(kappa);
        double r13 = sqrt(r1 + r2 + r5);
        double r14 = 1/tan(kappa);
        double r15 = 1/cos(kappa);
        double r16 = pow(r11 - z*r12 - x*r10,2);
        double r17 = pow(2*depth + z + r13,2);
        double r18 = pow(-z + r11,2);
        double r19 = pow(r14,2);
        double r20 = pow(r12,2);
        double r21 = pow(x*r12 - z*r10,2);
        double r22 = pow(r10,2);
        double r23 = pow(r13 + (2*depth + z)*r12 - x*r10,2);
        double r24 = pow(x*r12 + (2*depth + z)*r10,2);
        double r25 = pow(r2*r12 + x*(x*r12 - z*r10),2);
        double r26 = pow(r2*r12 + x*(x*r12 + (2*depth + z)*r10),2);
        double uxx=(burgers_vector.c[2]*(y*r10*((-1 + (x*r10)/r11)/(r11*(r11 - z*r12 - x*r10)) - ((x/r11 - r10)*(-x + r11*r10))/(r11*r16) - (x*(-x + r11*r10))/(r4*(r11 - z*r12 - x*r10)) + (-1 + (x*r10)/r13)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - ((x/r13 - r10)*(-x + r13*r10))/(r13*r23) - (x*(-x + r13*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10))) + 2*(-((y*(depth + z)*((-2*depth*x)/r7 - x/(r13*r17)))/r13) + (x*y*(depth + z)*(depth/(r1 + r2 + r5) + 1/(2*depth + z + r13)))/r6 + (1 - 2*nu)*(-((x*y*(1 + depth/r13))/(r13*r17)) - (depth*x*y)/(r6*(2*depth + z + r13)) + (y*r12*(depth/r13 + r12)*(x/r13 - r10))/r23 + (depth*x*y*r12)/(r6*(r13 + (2*depth + z)*r12 - x*r10))) + (y*(depth + z)*r12*((-2*depth*x*(2*depth + z))/r7 - ((depth/r13 + r12)*(2*depth + z + r13*r12)*(x/r13 - r10))/r23 + (x*r12*(depth/r13 + r12))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - (depth*x*(2*depth + z + r13*r12))/(r6*(r13 + (2*depth + z)*r12 - x*r10))))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - (y*(depth + z)*r12*(x/r13 - r10)*((depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r13 + r12)*(2*depth + z + r13*r12))/(r13 + (2*depth + z)*r12 - x*r10)))/(r13*r23) - (x*y*(depth + z)*r12*((depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r13 + r12)*(2*depth + z + r13*r12))/(r13 + (2*depth + z)*r12 - x*r10)))/(r6*(r13 + (2*depth + z)*r12 - x*r10)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(r1*(-(x/((r1 + r2 + r3)*r18)) - x/(r4*(-z + r11)) - x/((r1 + r2 + r5)*r17) - x/(r6*(2*depth + z + r13))) + 2*x*(1/(r11*(-z + r11)) + 1/(r13*(2*depth + z + r13))) + (6*depth*r1*(depth + z)*r14)/r8 - (2*depth*(depth + z)*r14)/r6 + (2*(1 - 2*nu)*(-((r9*(nu + depth/r13))/(r13*r17)) - (depth*r9)/(r6*(2*depth + z + r13)) + (2*x*(nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu)*r14 - (depth*r1*r14)/r6 + (depth*r14)/r13))/(2*depth + z + r13) - (2*(1 - 2*nu)*x*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu)*x*r14 + (depth*x*r14)/r13))/(r13*r17) + (2*(depth + z)*((3*depth*r9)/r8 - (2*depth*x)/r6 + (r9*(2*nu + depth/r11))/((r1 + r2 + r5)*r17) + (depth*r9)/(r4*r13*(2*depth + z + r13)) + (r9*(2*nu + depth/r11))/(r6*(2*depth + z + r13)) - (2*x*(2*nu + depth/r11))/(r13*(2*depth + z + r13)) + ((1 - 2*nu)*r14)/r13 - (x*(depth + (1 - 2*nu)*x*r14))/r6))/(2*depth + z + r13) - (2*x*(depth + z)*(2*nu - (depth*r1)/r6 - (r1*(2*nu + depth/r11))/(r13*(2*depth + z + r13)) + (depth + (1 - 2*nu)*x*r14)/r13))/(r13*r17) + ((x*r12 - z*r10)*(-1 + (x*r10)/r11))/(r11*(r11 - z*r12 - x*r10)) + (r12*(-x + r11*r10))/(r11*(r11 - z*r12 - x*r10)) - ((x/r11 - r10)*(x*r12 - z*r10)*(-x + r11*r10))/(r11*r16) - (x*(x*r12 - z*r10)*(-x + r11*r10))/(r4*(r11 - z*r12 - x*r10)) + ((x*r12 + (2*depth + z)*r10)*(-1 + (x*r10)/r13))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) + (r12*(-x + r13*r10))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - ((x/r13 - r10)*(x*r12 + (2*depth + z)*r10)*(-x + r13*r10))/(r13*r23) - (x*(x*r12 + (2*depth + z)*r10)*(-x + r13*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10)) + 2*(1 - 2*nu)*((x*(nu + 2*(1 - nu)*r19))/(r13*(2*depth + z + r13)) - (r12*(1 + 2*(1 - nu)*r19)*(x/r13 - r10))/(r13 + (2*depth + z)*r12 - x*r10)) - (2*(1 - 2*nu)*r14*(r20 - (depth*r15*(-1 + (x*r10)/r13))/r13 + (depth*x*r15*(-x + r13*r10))/r6))/(r13 + (2*depth + z)*r12 - x*r10) + (2*(1 - 2*nu)*r14*(x/r13 - r10)*(r12*(x*r12 + (2*depth + z)*r10) - (depth*r15*(-x + r13*r10))/r13))/r23 + (2*(depth + z)*r14*((-3*depth*r1*(2*depth + z)*r15)/r8 + (depth*(2*depth + z)*r15)/r6 + ((-x + r13*r10)*(((2*depth + z + r13*r12)*(1 + (depth*r15)/r13)*(x/r13 - r10))/r23 + (depth*x*(2*depth + z + r13*r12)*r15)/(r6*(r13 + (2*depth + z)*r12 - x*r10)) - (x*r12*(1 + (depth*r15)/r13))/(r13*(r13 + (2*depth + z)*r12 - x*r10))))/r13 + ((-1 + (x*r10)/r13)*(2*(1 - nu)*r12 - ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/r13 - (x*(-x + r13*r10)*(2*(1 - nu)*r12 - ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/r6))/(r13 + (2*depth + z)*r12 - x*r10) - (2*(depth + z)*r14*(x/r13 - r10)*((depth*x*(2*depth + z)*r15)/r6 - r12*r10 + ((-x + r13*r10)*(2*(1 - nu)*r12 - ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/r13))/r23 + (-1 + 2*nu)*(x/(r11*(-z + r11)) + x/(r13*(2*depth + z + r13)) - r12*((x/r11 - r10)/(r11 - z*r12 - x*r10) + (x/r13 - r10)/(r13 + (2*depth + z)*r12 - x*r10)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(x*y*(-(x/((r1 + r2 + r3)*r18)) - x/(r4*(-z + r11)) - x/((r1 + r2 + r5)*r17) - x/(r6*(2*depth + z + r13)))) - y*(1/(r11*(-z + r11)) + 1/(r13*(2*depth + z + r13))) - y*r12*((-1 + (x*r10)/r11)/(r11*(r11 - z*r12 - x*r10)) - ((x/r11 - r10)*(-x + r11*r10))/(r11*r16) - (x*(-x + r11*r10))/(r4*(r11 - z*r12 - x*r10)) + (-1 + (x*r10)/r13)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - ((x/r13 - r10)*(-x + r13*r10))/(r13*r23) - (x*(-x + r13*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10))) + 2*(1 - nu)*((2*y)/(r1 + r2) - (y*r12)/(r21*(1 + r2/r21)) - (y*r12)/(r24*(1 + r2/r24)) + (-((y*r11*r10*(2*x*r12 - z*r10))/r25) + (x*y*r10)/(r11*(r2*r12 + x*(x*r12 - z*r10))))/(1 + (r2*(r1 + r2 + r3)*r22)/r25) + (-((y*r13*r10*(2*x*r12 + (2*depth + z)*r10))/r26) + (x*y*r10)/(r13*(r2*r12 + x*(x*r12 + (2*depth + z)*r10))))/(1 + (r2*(r1 + r2 + r5)*r22)/r26)) + 2*((y*(depth + z)*((-2*depth*r1)/r7 + depth/(r1 + r2 + r5) - (r1*(2*nu + depth/r13))/(r13*r17) - (depth*r1)/(r6*(2*depth + z + r13)) + (2*nu + depth/r13)/(2*depth + z + r13)))/(r13*(2*depth + z + r13)) - (3*depth*x*y*(depth + z)*r14)/r8 - (x*y*(depth + z)*((depth*x)/(r1 + r2 + r5) + (x*(2*nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu)*r14))/((r1 + r2 + r5)*r17) - (x*y*(depth + z)*((depth*x)/(r1 + r2 + r5) + (x*(2*nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu)*r14))/(r6*(2*depth + z + r13)) + ((1 - 2*nu)*y*((r1*(nu + depth/r13))/(r13*r17) + (depth*r1)/(r6*(2*depth + z + r13)) - (nu + depth/r13)/(2*depth + z + r13) + (depth*x*r14)/r6))/(2*depth + z + r13) - ((1 - 2*nu)*x*y*(-((x*(nu + depth/r13))/(2*depth + z + r13)) + (1 - 2*nu - depth/r13)*r14))/(r13*r17) - ((1 - 2*nu)*y*r12*(depth/r13 + r12)*r14*(x/r13 - r10))/r23 - (depth*(1 - 2*nu)*x*y*r12*r14)/(r6*(r13 + (2*depth + z)*r12 - x*r10)) + (y*(depth + z)*((2*depth*x*(2*depth + z)*r12*r14)/r7 + (r12*((x*r12*(-(depth/r13) + (1 - 2*nu)*r12)*r14)/r13 + (depth*x*(2*depth + z + r13*r12)*r14)/r6 + 2*(1 - nu)*r12*(-1 + (x*r10)/r13)))/(r13 + (2*depth + z)*r12 - x*r10) - (r12*(x/r13 - r10)*((-(depth/r13) + (1 - 2*nu)*r12)*(2*depth + z + r13*r12)*r14 + 2*(1 - nu)*r12*(-x + r13*r10)))/r23))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - (y*(depth + z)*(x/r13 - r10)*(-((depth*(2*depth + z)*r12*r14)/(r1 + r2 + r5)) + (r12*((-(depth/r13) + (1 - 2*nu)*r12)*(2*depth + z + r13*r12)*r14 + 2*(1 - nu)*r12*(-x + r13*r10)))/(r13 + (2*depth + z)*r12 - x*r10)))/(r13*r23) - (x*y*(depth + z)*(-((depth*(2*depth + z)*r12*r14)/(r1 + r2 + r5)) + (r12*((-(depth/r13) + (1 - 2*nu)*r12)*(2*depth + z + r13*r12)*r14 + 2*(1 - nu)*r12*(-x + r13*r10)))/(r13 + (2*depth + z)*r12 - x*r10)))/(r6*(r13 + (2*depth + z)*r12 - x*r10)) - 2*(1 - 2*nu)*(1 - nu)*r19*(y/(r1 + r2) - (y*r12)/(r24*(1 + r2/r24)) + (-((y*r13*r10*(2*x*r12 + (2*depth + z)*r10))/r26) + (x*y*r10)/(r13*(r2*r12 + x*(x*r12 + (2*depth + z)*r10))))/(1 + (r2*(r1 + r2 + r5)*r22)/r26)))))/(8.*PI*(1 - nu));
        return uxx;
    }
    double uxycalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(2*depth + z,2);
        double r5 = pow(r1 + r2 + r3,1.5);
        double r6 = pow(r1 + r2 + r4,1.5);
        double r7 = pow(r1 + r2 + r4,2);
        double r8 = pow(r1 + r2 + r4,2.5);
        double r9 = sin(kappa);
        double r10 = sqrt(r1 + r2 + r3);
        double r11 = cos(kappa);
        double r12 = sqrt(r1 + r2 + r4);
        double r13 = 1/tan(kappa);
        double r14 = 1/cos(kappa);
        double r15 = tan(kappa);
        double r16 = pow(r10 - z*r11 - x*r9,2);
        double r17 = pow(2*depth + z + r12,2);
        double r18 = pow(-z + r10,2);
        double r19 = pow(x*r11 - z*r9,2);
        double r20 = pow(r9,2);
        double r21 = pow(r13,2);
        double r22 = pow(r12 + (2*depth + z)*r11 - x*r9,2);
        double r23 = pow(x*r11 + (2*depth + z)*r9,2);
        double r24 = pow(r2*r11 + x*(x*r11 - z*r9),2);
        double r25 = pow(r2*r11 + x*(x*r11 + (2*depth + z)*r9),2);
        double uxy=(burgers_vector.c[2]*(y*r9*((y*r9)/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + (y*r9)/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(-x + r10*r9))/((r1 + r2 + r3)*r16) - (y*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) - (y*(-x + r12*r9))/((r1 + r2 + r4)*r22) - (y*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + r9*((-x + r10*r9)/(r10*(r10 - z*r11 - x*r9)) + (-x + r12*r9)/(r12*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(-((y*(depth + z)*((-2*depth*y)/r7 - y/(r12*r17)))/r12) + (r2*(depth + z)*(depth/(r1 + r2 + r4) + 1/(2*depth + z + r12)))/r6 - ((depth + z)*(depth/(r1 + r2 + r4) + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*(-((r2*(1 + depth/r12))/(r12*r17)) - (depth*r2)/(r6*(2*depth + z + r12)) + (1 + depth/r12)/(2*depth + z + r12) + (r2*r11*(depth/r12 + r11))/(r12*r22) + (depth*r2*r11)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) - (r11*(depth/r12 + r11))/(r12 + (2*depth + z)*r11 - x*r9)) + (y*(depth + z)*r11*((-2*depth*y*(2*depth + z))/r7 - (y*(depth/r12 + r11)*(2*depth + z + r12*r11))/(r12*r22) + (y*r11*(depth/r12 + r11))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (depth*y*(2*depth + z + r12*r11))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (r2*(depth + z)*r11*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/((r1 + r2 + r4)*r22) - (r2*(depth + z)*r11*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + ((depth + z)*r11*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(x*y*(-(y/((r1 + r2 + r3)*r18)) - y/(r5*(-z + r10)) - y/((r1 + r2 + r4)*r17) - y/(r6*(2*depth + z + r12)))) - x*(1/(r10*(-z + r10)) + 1/(r12*(2*depth + z + r12))) - y*r11*((y*r9)/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + (y*r9)/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(-x + r10*r9))/((r1 + r2 + r3)*r16) - (y*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) - (y*(-x + r12*r9))/((r1 + r2 + r4)*r22) - (y*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) - r11*((-x + r10*r9)/(r10*(r10 - z*r11 - x*r9)) + (-x + r12*r9)/(r12*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(1 - nu)*((-2*x)/(r1 + r2) + 1/((x*r11 - z*r9)*(1 + r2/r19)) + 1/((x*r11 + (2*depth + z)*r9)*(1 + r2/r23)) + ((-2*r2*r10*r11*r9)/r24 + (r2*r9)/(r10*(r2*r11 + x*(x*r11 - z*r9))) + (r10*r9)/(r2*r11 + x*(x*r11 - z*r9)))/(1 + (r2*(r1 + r2 + r3)*r20)/r24) + ((-2*r2*r12*r11*r9)/r25 + (r2*r9)/(r12*(r2*r11 + x*(x*r11 + (2*depth + z)*r9))) + (r12*r9)/(r2*r11 + x*(x*r11 + (2*depth + z)*r9)))/(1 + (r2*(r1 + r2 + r4)*r20)/r25)) + 2*((y*(depth + z)*((-2*depth*x*y)/r7 - (x*y*(2*nu + depth/r12))/(r12*r17) - (depth*x*y)/(r6*(2*depth + z + r12))))/(r12*(2*depth + z + r12)) - (3*depth*r2*(depth + z)*r13)/r8 + (depth*(depth + z)*r13)/r6 - (r2*(depth + z)*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/((r1 + r2 + r4)*r17) - (r2*(depth + z)*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/(r6*(2*depth + z + r12)) + ((depth + z)*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/(r12*(2*depth + z + r12)) + ((1 - 2*nu)*y*((x*y*(nu + depth/r12))/(r12*r17) + (depth*x*y)/(r6*(2*depth + z + r12)) + (depth*y*r13)/r6))/(2*depth + z + r12) - ((1 - 2*nu)*r2*(-((x*(nu + depth/r12))/(2*depth + z + r12)) + (1 - 2*nu - depth/r12)*r13))/(r12*r17) + ((1 - 2*nu)*(-((x*(nu + depth/r12))/(2*depth + z + r12)) + (1 - 2*nu - depth/r12)*r13))/(2*depth + z + r12) - ((1 - 2*nu)*r2*r11*(depth/r12 + r11)*r13)/(r12*r22) - (depth*(1 - 2*nu)*r2*r11*r13)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + ((1 - 2*nu)*r11*(depth/r12 + r11)*r13)/(r12 + (2*depth + z)*r11 - x*r9) + (y*(depth + z)*((2*depth*y*(2*depth + z)*r11*r13)/r7 + (r11*((y*r11*(-(depth/r12) + (1 - 2*nu)*r11)*r13)/r12 + (depth*y*(2*depth + z + r12*r11)*r13)/r6 + (2*(1 - nu)*y*r11*r9)/r12))/(r12 + (2*depth + z)*r11 - x*r9) - (y*r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12*r22)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (r2*(depth + z)*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/((r1 + r2 + r4)*r22) - (r2*(depth + z)*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + ((depth + z)*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - 2*(1 - 2*nu)*(1 - nu)*r21*(-(x/(r1 + r2)) + 1/((x*r11 + (2*depth + z)*r9)*(1 + r2/r23)) + ((-2*r2*r12*r11*r9)/r25 + (r2*r9)/(r12*(r2*r11 + x*(x*r11 + (2*depth + z)*r9))) + (r12*r9)/(r2*r11 + x*(x*r11 + (2*depth + z)*r9)))/(1 + (r2*(r1 + r2 + r4)*r20)/r25)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(r1*(-(y/((r1 + r2 + r3)*r18)) - y/(r5*(-z + r10)) - y/((r1 + r2 + r4)*r17) - y/(r6*(2*depth + z + r12))) + (6*depth*x*y*(depth + z)*r13)/r8 + (2*(1 - 2*nu)*(-((r1*y*(nu + depth/r12))/(r12*r17)) - (depth*r1*y)/(r6*(2*depth + z + r12)) - (depth*x*y*r13)/r6))/(2*depth + z + r12) - (2*(1 - 2*nu)*y*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*x*r13 + (depth*x*r13)/r12))/(r12*r17) + (2*(depth + z)*((3*depth*r1*y)/r8 + (r1*y*(2*nu + depth/r10))/((r1 + r2 + r4)*r17) + (depth*r1*y)/(r5*r12*(2*depth + z + r12)) + (r1*y*(2*nu + depth/r10))/(r6*(2*depth + z + r12)) - (y*(depth + (1 - 2*nu)*x*r13))/r6))/(2*depth + z + r12) - (2*y*(depth + z)*(2*nu - (depth*r1)/r6 - (r1*(2*nu + depth/r10))/(r12*(2*depth + z + r12)) + (depth + (1 - 2*nu)*x*r13)/r12))/(r12*r17) + (y*r9*(x*r11 - z*r9))/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + (y*r9*(x*r11 + (2*depth + z)*r9))/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(x*r11 - z*r9)*(-x + r10*r9))/((r1 + r2 + r3)*r16) - (y*(x*r11 - z*r9)*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) - (y*(x*r11 + (2*depth + z)*r9)*(-x + r12*r9))/((r1 + r2 + r4)*r22) - (y*(x*r11 + (2*depth + z)*r9)*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + 2*(1 - 2*nu)*((y*(nu + 2*(1 - nu)*r21))/(r12*(2*depth + z + r12)) - (y*r11*(1 + 2*(1 - nu)*r21))/(r12*(r12 + (2*depth + z)*r11 - x*r9))) + (2*(1 - 2*nu)*y*r13*(r11*(x*r11 + (2*depth + z)*r9) - (depth*r14*(-x + r12*r9))/r12))/(r12*r22) + (-1 + 2*nu)*(y/(r10*(-z + r10)) + y/(r12*(2*depth + z + r12)) - r11*(y/(r10*(r10 - z*r11 - x*r9)) + y/(r12*(r12 + (2*depth + z)*r11 - x*r9)))) + (2*(depth + z)*r13*((-3*depth*x*y*(2*depth + z)*r14)/r8 + ((-x + r12*r9)*((y*(2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12*r22) + (depth*y*(2*depth + z + r12*r11)*r14)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) - (y*r11*(1 + (depth*r14)/r12))/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/r12 + (y*r9*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r1 + r2 + r4) - (y*(-x + r12*r9)*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r6))/(r12 + (2*depth + z)*r11 - x*r9) - (2*y*(depth + z)*r13*((depth*x*(2*depth + z)*r14)/r6 - r11*r9 + ((-x + r12*r9)*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r12))/(r12*r22) - (2*(1 - 2*nu)*r13*((depth*y*r14*(-x + r12*r9))/r6 - (depth*y*r15)/(r1 + r2 + r4)))/(r12 + (2*depth + z)*r11 - x*r9)))/(8.*PI*(1 - nu));
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(2*depth + z,2);
        double r5 = pow(r1 + r2 + r3,1.5);
        double r6 = pow(r1 + r2 + r4,1.5);
        double r7 = pow(r1 + r2 + r4,2);
        double r8 = pow(r1 + r2 + r4,2.5);
        double r9 = sin(kappa);
        double r10 = sqrt(r1 + r2 + r3);
        double r11 = cos(kappa);
        double r12 = sqrt(r1 + r2 + r4);
        double r13 = 1/tan(kappa);
        double r14 = 1/cos(kappa);
        double r15 = tan(kappa);
        double r16 = pow(r10 - z*r11 - x*r9,2);
        double r17 = pow(2*depth + z + r12,2);
        double r18 = pow(-z + r10,2);
        double r19 = pow(x*r11 - z*r9,2);
        double r20 = pow(r9,2);
        double r21 = pow(r13,2);
        double r22 = pow(r12 + (2*depth + z)*r11 - x*r9,2);
        double r23 = pow(x*r11 + (2*depth + z)*r9,2);
        double r24 = pow(r2*r11 + x*(x*r11 - z*r9),2);
        double r25 = pow(r2*r11 + x*(x*r11 + (2*depth + z)*r9),2);
        double uxz=(burgers_vector.c[2]*(y*r9*((z*r9)/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + ((2*depth + z)*r9)/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - ((z/r10 - r11)*(-x + r10*r9))/(r10*r16) - (z*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) - (((2*depth + z)/r12 + r11)*(-x + r12*r9))/(r12*r22) - ((2*depth + z)*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(-((y*(depth + z)*((-2*depth*(2*depth + z))/r7 - (1 + (2*depth + z)/r12)/r17))/r12) + (y*(depth + z)*(2*depth + z)*(depth/(r1 + r2 + r4) + 1/(2*depth + z + r12)))/r6 - (y*(depth/(r1 + r2 + r4) + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*(-((y*(1 + depth/r12)*(1 + (2*depth + z)/r12))/r17) - (depth*y*(2*depth + z))/(r6*(2*depth + z + r12)) + (y*r11*(depth/r12 + r11)*((2*depth + z)/r12 + r11))/r22 + (depth*y*(2*depth + z)*r11)/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + (y*(depth + z)*r11*((-2*depth*r4)/r7 + depth/(r1 + r2 + r4) - ((depth/r12 + r11)*((2*depth + z)/r12 + r11)*(2*depth + z + r12*r11))/r22 + ((depth/r12 + r11)*(1 + ((2*depth + z)*r11)/r12))/(r12 + (2*depth + z)*r11 - x*r9) - (depth*(2*depth + z)*(2*depth + z + r12*r11))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(depth + z)*r11*((2*depth + z)/r12 + r11)*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*r22) - (y*(depth + z)*(2*depth + z)*r11*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*r11*((depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r12 + r11)*(2*depth + z + r12*r11))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(x*y*(-((-1 + z/r10)/(r10*r18)) - z/(r5*(-z + r10)) - (1 + (2*depth + z)/r12)/(r12*r17) - (2*depth + z)/(r6*(2*depth + z + r12)))) - y*r11*((z*r9)/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + ((2*depth + z)*r9)/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - ((z/r10 - r11)*(-x + r10*r9))/(r10*r16) - (z*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) - (((2*depth + z)/r12 + r11)*(-x + r12*r9))/(r12*r22) - ((2*depth + z)*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(1 - nu)*((y*r9)/(r19*(1 + r2/r19)) - (y*r9)/(r23*(1 + r2/r23)) + ((x*y*r10*r20)/r24 + (y*z*r9)/(r10*(r2*r11 + x*(x*r11 - z*r9))))/(1 + (r2*(r1 + r2 + r3)*r20)/r24) + (-((x*y*r12*r20)/r25) + (y*(2*depth + z)*r9)/(r12*(r2*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r2*(r1 + r2 + r4)*r20)/r25)) + 2*((y*(depth + z)*((-2*depth*x*(2*depth + z))/r7 - (x*(2*nu + depth/r12)*(1 + (2*depth + z)/r12))/r17 - (depth*x*(2*depth + z))/(r6*(2*depth + z + r12))))/(r12*(2*depth + z + r12)) - (3*depth*y*(depth + z)*(2*depth + z)*r13)/r8 + (depth*y*r13)/r6 - (y*(depth + z)*(1 + (2*depth + z)/r12)*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/(r12*r17) - (y*(depth + z)*(2*depth + z)*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/(r6*(2*depth + z + r12)) + (y*((depth*x)/(r1 + r2 + r4) + (x*(2*nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*r13))/(r12*(2*depth + z + r12)) + ((1 - 2*nu)*y*((x*(nu + depth/r12)*(1 + (2*depth + z)/r12))/r17 + (depth*x*(2*depth + z))/(r6*(2*depth + z + r12)) + (depth*(2*depth + z)*r13)/r6))/(2*depth + z + r12) - ((1 - 2*nu)*y*(1 + (2*depth + z)/r12)*(-((x*(nu + depth/r12))/(2*depth + z + r12)) + (1 - 2*nu - depth/r12)*r13))/r17 - ((1 - 2*nu)*y*r11*(depth/r12 + r11)*((2*depth + z)/r12 + r11)*r13)/r22 - (depth*(1 - 2*nu)*y*(2*depth + z)*r11*r13)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*(depth + z)*((2*depth*r4*r11*r13)/r7 - (depth*r11*r13)/(r1 + r2 + r4) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(1 + ((2*depth + z)*r11)/r12)*r13 + (depth*(2*depth + z)*(2*depth + z + r12*r11)*r13)/r6 + (2*(1 - nu)*(2*depth + z)*r11*r9)/r12))/(r12 + (2*depth + z)*r11 - x*r9) - (r11*((2*depth + z)/r12 + r11)*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/r22))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(depth + z)*((2*depth + z)/r12 + r11)*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*r22) - (y*(depth + z)*(2*depth + z)*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*(-((depth*(2*depth + z)*r11*r13)/(r1 + r2 + r4)) + (r11*((-(depth/r12) + (1 - 2*nu)*r11)*(2*depth + z + r12*r11)*r13 + 2*(1 - nu)*r11*(-x + r12*r9)))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - 2*(1 - 2*nu)*(1 - nu)*r21*(-((y*r9)/(r23*(1 + r2/r23))) + (-((x*y*r12*r20)/r25) + (y*(2*depth + z)*r9)/(r12*(r2*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r2*(r1 + r2 + r4)*r20)/r25)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(r1*(-((-1 + z/r10)/(r10*r18)) - z/(r5*(-z + r10)) - (1 + (2*depth + z)/r12)/(r12*r17) - (2*depth + z)/(r6*(2*depth + z + r12))) + (6*depth*x*(depth + z)*(2*depth + z)*r13)/r8 - (2*depth*x*r13)/r6 + (2*(1 - 2*nu)*(nu - (r1*(nu + depth/r12)*(1 + (2*depth + z)/r12))/r17 - (depth*r1*(2*depth + z))/(r6*(2*depth + z + r12)) - (depth*x*(2*depth + z)*r13)/r6))/(2*depth + z + r12) - (2*(1 - 2*nu)*(1 + (2*depth + z)/r12)*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu)*x*r13 + (depth*x*r13)/r12))/r17 + (2*(depth + z)*((3*depth*r1*(2*depth + z))/r8 + (r1*(2*nu + depth/r10)*(1 + (2*depth + z)/r12))/(r12*r17) + (depth*r1*z)/(r5*r12*(2*depth + z + r12)) + (r1*(2*depth + z)*(2*nu + depth/r10))/(r6*(2*depth + z + r12)) - ((2*depth + z)*(depth + (1 - 2*nu)*x*r13))/r6))/(2*depth + z + r12) - (2*(depth + z)*(1 + (2*depth + z)/r12)*(2*nu - (depth*r1)/r6 - (r1*(2*nu + depth/r10))/(r12*(2*depth + z + r12)) + (depth + (1 - 2*nu)*x*r13)/r12))/r17 + (2*(2*nu - (depth*r1)/r6 - (r1*(2*nu + depth/r10))/(r12*(2*depth + z + r12)) + (depth + (1 - 2*nu)*x*r13)/r12))/(2*depth + z + r12) + (z*r9*(x*r11 - z*r9))/((r1 + r2 + r3)*(r10 - z*r11 - x*r9)) + ((2*depth + z)*r9*(x*r11 + (2*depth + z)*r9))/((r1 + r2 + r4)*(r12 + (2*depth + z)*r11 - x*r9)) - (r9*(-x + r10*r9))/(r10*(r10 - z*r11 - x*r9)) - ((z/r10 - r11)*(x*r11 - z*r9)*(-x + r10*r9))/(r10*r16) - (z*(x*r11 - z*r9)*(-x + r10*r9))/(r5*(r10 - z*r11 - x*r9)) + (r9*(-x + r12*r9))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (((2*depth + z)/r12 + r11)*(x*r11 + (2*depth + z)*r9)*(-x + r12*r9))/(r12*r22) - ((2*depth + z)*(x*r11 + (2*depth + z)*r9)*(-x + r12*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + 2*(1 - 2*nu)*(((1 + (2*depth + z)/r12)*(nu + 2*(1 - nu)*r21))/(2*depth + z + r12) - (r11*((2*depth + z)/r12 + r11)*(1 + 2*(1 - nu)*r21))/(r12 + (2*depth + z)*r11 - x*r9)) + (2*(1 - 2*nu)*((2*depth + z)/r12 + r11)*r13*(r11*(x*r11 + (2*depth + z)*r9) - (depth*r14*(-x + r12*r9))/r12))/r22 + (-1 + 2*nu)*((-1 + z/r10)/(-z + r10) + (1 + (2*depth + z)/r12)/(2*depth + z + r12) - r11*((z/r10 - r11)/(r10 - z*r11 - x*r9) + ((2*depth + z)/r12 + r11)/(r12 + (2*depth + z)*r11 - x*r9))) + (2*(depth + z)*r13*((-3*depth*x*r4*r14)/r8 + (depth*x*r14)/r6 + ((-x + r12*r9)*((((2*depth + z)/r12 + r11)*(2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/r22 + (depth*(2*depth + z)*(2*depth + z + r12*r11)*r14)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) - ((1 + ((2*depth + z)*r11)/r12)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r12 + ((2*depth + z)*r9*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r1 + r2 + r4) - ((2*depth + z)*(-x + r12*r9)*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r6))/(r12 + (2*depth + z)*r11 - x*r9) - (2*(depth + z)*((2*depth + z)/r12 + r11)*r13*((depth*x*(2*depth + z)*r14)/r6 - r11*r9 + ((-x + r12*r9)*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r12))/r22 + (2*r13*((depth*x*(2*depth + z)*r14)/r6 - r11*r9 + ((-x + r12*r9)*(2*(1 - nu)*r11 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/r12))/(r12 + (2*depth + z)*r11 - x*r9) - (2*(1 - 2*nu)*r13*(r11*r9 + (depth*(2*depth + z)*r14*(-x + r12*r9))/r6 - (depth*(2*depth + z)*r15)/(r1 + r2 + r4)))/(r12 + (2*depth + z)*r11 - x*r9)))/(8.*PI*(1 - nu));
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override
    {
        double r1 = pow(y,2);
        double r2 = pow(x,2);
        double r3 = pow(z,2);
        double r4 = pow(r2 + r1 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r2 + r1 + r5,1.5);
        double r7 = pow(r2 + r1 + r5,2);
        double r8 = pow(r2 + r1 + r5,2.5);
        double r9 = sin(kappa);
        double r10 = sqrt(r2 + r1 + r3);
        double r11 = cos(kappa);
        double r12 = sqrt(r2 + r1 + r5);
        double r13 = 1/tan(kappa);
        double r14 = 1/cos(kappa);
        double r15 = 1/sin(kappa);
        double r16 = pow(r10 - z*r11 - x*r9,2);
        double r17 = pow(2*depth + z + r12,2);
        double r18 = pow(-z + r10,2);
        double r19 = pow(r13,2);
        double r20 = pow(r11,2);
        double r21 = pow(x*r11 - z*r9,2);
        double r22 = pow(r9,2);
        double r23 = pow(r12 + (2*depth + z)*r11 - x*r9,2);
        double r24 = pow(x*r11 + (2*depth + z)*r9,2);
        double r25 = pow(r1*r11 + x*(x*r11 - z*r9),2);
        double r26 = pow(r1*r11 + x*(x*r11 + (2*depth + z)*r9),2);
        double uyx=(burgers_vector.c[2]*(-(r1*r9*(-((x/r10 - r9)/(r10*r16)) - x/(r4*(r10 - z*r11 - x*r9)) - (x/r12 - r9)/(r12*r23) - x/(r6*(r12 + (2*depth + z)*r11 - x*r9)))) + (1 - 2*nu)*r9*((x/r10 - r9)/(r10 - z*r11 - x*r9) + (x/r12 - r9)/(r12 + (2*depth + z)*r11 - x*r9)) + 2*((x*(depth + z)*((-2*depth*x)/r7 - x/(r12*r17)))/r12 - (r2*(depth + z)*(depth/(r2 + r1 + r5) + 1/(2*depth + z + r12)))/r6 + ((depth + z)*(depth/(r2 + r1 + r5) + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*((r2*(1 + depth/r12))/(r12*r17) + (depth*r2)/(r6*(2*depth + z + r12)) - (1 + depth/r12)/(2*depth + z + r12) + (r11*(depth/r12 + r11))/(r12 + (2*depth + z)*r11 - x*r9) - ((x/r12 - r9)*r9)/(r12 + (2*depth + z)*r11 - x*r9) - ((depth/r12 + r11)*(x/r12 - r9)*(x*r11 + (2*depth + z)*r9))/r23 - (depth*x*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) - ((depth + z)*(((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*r11)/r12 + (depth*x*r9)/r6 - (2*depth*x*(2*depth + z)*(x*r11 + (2*depth + z)*r9))/r8 - (x*(1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r6 - (-((depth*r11*(2*depth + z + r12*r11))/r12) - (depth*x*r11*(x*r11 + (2*depth + z)*r9))/(r2 + r1 + r5) + (depth*x*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r6)/(r12*(r12 + (2*depth + z)*r11 - x*r9)) + ((x/r12 - r9)*(r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12))/(r12*r23) + (x*(r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12 + (2*depth + z)*r11 - x*r9) + ((depth + z)*(x/r12 - r9)*((-(depth/r12) + r11)*r9 + ((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/r23)))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(r1*(-(x/((r2 + r1 + r3)*r18)) - x/(r4*(-z + r10)) - x/((r2 + r1 + r5)*r17) - x/(r6*(2*depth + z + r12)) - r11*(-((x/r10 - r9)/(r10*r16)) - x/(r4*(r10 - z*r11 - x*r9)) - (x/r12 - r9)/(r12*r23) - x/(r6*(r12 + (2*depth + z)*r11 - x*r9))))) + (1 - 2*nu)*(x/(r10*(-z + r10)) + x/(r12*(2*depth + z + r12)) - r11*((x/r10 - r9)/(r10 - z*r11 - x*r9) + (x/r12 - r9)/(r12 + (2*depth + z)*r11 - x*r9))) + 2*((3*depth*r2*(depth + z)*r13)/r8 - (depth*(depth + z)*r13)/r6 - ((1 - 2*nu)*(-((x*r1*(nu + depth/r12))/(r12*r17)) - (depth*x*r1)/(r6*(2*depth + z + r12)) + (depth*r2*r13)/r6 + (1 - 2*nu - depth/r12)*r13))/(2*depth + z + r12) + ((1 - 2*nu)*x*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r12))/(2*depth + z + r12) + x*(1 - 2*nu - depth/r12)*r13))/(r12*r17) + ((depth + z)*((-3*depth*x*r1)/r8 - (x*r1*(2*nu + depth/r12))/((r2 + r1 + r5)*r17) - (depth*x*r1)/(r7*(2*depth + z + r12)) - (x*r1*(2*nu + depth/r12))/(r6*(2*depth + z + r12)) + ((1 - 2*nu)*r13)/r12 - (x*(-depth + (1 - 2*nu)*x*r13))/r6))/(2*depth + z + r12) - (x*(depth + z)*(-2*nu + (depth*r1)/r6 + (r1*(2*nu + depth/r12))/(r12*(2*depth + z + r12)) + (-depth + (1 - 2*nu)*x*r13)/r12))/(r12*r17) - ((1 - 2*nu)*r11*(depth/r12 + r11)*r13)/(r12 + (2*depth + z)*r11 - x*r9) + ((1 - 2*nu)*(depth/r12 + r11)*r13*(x/r12 - r9)*(x*r11 + (2*depth + z)*r9))/r23 + (depth*(1 - 2*nu)*x*r13*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (1 - 2*nu)*((x*(-nu + 2*(1 - nu)*r19))/(r12*(2*depth + z + r12)) - (r11*(1 - 2*nu + 2*(1 - nu)*r19)*(x/r12 - r9))/(r12 + (2*depth + z)*r11 - x*r9)) + ((depth + z)*((depth*(2*depth + z)*r11*r13)/r6 - ((1 - 2*nu)*r11*r13)/r12 - (3*depth*x*(2*depth + z)*r13*(x*r11 + (2*depth + z)*r9))/r8 + (x*(depth*r11 + (1 - 2*nu)*r13*(x*r11 + (2*depth + z)*r9)))/r6 - (-((depth*r11*(2*depth + z + r12*r11)*r13)/r12) - (depth*x*r11*r13*(x*r11 + (2*depth + z)*r9))/(r2 + r1 + r5) + (depth*x*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r6)/(r12*(r12 + (2*depth + z)*r11 - x*r9)) + ((x/r12 - r9)*(r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12))/(r12*r23) + (x*(r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12 + (2*depth + z)*r11 - x*r9) - ((depth + z)*(x/r12 - r9)*(r20 + (depth*(2*depth + z)*r13*(x*r11 + (2*depth + z)*r9))/r6 - (depth*r11 + (1 - 2*nu)*r13*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/r23)))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(x*y*(-(x/((r2 + r1 + r3)*r18)) - x/(r4*(-z + r10)) - x/((r2 + r1 + r5)*r17) - x/(r6*(2*depth + z + r12))) + y*(1/(r10*(-z + r10)) + 1/(r12*(2*depth + z + r12))) - y*(r11/(r10*(r10 - z*r11 - x*r9)) + r11/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - ((x/r10 - r9)*(x*r11 - z*r9))/(r10*r16) - (x*(x*r11 - z*r9))/(r4*(r10 - z*r11 - x*r9)) - ((x/r12 - r9)*(x*r11 + (2*depth + z)*r9))/(r12*r23) - (x*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(1 - nu)*((2*y)/(r2 + r1) - (y*r11)/(r21*(1 + r1/r21)) - (y*r11)/(r24*(1 + r1/r24)) + (-((y*r10*r9*(2*x*r11 - z*r9))/r25) + (x*y*r9)/(r10*(r1*r11 + x*(x*r11 - z*r9))))/(1 + (r1*(r2 + r1 + r3)*r22)/r25) + (-((y*r12*r9*(2*x*r11 + (2*depth + z)*r9))/r26) + (x*y*r9)/(r12*(r1*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r1*(r2 + r1 + r5)*r22)/r26)) + 2*((y*(depth + z)*((2*nu*r2)/(r12*r17) - (2*nu)/(2*depth + z + r12) - (depth*x*(-(x/r6) - x/(r12*r17)))/r12 + (depth*r2*(1/r12 + 1/(2*depth + z + r12)))/r6 - (depth*(1/r12 + 1/(2*depth + z + r12)))/r12))/(r12*(2*depth + z + r12)) + (3*depth*x*y*(depth + z)*r13)/r8 - (x*y*(depth + z)*((-2*nu*x)/(2*depth + z + r12) - (depth*x*(1/r12 + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*r13))/((r2 + r1 + r5)*r17) - (x*y*(depth + z)*((-2*nu*x)/(2*depth + z + r12) - (depth*x*(1/r12 + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*r13))/(r6*(2*depth + z + r12)) + ((1 - 2*nu)*y*(-((r2*(nu + depth/r12))/(r12*r17)) - (depth*r2)/(r6*(2*depth + z + r12)) + (nu + depth/r12)/(2*depth + z + r12) - (depth*x*r13)/r6))/(2*depth + z + r12) - ((1 - 2*nu)*x*y*((x*(nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu + depth/r12)*r13))/(r12*r17) + ((1 - 2*nu)*y*r13*(1 + (depth*r14)/r12)*(x/r12 - r9))/r23 + (depth*(1 - 2*nu)*x*y*r15)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*(depth + z)*r13*((-2*depth*x*(2*depth + z)*r14)/r7 - ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12)*(x/r12 - r9))/r23 - (depth*x*(2*depth + z + r12*r11)*r14)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (x*r11*(1 + (depth*r14)/r12))/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(depth + z)*r13*(x/r12 - r9)*(-2*(1 - nu)*r11 + (depth*(2*depth + z)*r14)/(r2 + r1 + r5) + ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*r23) - (x*y*(depth + z)*r13*(-2*(1 - nu)*r11 + (depth*(2*depth + z)*r14)/(r2 + r1 + r5) + ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + 2*(1 - 2*nu)*(1 - nu)*r19*(y/(r2 + r1) - (y*r11)/(r24*(1 + r1/r24)) + (-((y*r12*r9*(2*x*r11 + (2*depth + z)*r9))/r26) + (x*y*r9)/(r12*(r1*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r1*(r2 + r1 + r5)*r22)/r26)))))/(8.*PI*(1 - nu));
        return uyx;
    }
    double uyycalc(double x, double y, double z) override
    {
        double r1 = pow(y,2);
        double r2 = pow(x,2);
        double r3 = pow(z,2);
        double r4 = pow(r2 + r1 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r2 + r1 + r5,1.5);
        double r7 = pow(r2 + r1 + r5,2);
        double r8 = pow(r2 + r1 + r5,2.5);
        double r9 = pow(y,3);
        double r10 = sin(kappa);
        double r11 = sqrt(r2 + r1 + r3);
        double r12 = cos(kappa);
        double r13 = sqrt(r2 + r1 + r5);
        double r14 = 1/tan(kappa);
        double r15 = 1/cos(kappa);
        double r16 = 1/sin(kappa);
        double r17 = pow(r11 - z*r12 - x*r10,2);
        double r18 = pow(2*depth + z + r13,2);
        double r19 = pow(-z + r11,2);
        double r20 = pow(r14,2);
        double r21 = pow(r12,2);
        double r22 = pow(x*r12 - z*r10,2);
        double r23 = pow(r10,2);
        double r24 = pow(r13 + (2*depth + z)*r12 - x*r10,2);
        double r25 = pow(x*r12 + (2*depth + z)*r10,2);
        double r26 = pow(r1*r12 + x*(x*r12 - z*r10),2);
        double r27 = pow(r1*r12 + x*(x*r12 + (2*depth + z)*r10),2);
        double uyy=(burgers_vector.c[2]*(-(r1*r10*(-(y/((r2 + r1 + r3)*r17)) - y/(r4*(r11 - z*r12 - x*r10)) - y/((r2 + r1 + r5)*r24) - y/(r6*(r13 + (2*depth + z)*r12 - x*r10)))) - 2*y*r10*(1/(r11*(r11 - z*r12 - x*r10)) + 1/(r13*(r13 + (2*depth + z)*r12 - x*r10))) + (1 - 2*nu)*r10*(y/(r11*(r11 - z*r12 - x*r10)) + y/(r13*(r13 + (2*depth + z)*r12 - x*r10))) + 2*((x*(depth + z)*((-2*depth*y)/r7 - y/(r13*r18)))/r13 - (x*y*(depth + z)*(depth/(r2 + r1 + r5) + 1/(2*depth + z + r13)))/r6 + (1 - 2*nu)*((x*y*(1 + depth/r13))/(r13*r18) + (depth*x*y)/(r6*(2*depth + z + r13)) - (y*r10)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - (y*(depth/r13 + r12)*(x*r12 + (2*depth + z)*r10))/(r13*r24) - (depth*y*(x*r12 + (2*depth + z)*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10))) - ((depth + z)*((depth*y*r10)/r6 - (2*depth*y*(2*depth + z)*(x*r12 + (2*depth + z)*r10))/r8 - (y*(1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r12 + (2*depth + z)*r10))/r6 - (2*y*r12*r10 - (depth*y*r12*(x*r12 + (2*depth + z)*r10))/(r2 + r1 + r5) + (depth*y*(2*depth + z + r13*r12)*(x*r12 + (2*depth + z)*r10))/r6)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) + (y*(r1*r12*r10 - (depth*(2*depth + z + r13*r12)*(x*r12 + (2*depth + z)*r10))/r13))/((r2 + r1 + r5)*r24) + (y*(r1*r12*r10 - (depth*(2*depth + z + r13*r12)*(x*r12 + (2*depth + z)*r10))/r13))/(r6*(r13 + (2*depth + z)*r12 - x*r10))))/(r13 + (2*depth + z)*r12 - x*r10) + (y*(depth + z)*((-(depth/r13) + r12)*r10 + ((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r12 + (2*depth + z)*r10))/r13 - (r1*r12*r10 - (depth*(2*depth + z + r13*r12)*(x*r12 + (2*depth + z)*r10))/r13)/(r13*(r13 + (2*depth + z)*r12 - x*r10))))/(r13*r24))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(r1*(-(y/((r2 + r1 + r3)*r19)) - y/(r4*(-z + r11)) - y/((r2 + r1 + r5)*r18) - y/(r6*(2*depth + z + r13)) - r12*(-(y/((r2 + r1 + r3)*r17)) - y/(r4*(r11 - z*r12 - x*r10)) - y/((r2 + r1 + r5)*r24) - y/(r6*(r13 + (2*depth + z)*r12 - x*r10))))) - 2*y*(1/(r11*(-z + r11)) + 1/(r13*(2*depth + z + r13)) - r12*(1/(r11*(r11 - z*r12 - x*r10)) + 1/(r13*(r13 + (2*depth + z)*r12 - x*r10)))) + (1 - 2*nu)*(y/(r11*(-z + r11)) + y/(r13*(2*depth + z + r13)) - r12*(y/(r11*(r11 - z*r12 - x*r10)) + y/(r13*(r13 + (2*depth + z)*r12 - x*r10)))) + 2*((3*depth*x*y*(depth + z)*r14)/r8 - ((1 - 2*nu)*(-((r9*(nu + depth/r13))/(r13*r18)) - (depth*r9)/(r6*(2*depth + z + r13)) + (2*y*(nu + depth/r13))/(2*depth + z + r13) + (depth*x*y*r14)/r6))/(2*depth + z + r13) + ((1 - 2*nu)*y*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r13))/(2*depth + z + r13) + x*(1 - 2*nu - depth/r13)*r14))/(r13*r18) + ((depth + z)*((-3*depth*r9)/r8 + (2*depth*y)/r6 - (r9*(2*nu + depth/r13))/((r2 + r1 + r5)*r18) - (depth*r9)/(r7*(2*depth + z + r13)) - (r9*(2*nu + depth/r13))/(r6*(2*depth + z + r13)) + (2*y*(2*nu + depth/r13))/(r13*(2*depth + z + r13)) - (y*(-depth + (1 - 2*nu)*x*r14))/r6))/(2*depth + z + r13) - (y*(depth + z)*(-2*nu + (depth*r1)/r6 + (r1*(2*nu + depth/r13))/(r13*(2*depth + z + r13)) + (-depth + (1 - 2*nu)*x*r14)/r13))/(r13*r18) + ((1 - 2*nu)*y*(depth/r13 + r12)*r14*(x*r12 + (2*depth + z)*r10))/(r13*r24) + (depth*(1 - 2*nu)*y*r14*(x*r12 + (2*depth + z)*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10)) + (1 - 2*nu)*((y*(-nu + 2*(1 - nu)*r20))/(r13*(2*depth + z + r13)) - (y*r12*(1 - 2*nu + 2*(1 - nu)*r20))/(r13*(r13 + (2*depth + z)*r12 - x*r10))) + ((depth + z)*((-3*depth*y*(2*depth + z)*r14*(x*r12 + (2*depth + z)*r10))/r8 + (y*(depth*r12 + (1 - 2*nu)*r14*(x*r12 + (2*depth + z)*r10)))/r6 - (2*y*r21 - (depth*y*r12*r14*(x*r12 + (2*depth + z)*r10))/(r2 + r1 + r5) + (depth*y*(2*depth + z + r13*r12)*r14*(x*r12 + (2*depth + z)*r10))/r6)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) + (y*(r1*r21 - (depth*(2*depth + z + r13*r12)*r14*(x*r12 + (2*depth + z)*r10))/r13))/((r2 + r1 + r5)*r24) + (y*(r1*r21 - (depth*(2*depth + z + r13*r12)*r14*(x*r12 + (2*depth + z)*r10))/r13))/(r6*(r13 + (2*depth + z)*r12 - x*r10))))/(r13 + (2*depth + z)*r12 - x*r10) - (y*(depth + z)*(r21 + (depth*(2*depth + z)*r14*(x*r12 + (2*depth + z)*r10))/r6 - (depth*r12 + (1 - 2*nu)*r14*(x*r12 + (2*depth + z)*r10))/r13 - (r1*r21 - (depth*(2*depth + z + r13*r12)*r14*(x*r12 + (2*depth + z)*r10))/r13)/(r13*(r13 + (2*depth + z)*r12 - x*r10))))/(r13*r24))))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(x*y*(-(y/((r2 + r1 + r3)*r19)) - y/(r4*(-z + r11)) - y/((r2 + r1 + r5)*r18) - y/(r6*(2*depth + z + r13))) + x*(1/(r11*(-z + r11)) + 1/(r13*(2*depth + z + r13))) - (x*r12 - z*r10)/(r11*(r11 - z*r12 - x*r10)) - (x*r12 + (2*depth + z)*r10)/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - y*(-((y*(x*r12 - z*r10))/((r2 + r1 + r3)*r17)) - (y*(x*r12 - z*r10))/(r4*(r11 - z*r12 - x*r10)) - (y*(x*r12 + (2*depth + z)*r10))/((r2 + r1 + r5)*r24) - (y*(x*r12 + (2*depth + z)*r10))/(r6*(r13 + (2*depth + z)*r12 - x*r10))) + 2*(1 - nu)*((-2*x)/(r2 + r1) + 1/((x*r12 - z*r10)*(1 + r1/r22)) + 1/((x*r12 + (2*depth + z)*r10)*(1 + r1/r25)) + ((-2*r1*r11*r12*r10)/r26 + (r1*r10)/(r11*(r1*r12 + x*(x*r12 - z*r10))) + (r11*r10)/(r1*r12 + x*(x*r12 - z*r10)))/(1 + (r1*(r2 + r1 + r3)*r23)/r26) + ((-2*r1*r13*r12*r10)/r27 + (r1*r10)/(r13*(r1*r12 + x*(x*r12 + (2*depth + z)*r10))) + (r13*r10)/(r1*r12 + x*(x*r12 + (2*depth + z)*r10)))/(1 + (r1*(r2 + r1 + r5)*r23)/r27)) + 2*((y*(depth + z)*((2*nu*x*y)/(r13*r18) - (depth*x*(-(y/r6) - y/(r13*r18)))/r13 + (depth*x*y*(1/r13 + 1/(2*depth + z + r13)))/r6))/(r13*(2*depth + z + r13)) + (3*depth*r1*(depth + z)*r14)/r8 - (depth*(depth + z)*r14)/r6 - (r1*(depth + z)*((-2*nu*x)/(2*depth + z + r13) - (depth*x*(1/r13 + 1/(2*depth + z + r13)))/r13 + (1 - 2*nu)*r14))/((r2 + r1 + r5)*r18) - (r1*(depth + z)*((-2*nu*x)/(2*depth + z + r13) - (depth*x*(1/r13 + 1/(2*depth + z + r13)))/r13 + (1 - 2*nu)*r14))/(r6*(2*depth + z + r13)) + ((depth + z)*((-2*nu*x)/(2*depth + z + r13) - (depth*x*(1/r13 + 1/(2*depth + z + r13)))/r13 + (1 - 2*nu)*r14))/(r13*(2*depth + z + r13)) + ((1 - 2*nu)*y*(-((x*y*(nu + depth/r13))/(r13*r18)) - (depth*x*y)/(r6*(2*depth + z + r13)) - (depth*y*r14)/r6))/(2*depth + z + r13) - ((1 - 2*nu)*r1*((x*(nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu + depth/r13)*r14))/(r13*r18) + ((1 - 2*nu)*((x*(nu + depth/r13))/(2*depth + z + r13) + (-1 + 2*nu + depth/r13)*r14))/(2*depth + z + r13) + ((1 - 2*nu)*r1*r14*(1 + (depth*r15)/r13))/(r13*r24) + (depth*(1 - 2*nu)*r1*r16)/(r6*(r13 + (2*depth + z)*r12 - x*r10)) - ((1 - 2*nu)*r14*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10) + (y*(depth + z)*r14*((-2*depth*y*(2*depth + z)*r15)/r7 - (y*(2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13*r24) - (depth*y*(2*depth + z + r13*r12)*r15)/(r6*(r13 + (2*depth + z)*r12 - x*r10)) + (y*r12*(1 + (depth*r15)/r13))/(r13*(r13 + (2*depth + z)*r12 - x*r10))))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) - (r1*(depth + z)*r14*(-2*(1 - nu)*r12 + (depth*(2*depth + z)*r15)/(r2 + r1 + r5) + ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/((r2 + r1 + r5)*r24) - (r1*(depth + z)*r14*(-2*(1 - nu)*r12 + (depth*(2*depth + z)*r15)/(r2 + r1 + r5) + ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/(r6*(r13 + (2*depth + z)*r12 - x*r10)) + ((depth + z)*r14*(-2*(1 - nu)*r12 + (depth*(2*depth + z)*r15)/(r2 + r1 + r5) + ((2*depth + z + r13*r12)*(1 + (depth*r15)/r13))/(r13 + (2*depth + z)*r12 - x*r10)))/(r13*(r13 + (2*depth + z)*r12 - x*r10)) + 2*(1 - 2*nu)*(1 - nu)*r20*(-(x/(r2 + r1)) + 1/((x*r12 + (2*depth + z)*r10)*(1 + r1/r25)) + ((-2*r1*r13*r12*r10)/r27 + (r1*r10)/(r13*(r1*r12 + x*(x*r12 + (2*depth + z)*r10))) + (r13*r10)/(r1*r12 + x*(x*r12 + (2*depth + z)*r10)))/(1 + (r1*(r2 + r1 + r5)*r23)/r27)))))/(8.*PI*(1 - nu));
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override
    {
        double r1 = pow(y,2);
        double r2 = pow(x,2);
        double r3 = pow(z,2);
        double r4 = pow(r2 + r1 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r2 + r1 + r5,1.5);
        double r7 = pow(r2 + r1 + r5,2);
        double r8 = pow(r2 + r1 + r5,2.5);
        double r9 = sin(kappa);
        double r10 = sqrt(r2 + r1 + r3);
        double r11 = cos(kappa);
        double r12 = sqrt(r2 + r1 + r5);
        double r13 = 1/tan(kappa);
        double r14 = 1/cos(kappa);
        double r15 = 1/sin(kappa);
        double r16 = pow(r10 - z*r11 - x*r9,2);
        double r17 = pow(2*depth + z + r12,2);
        double r18 = pow(-z + r10,2);
        double r19 = pow(r13,2);
        double r20 = pow(r11,2);
        double r21 = pow(x*r11 - z*r9,2);
        double r22 = pow(r9,2);
        double r23 = pow(r12 + (2*depth + z)*r11 - x*r9,2);
        double r24 = pow(x*r11 + (2*depth + z)*r9,2);
        double r25 = pow(r1*r11 + x*(x*r11 - z*r9),2);
        double r26 = pow(r1*r11 + x*(x*r11 + (2*depth + z)*r9),2);
        double uyz=(burgers_vector.c[2]*(-(r1*r9*(-((z/r10 - r11)/(r10*r16)) - z/(r4*(r10 - z*r11 - x*r9)) - ((2*depth + z)/r12 + r11)/(r12*r23) - (2*depth + z)/(r6*(r12 + (2*depth + z)*r11 - x*r9)))) + (1 - 2*nu)*r9*((z/r10 - r11)/(r10 - z*r11 - x*r9) + ((2*depth + z)/r12 + r11)/(r12 + (2*depth + z)*r11 - x*r9)) + 2*((x*(depth + z)*((-2*depth*(2*depth + z))/r7 - (1 + (2*depth + z)/r12)/r17))/r12 - (x*(depth + z)*(2*depth + z)*(depth/(r2 + r1 + r5) + 1/(2*depth + z + r12)))/r6 + (x*(depth/(r2 + r1 + r5) + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*((x*(1 + depth/r12)*(1 + (2*depth + z)/r12))/r17 + (depth*x*(2*depth + z))/(r6*(2*depth + z + r12)) + ((depth/r12 + r11)*r9)/(r12 + (2*depth + z)*r11 - x*r9) - (((2*depth + z)/r12 + r11)*r9)/(r12 + (2*depth + z)*r11 - x*r9) - ((depth/r12 + r11)*((2*depth + z)/r12 + r11)*(x*r11 + (2*depth + z)*r9))/r23 - (depth*(2*depth + z)*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) - ((depth + z)*((depth*(2*depth + z)*r9)/r6 + ((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*r9)/r12 + (((-2*depth*r5)/r7 + depth/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r12 - ((2*depth + z)*(1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r6 - (-((depth*(2*depth + z + r12*r11)*r9)/r12) - (depth*(1 + ((2*depth + z)*r11)/r12)*(x*r11 + (2*depth + z)*r9))/r12 + (depth*(2*depth + z)*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r6)/(r12*(r12 + (2*depth + z)*r11 - x*r9)) + (((2*depth + z)/r12 + r11)*(r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12))/(r12*r23) + ((2*depth + z)*(r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12 + (2*depth + z)*r11 - x*r9) + ((depth + z)*((2*depth + z)/r12 + r11)*((-(depth/r12) + r11)*r9 + ((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/r23 - ((-(depth/r12) + r11)*r9 + ((1 + (depth*(2*depth + z))/(r2 + r1 + r5))*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r11*r9 - (depth*(2*depth + z + r12*r11)*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9)))/(r12 + (2*depth + z)*r11 - x*r9))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(-(r1*(-((-1 + z/r10)/(r10*r18)) - z/(r4*(-z + r10)) - (1 + (2*depth + z)/r12)/(r12*r17) - (2*depth + z)/(r6*(2*depth + z + r12)) - r11*(-((z/r10 - r11)/(r10*r16)) - z/(r4*(r10 - z*r11 - x*r9)) - ((2*depth + z)/r12 + r11)/(r12*r23) - (2*depth + z)/(r6*(r12 + (2*depth + z)*r11 - x*r9))))) + (1 - 2*nu)*((-1 + z/r10)/(-z + r10) + (1 + (2*depth + z)/r12)/(2*depth + z + r12) - r11*((z/r10 - r11)/(r10 - z*r11 - x*r9) + ((2*depth + z)/r12 + r11)/(r12 + (2*depth + z)*r11 - x*r9))) + 2*((3*depth*x*(depth + z)*(2*depth + z)*r13)/r8 - (depth*x*r13)/r6 - ((1 - 2*nu)*(nu - (r1*(nu + depth/r12)*(1 + (2*depth + z)/r12))/r17 - (depth*r1*(2*depth + z))/(r6*(2*depth + z + r12)) + (depth*x*(2*depth + z)*r13)/r6))/(2*depth + z + r12) + ((1 - 2*nu)*(1 + (2*depth + z)/r12)*(-depth + nu*(2*depth + z) + (r1*(nu + depth/r12))/(2*depth + z + r12) + x*(1 - 2*nu - depth/r12)*r13))/r17 + ((depth + z)*((-3*depth*r1*(2*depth + z))/r8 - (r1*(2*nu + depth/r12)*(1 + (2*depth + z)/r12))/(r12*r17) - (depth*r1*(2*depth + z))/(r7*(2*depth + z + r12)) - (r1*(2*depth + z)*(2*nu + depth/r12))/(r6*(2*depth + z + r12)) - ((2*depth + z)*(-depth + (1 - 2*nu)*x*r13))/r6))/(2*depth + z + r12) - ((depth + z)*(1 + (2*depth + z)/r12)*(-2*nu + (depth*r1)/r6 + (r1*(2*nu + depth/r12))/(r12*(2*depth + z + r12)) + (-depth + (1 - 2*nu)*x*r13)/r12))/r17 + (-2*nu + (depth*r1)/r6 + (r1*(2*nu + depth/r12))/(r12*(2*depth + z + r12)) + (-depth + (1 - 2*nu)*x*r13)/r12)/(2*depth + z + r12) - ((1 - 2*nu)*r11*(depth/r12 + r11))/(r12 + (2*depth + z)*r11 - x*r9) + ((1 - 2*nu)*(depth/r12 + r11)*((2*depth + z)/r12 + r11)*r13*(x*r11 + (2*depth + z)*r9))/r23 + (depth*(1 - 2*nu)*(2*depth + z)*r13*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (1 - 2*nu)*(((1 + (2*depth + z)/r12)*(-nu + 2*(1 - nu)*r19))/(2*depth + z + r12) - (r11*((2*depth + z)/r12 + r11)*(1 - 2*nu + 2*(1 - nu)*r19))/(r12 + (2*depth + z)*r11 - x*r9)) + ((depth + z)*((depth*(2*depth + z)*r11)/r6 - ((1 - 2*nu)*r11)/r12 - (3*depth*r5*r13*(x*r11 + (2*depth + z)*r9))/r8 + (depth*r13*(x*r11 + (2*depth + z)*r9))/r6 + ((2*depth + z)*(depth*r11 + (1 - 2*nu)*r13*(x*r11 + (2*depth + z)*r9)))/r6 - (-((depth*r11*(2*depth + z + r12*r11))/r12) - (depth*(1 + ((2*depth + z)*r11)/r12)*r13*(x*r11 + (2*depth + z)*r9))/r12 + (depth*(2*depth + z)*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r6)/(r12*(r12 + (2*depth + z)*r11 - x*r9)) + (((2*depth + z)/r12 + r11)*(r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12))/(r12*r23) + ((2*depth + z)*(r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12))/(r6*(r12 + (2*depth + z)*r11 - x*r9))))/(r12 + (2*depth + z)*r11 - x*r9) - ((depth + z)*((2*depth + z)/r12 + r11)*(r20 + (depth*(2*depth + z)*r13*(x*r11 + (2*depth + z)*r9))/r6 - (depth*r11 + (1 - 2*nu)*r13*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9))))/r23 + (r20 + (depth*(2*depth + z)*r13*(x*r11 + (2*depth + z)*r9))/r6 - (depth*r11 + (1 - 2*nu)*r13*(x*r11 + (2*depth + z)*r9))/r12 - (r1*r20 - (depth*(2*depth + z + r12*r11)*r13*(x*r11 + (2*depth + z)*r9))/r12)/(r12*(r12 + (2*depth + z)*r11 - x*r9)))/(r12 + (2*depth + z)*r11 - x*r9))))/(8.*PI*(1 - nu)) + (burgers_vector.c[1]*(x*y*(-((-1 + z/r10)/(r10*r18)) - z/(r4*(-z + r10)) - (1 + (2*depth + z)/r12)/(r12*r17) - (2*depth + z)/(r6*(2*depth + z + r12))) - y*(-(r9/(r10*(r10 - z*r11 - x*r9))) + r9/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - ((z/r10 - r11)*(x*r11 - z*r9))/(r10*r16) - (z*(x*r11 - z*r9))/(r4*(r10 - z*r11 - x*r9)) - (((2*depth + z)/r12 + r11)*(x*r11 + (2*depth + z)*r9))/(r12*r23) - ((2*depth + z)*(x*r11 + (2*depth + z)*r9))/(r6*(r12 + (2*depth + z)*r11 - x*r9))) + 2*(1 - nu)*((y*r9)/(r21*(1 + r1/r21)) - (y*r9)/(r24*(1 + r1/r24)) + ((x*y*r10*r22)/r25 + (y*z*r9)/(r10*(r1*r11 + x*(x*r11 - z*r9))))/(1 + (r1*(r2 + r1 + r3)*r22)/r25) + (-((x*y*r12*r22)/r26) + (y*(2*depth + z)*r9)/(r12*(r1*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r1*(r2 + r1 + r5)*r22)/r26)) + 2*((y*(depth + z)*((2*nu*x*(1 + (2*depth + z)/r12))/r17 - (depth*x*(-((2*depth + z)/r6) - (1 + (2*depth + z)/r12)/r17))/r12 + (depth*x*(2*depth + z)*(1/r12 + 1/(2*depth + z + r12)))/r6))/(r12*(2*depth + z + r12)) + (3*depth*y*(depth + z)*(2*depth + z)*r13)/r8 - (depth*y*r13)/r6 - (y*(depth + z)*(1 + (2*depth + z)/r12)*((-2*nu*x)/(2*depth + z + r12) - (depth*x*(1/r12 + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*r13))/(r12*r17) - (y*(depth + z)*(2*depth + z)*((-2*nu*x)/(2*depth + z + r12) - (depth*x*(1/r12 + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*r13))/(r6*(2*depth + z + r12)) + (y*((-2*nu*x)/(2*depth + z + r12) - (depth*x*(1/r12 + 1/(2*depth + z + r12)))/r12 + (1 - 2*nu)*r13))/(r12*(2*depth + z + r12)) + ((1 - 2*nu)*y*(-((x*(nu + depth/r12)*(1 + (2*depth + z)/r12))/r17) - (depth*x*(2*depth + z))/(r6*(2*depth + z + r12)) - (depth*(2*depth + z)*r13)/r6))/(2*depth + z + r12) - ((1 - 2*nu)*y*(1 + (2*depth + z)/r12)*((x*(nu + depth/r12))/(2*depth + z + r12) + (-1 + 2*nu + depth/r12)*r13))/r17 + ((1 - 2*nu)*y*((2*depth + z)/r12 + r11)*r13*(1 + (depth*r14)/r12))/r23 + (depth*(1 - 2*nu)*y*(2*depth + z)*r15)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*(depth + z)*r13*((-2*depth*r5*r14)/r7 + (depth*r14)/(r2 + r1 + r5) - (((2*depth + z)/r12 + r11)*(2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/r23 - (depth*(2*depth + z)*(2*depth + z + r12*r11)*r14)/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + ((1 + ((2*depth + z)*r11)/r12)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) - (y*(depth + z)*((2*depth + z)/r12 + r11)*r13*(-2*(1 - nu)*r11 + (depth*(2*depth + z)*r14)/(r2 + r1 + r5) + ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*r23) - (y*(depth + z)*(2*depth + z)*r13*(-2*(1 - nu)*r11 + (depth*(2*depth + z)*r14)/(r2 + r1 + r5) + ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r6*(r12 + (2*depth + z)*r11 - x*r9)) + (y*r13*(-2*(1 - nu)*r11 + (depth*(2*depth + z)*r14)/(r2 + r1 + r5) + ((2*depth + z + r12*r11)*(1 + (depth*r14)/r12))/(r12 + (2*depth + z)*r11 - x*r9)))/(r12*(r12 + (2*depth + z)*r11 - x*r9)) + 2*(1 - 2*nu)*(1 - nu)*r19*(-((y*r9)/(r24*(1 + r1/r24))) + (-((x*y*r12*r22)/r26) + (y*(2*depth + z)*r9)/(r12*(r1*r11 + x*(x*r11 + (2*depth + z)*r9))))/(1 + (r1*(r2 + r1 + r5)*r22)/r26)))))/(8.*PI*(1 - nu));
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(2*depth + z,2);
        double r5 = pow(r1 + r2 + r3,1.5);
        double r6 = pow(r1 + r2 + r4,1.5);
        double r7 = pow(r1 + r2 + r4,2);
        double r8 = sqrt(r1 + r2 + r3);
        double r9 = sqrt(r1 + r2 + r4);
        double r10 = cos(kappa);
        double r11 = sin(kappa);
        double r12 = 1/tan(kappa);
        double r13 = pow(r8 - z*r10 - x*r11,2);
        double r14 = pow(2*depth + z + r9,2);
        double r15 = pow(x*r10 - z*r11,2);
        double r16 = pow(r11,2);
        double r17 = pow(r9 + (2*depth + z)*r10 - x*r11,2);
        double r18 = pow(x*r10 + (2*depth + z)*r11,2);
        double r19 = pow(r2*r10 + x*(x*r10 - z*r11),2);
        double r20 = pow(r2*r10 + x*(x*r10 + (2*depth + z)*r11),2);
        double uzx=(burgers_vector.c[1]*(-(1/r8) + 1/r9 - x*(-(x/r5) + x/r6) + (r10*(-z + r8*r10))/(r8*(r8 - z*r10 - x*r11)) - (r10*(2*depth + z + r9*r10))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) - ((-z + r8*r10)*(x/r8 - r11)*(x*r10 - z*r11))/(r8*r13) + (x*r10*(x*r10 - z*r11))/((r1 + r2 + r3)*(r8 - z*r10 - x*r11)) - (x*(-z + r8*r10)*(x*r10 - z*r11))/(r5*(r8 - z*r10 - x*r11)) + ((2*depth + z + r9*r10)*(x/r9 - r11)*(x*r10 + (2*depth + z)*r11))/(r9*r17) - (x*r10*(x*r10 + (2*depth + z)*r11))/((r1 + r2 + r4)*(r9 + (2*depth + z)*r10 - x*r11)) + (x*(2*depth + z + r9*r10)*(x*r10 + (2*depth + z)*r11))/(r6*(r9 + (2*depth + z)*r10 - x*r11)) + (-1 + 2*nu)*r11*((x/r8 - r11)/(r8 - z*r10 - x*r11) - (x/r9 - r11)/(r9 + (2*depth + z)*r10 - x*r11)) + 2*((2*(1 - nu)*r1*(2*nu + depth/r9))/(r9*r14) + (2*depth*(1 - nu)*r1)/(r6*(2*depth + z + r9)) - (2*(1 - nu)*(2*nu + depth/r9))/(2*depth + z + r9) + ((depth + z)*((2*depth*r1)/r7 - depth/(r1 + r2 + r4) + (2*nu*r1)/(r9*r14) - (2*nu)/(2*depth + z + r9)))/r9 - (x*(depth + z)*(-((depth*x)/(r1 + r2 + r4)) - (2*nu*x)/(2*depth + z + r9) + (1 - 2*nu)*r12))/r6 + (2*(1 - nu)*r10*(depth/r9 + r10))/(r9 + (2*depth + z)*r10 - x*r11) - (2*(1 - nu)*(depth/r9 + r10)*(x/r9 - r11)*(x*r10 + (2*depth + z)*r11))/r17 - (2*depth*(1 - nu)*x*(x*r10 + (2*depth + z)*r11))/(r6*(r9 + (2*depth + z)*r10 - x*r11)) - 2*(1 - 2*nu)*(1 - nu)*r12*(x/(r9*(2*depth + z + r9)) - (r10*(x/r9 - r11))/(r9 + (2*depth + z)*r10 - x*r11)) - ((depth + z)*(((2*depth + z + r9*r10)*r12*(((2*depth + z + r9*r10)*(x/r9 - r11))/r17 - (x*r10)/(r9*(r9 + (2*depth + z)*r10 - x*r11))))/r9 + (x*r10*r12*(2*(1 - nu)*r10 - (2*depth + z + r9*r10)/(r9 + (2*depth + z)*r10 - x*r11)))/(r1 + r2 + r4) - (x*(2*depth + z + r9*r10)*r12*(2*(1 - nu)*r10 - (2*depth + z + r9*r10)/(r9 + (2*depth + z)*r10 - x*r11)))/r6 + (depth*(-(((2*depth + z)*r10)/(r1 + r2 + r4)) - (r10*(2*depth + z + r9*r10))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) + (2*x*(2*depth + z)*(x*r10 + (2*depth + z)*r11))/r7 + ((2*depth + z + r9*r10)*(x/r9 - r11)*(x*r10 + (2*depth + z)*r11))/(r9*r17) - (x*r10*(x*r10 + (2*depth + z)*r11))/((r1 + r2 + r4)*(r9 + (2*depth + z)*r10 - x*r11)) + (x*(2*depth + z + r9*r10)*(x*r10 + (2*depth + z)*r11))/(r6*(r9 + (2*depth + z)*r10 - x*r11))))/r9 - (depth*x*(r11 - ((2*depth + z)*(x*r10 + (2*depth + z)*r11))/(r1 + r2 + r4) - ((2*depth + z + r9*r10)*(x*r10 + (2*depth + z)*r11))/(r9*(r9 + (2*depth + z)*r10 - x*r11))))/r6))/(r9 + (2*depth + z)*r10 - x*r11) + ((depth + z)*(x/r9 - r11)*(r10*r11 + ((2*depth + z + r9*r10)*r12*(2*(1 - nu)*r10 - (2*depth + z + r9*r10)/(r9 + (2*depth + z)*r10 - x*r11)))/r9 + (depth*(r11 - ((2*depth + z)*(x*r10 + (2*depth + z)*r11))/(r1 + r2 + r4) - ((2*depth + z + r9*r10)*(x*r10 + (2*depth + z)*r11))/(r9*(r9 + (2*depth + z)*r10 - x*r11))))/r9))/r17)))/(8.*PI*(1 - nu)) + (burgers_vector.c[2]*(y*r11*(-(((-z + r8*r10)*(x/r8 - r11))/(r8*r13)) + (x*r10)/((r1 + r2 + r3)*(r8 - z*r10 - x*r11)) - (x*(-z + r8*r10))/(r5*(r8 - z*r10 - x*r11)) + ((2*depth + z + r9*r10)*(x/r9 - r11))/(r9*r17) - (x*r10)/((r1 + r2 + r4)*(r9 + (2*depth + z)*r10 - x*r11)) + (x*(2*depth + z + r9*r10))/(r6*(r9 + (2*depth + z)*r10 - x*r11))) + 2*(1 - nu)*(-((y*r10)/(r15*(1 + r2/r15))) + (y*r10)/(r18*(1 + r2/r18)) + (-((y*r8*r11*(2*x*r10 - z*r11))/r19) + (x*y*r11)/(r8*(r2*r10 + x*(x*r10 - z*r11))))/(1 + (r2*(r1 + r2 + r3)*r16)/r19) - (-((y*r9*r11*(2*x*r10 + (2*depth + z)*r11))/r20) + (x*y*r11)/(r9*(r2*r10 + x*(x*r10 + (2*depth + z)*r11))))/(1 + (r2*(r1 + r2 + r4)*r16)/r20)) + 2*((y*(depth + z)*r11*((-2*depth*x*(2*depth + z))/r7 - ((depth/r9 + r10)*(2*depth + z + r9*r10)*(x/r9 - r11))/r17 + (x*r10*(depth/r9 + r10))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) - (depth*x*(2*depth + z + r9*r10))/(r6*(r9 + (2*depth + z)*r10 - x*r11))))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) - (y*(depth + z)*(x/r9 - r11)*r11*(1 + (depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r9 + r10)*(2*depth + z + r9*r10))/(r9 + (2*depth + z)*r10 - x*r11)))/(r9*r17) - (x*y*(depth + z)*r11*(1 + (depth*(2*depth + z))/(r1 + r2 + r4) + ((depth/r9 + r10)*(2*depth + z + r9*r10))/(r9 + (2*depth + z)*r10 - x*r11)))/(r6*(r9 + (2*depth + z)*r10 - x*r11)) + 2*(1 - nu)*(y/(r1 + r2) - (y*(depth/r9 + r10)*(x/r9 - r11)*r11)/r17 - (depth*x*y*r11)/(r6*(r9 + (2*depth + z)*r10 - x*r11)) - (y*r10)/(r18*(1 + r2/r18)) + (-((y*r9*r11*(2*x*r10 + (2*depth + z)*r11))/r20) + (x*y*r11)/(r9*(r2*r10 + x*(x*r10 + (2*depth + z)*r11))))/(1 + (r2*(r1 + r2 + r4)*r16)/r20)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(y*(-(x/r5) + x/r6 - r10*(-(((-z + r8*r10)*(x/r8 - r11))/(r8*r13)) + (x*r10)/((r1 + r2 + r3)*(r8 - z*r10 - x*r11)) - (x*(-z + r8*r10))/(r5*(r8 - z*r10 - x*r11)) + ((2*depth + z + r9*r10)*(x/r9 - r11))/(r9*r17) - (x*r10)/((r1 + r2 + r4)*(r9 + (2*depth + z)*r10 - x*r11)) + (x*(2*depth + z + r9*r10))/(r6*(r9 + (2*depth + z)*r10 - x*r11)))) + 2*((y*(depth + z)*((-2*depth*x)/r7 - (2*nu*x)/(r9*r14)))/r9 - (x*y*(depth + z)*(depth/(r1 + r2 + r4) + (2*nu)/(2*depth + z + r9)))/r6 + (y*(depth + z)*r10*((2*depth*x*(2*depth + z))/r7 + ((depth/r9 + r10)*(2*depth + z + r9*r10)*(x/r9 - r11))/r17 - (x*r10*(depth/r9 + r10))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) + (depth*x*(2*depth + z + r9*r10))/(r6*(r9 + (2*depth + z)*r10 - x*r11))))/(r9*(r9 + (2*depth + z)*r10 - x*r11)) - (y*(depth + z)*r10*(x/r9 - r11)*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r4) - ((depth/r9 + r10)*(2*depth + z + r9*r10))/(r9 + (2*depth + z)*r10 - x*r11)))/(r9*r17) - (x*y*(depth + z)*r10*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r4) - ((depth/r9 + r10)*(2*depth + z + r9*r10))/(r9 + (2*depth + z)*r10 - x*r11)))/(r6*(r9 + (2*depth + z)*r10 - x*r11)) + 2*(1 - nu)*(-((x*y*(2*nu + depth/r9))/(r9*r14)) - (depth*x*y)/(r6*(2*depth + z + r9)) + (y*r10*(depth/r9 + r10)*(x/r9 - r11))/r17 + (depth*x*y*r10)/(r6*(r9 + (2*depth + z)*r10 - x*r11)) + (1 - 2*nu)*r12*(y/(r1 + r2) - (y*r10)/(r18*(1 + r2/r18)) + (-((y*r9*r11*(2*x*r10 + (2*depth + z)*r11))/r20) + (x*y*r11)/(r9*(r2*r10 + x*(x*r10 + (2*depth + z)*r11))))/(1 + (r2*(r1 + r2 + r4)*r16)/r20))))))/(8.*PI*(1 - nu));
        return uzx;
    }
    double uzycalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(r1 + r2 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r1 + r2 + r5,1.5);
        double r7 = pow(r1 + r2 + r5,2);
        double r8 = sqrt(r1 + r2 + r3);
        double r9 = cos(kappa);
        double r10 = sin(kappa);
        double r11 = sqrt(r1 + r2 + r5);
        double r12 = 1/tan(kappa);
        double r13 = pow(r8 - z*r9 - x*r10,2);
        double r14 = pow(2*depth + z + r11,2);
        double r15 = pow(x*r9 - z*r10,2);
        double r16 = pow(r10,2);
        double r17 = pow(r11 + (2*depth + z)*r9 - x*r10,2);
        double r18 = pow(x*r9 + (2*depth + z)*r10,2);
        double r19 = pow(r2*r9 + x*(x*r9 - z*r10),2);
        double r20 = pow(r2*r9 + x*(x*r9 + (2*depth + z)*r10),2);
        double uzy=(burgers_vector.c[1]*(-(x*(-(y/r4) + y/r6)) - (y*(-z + r8*r9)*(x*r9 - z*r10))/((r1 + r2 + r3)*r13) + (y*r9*(x*r9 - z*r10))/((r1 + r2 + r3)*(r8 - z*r9 - x*r10)) - (y*(-z + r8*r9)*(x*r9 - z*r10))/(r4*(r8 - z*r9 - x*r10)) + (y*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/((r1 + r2 + r5)*r17) - (y*r9*(x*r9 + (2*depth + z)*r10))/((r1 + r2 + r5)*(r11 + (2*depth + z)*r9 - x*r10)) + (y*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + (-1 + 2*nu)*r10*(y/(r8*(r8 - z*r9 - x*r10)) - y/(r11*(r11 + (2*depth + z)*r9 - x*r10))) + 2*((2*(1 - nu)*x*y*(2*nu + depth/r11))/(r11*r14) + (2*depth*(1 - nu)*x*y)/(r6*(2*depth + z + r11)) + ((depth + z)*((2*depth*x*y)/r7 + (2*nu*x*y)/(r11*r14)))/r11 - (y*(depth + z)*(-((depth*x)/(r1 + r2 + r5)) - (2*nu*x)/(2*depth + z + r11) + (1 - 2*nu)*r12))/r6 - (2*(1 - nu)*y*(depth/r11 + r9)*(x*r9 + (2*depth + z)*r10))/(r11*r17) - (2*depth*(1 - nu)*y*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) - 2*(1 - 2*nu)*(1 - nu)*r12*(y/(r11*(2*depth + z + r11)) - (y*r9)/(r11*(r11 + (2*depth + z)*r9 - x*r10))) - ((depth + z)*(((2*depth + z + r11*r9)*r12*((y*(2*depth + z + r11*r9))/(r11*r17) - (y*r9)/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r11 + (y*r9*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/(r1 + r2 + r5) - (y*(2*depth + z + r11*r9)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r6 + (depth*((2*y*(2*depth + z)*(x*r9 + (2*depth + z)*r10))/r7 + (y*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/((r1 + r2 + r5)*r17) - (y*r9*(x*r9 + (2*depth + z)*r10))/((r1 + r2 + r5)*(r11 + (2*depth + z)*r9 - x*r10)) + (y*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/r11 - (depth*y*(r10 - ((2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r1 + r2 + r5) - ((2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r6))/(r11 + (2*depth + z)*r9 - x*r10) + (y*(depth + z)*(r9*r10 + ((2*depth + z + r11*r9)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r11 + (depth*(r10 - ((2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r1 + r2 + r5) - ((2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r11))/(r11*r17))))/(8.*PI*(1 - nu)) + (burgers_vector.c[2]*(y*r10*(-((y*(-z + r8*r9))/((r1 + r2 + r3)*r13)) + (y*r9)/((r1 + r2 + r3)*(r8 - z*r9 - x*r10)) - (y*(-z + r8*r9))/(r4*(r8 - z*r9 - x*r10)) + (y*(2*depth + z + r11*r9))/((r1 + r2 + r5)*r17) - (y*r9)/((r1 + r2 + r5)*(r11 + (2*depth + z)*r9 - x*r10)) + (y*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))) + r10*((-z + r8*r9)/(r8*(r8 - z*r9 - x*r10)) - (2*depth + z + r11*r9)/(r11*(r11 + (2*depth + z)*r9 - x*r10))) + 2*(1 - nu)*(1/((x*r9 - z*r10)*(1 + r2/r15)) - 1/((x*r9 + (2*depth + z)*r10)*(1 + r2/r18)) + ((-2*r2*r8*r9*r10)/r19 + (r2*r10)/(r8*(r2*r9 + x*(x*r9 - z*r10))) + (r8*r10)/(r2*r9 + x*(x*r9 - z*r10)))/(1 + (r2*(r1 + r2 + r3)*r16)/r19) - ((-2*r2*r11*r9*r10)/r20 + (r2*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))) + (r11*r10)/(r2*r9 + x*(x*r9 + (2*depth + z)*r10)))/(1 + (r2*(r1 + r2 + r5)*r16)/r20)) + 2*((y*(depth + z)*r10*((-2*depth*y*(2*depth + z))/r7 - (y*(depth/r11 + r9)*(2*depth + z + r11*r9))/(r11*r17) + (y*r9*(depth/r11 + r9))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - (depth*y*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - (r2*(depth + z)*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/((r1 + r2 + r5)*r17) - (r2*(depth + z)*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + ((depth + z)*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + 2*(1 - nu)*(-(x/(r1 + r2)) - (r2*(depth/r11 + r9)*r10)/(r11*r17) - (depth*r2*r10)/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + ((depth/r11 + r9)*r10)/(r11 + (2*depth + z)*r9 - x*r10) + 1/((x*r9 + (2*depth + z)*r10)*(1 + r2/r18)) + ((-2*r2*r11*r9*r10)/r20 + (r2*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))) + (r11*r10)/(r2*r9 + x*(x*r9 + (2*depth + z)*r10)))/(1 + (r2*(r1 + r2 + r5)*r16)/r20)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(1/r8 - 1/r11 - r9*((-z + r8*r9)/(r8*(r8 - z*r9 - x*r10)) - (2*depth + z + r11*r9)/(r11*(r11 + (2*depth + z)*r9 - x*r10))) + y*(-(y/r4) + y/r6 - r9*(-((y*(-z + r8*r9))/((r1 + r2 + r3)*r13)) + (y*r9)/((r1 + r2 + r3)*(r8 - z*r9 - x*r10)) - (y*(-z + r8*r9))/(r4*(r8 - z*r9 - x*r10)) + (y*(2*depth + z + r11*r9))/((r1 + r2 + r5)*r17) - (y*r9)/((r1 + r2 + r5)*(r11 + (2*depth + z)*r9 - x*r10)) + (y*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10)))) + 2*((y*(depth + z)*((-2*depth*y)/r7 - (2*nu*y)/(r11*r14)))/r11 - (r2*(depth + z)*(depth/(r1 + r2 + r5) + (2*nu)/(2*depth + z + r11)))/r6 + ((depth + z)*(depth/(r1 + r2 + r5) + (2*nu)/(2*depth + z + r11)))/r11 + (y*(depth + z)*r9*((2*depth*y*(2*depth + z))/r7 + (y*(depth/r11 + r9)*(2*depth + z + r11*r9))/(r11*r17) - (y*r9*(depth/r11 + r9))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + (depth*y*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - (r2*(depth + z)*r9*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/((r1 + r2 + r5)*r17) - (r2*(depth + z)*r9*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + ((depth + z)*r9*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + 2*(1 - nu)*(-((r2*(2*nu + depth/r11))/(r11*r14)) - (depth*r2)/(r6*(2*depth + z + r11)) + (2*nu + depth/r11)/(2*depth + z + r11) + (r2*r9*(depth/r11 + r9))/(r11*r17) + (depth*r2*r9)/(r6*(r11 + (2*depth + z)*r9 - x*r10)) - (r9*(depth/r11 + r9))/(r11 + (2*depth + z)*r9 - x*r10) + (1 - 2*nu)*r12*(-(x/(r1 + r2)) + 1/((x*r9 + (2*depth + z)*r10)*(1 + r2/r18)) + ((-2*r2*r11*r9*r10)/r20 + (r2*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))) + (r11*r10)/(r2*r9 + x*(x*r9 + (2*depth + z)*r10)))/(1 + (r2*(r1 + r2 + r5)*r16)/r20))))))/(8.*PI*(1 - nu));
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override
    {
        double r1 = pow(x,2);
        double r2 = pow(y,2);
        double r3 = pow(z,2);
        double r4 = pow(r1 + r2 + r3,1.5);
        double r5 = pow(2*depth + z,2);
        double r6 = pow(r1 + r2 + r5,1.5);
        double r7 = pow(r1 + r2 + r5,2);
        double r8 = sqrt(r1 + r2 + r3);
        double r9 = cos(kappa);
        double r10 = sin(kappa);
        double r11 = sqrt(r1 + r2 + r5);
        double r12 = 1/tan(kappa);
        double r13 = pow(r8 - z*r9 - x*r10,2);
        double r14 = pow(2*depth + z + r11,2);
        double r15 = pow(x*r9 - z*r10,2);
        double r16 = pow(r10,2);
        double r17 = pow(r11 + (2*depth + z)*r9 - x*r10,2);
        double r18 = pow(x*r9 + (2*depth + z)*r10,2);
        double r19 = pow(r2*r9 + x*(x*r9 - z*r10),2);
        double r20 = pow(r2*r9 + x*(x*r9 + (2*depth + z)*r10),2);
        double uzz=(burgers_vector.c[1]*(-(x*(-(z/r4) + (2*depth + z)/r6)) - ((-z + r8*r9)*r10)/(r8*(r8 - z*r9 - x*r10)) - ((2*depth + z + r11*r9)*r10)/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - ((z/r8 - r9)*(-z + r8*r9)*(x*r9 - z*r10))/(r8*r13) + ((-1 + (z*r9)/r8)*(x*r9 - z*r10))/(r8*(r8 - z*r9 - x*r10)) - (z*(-z + r8*r9)*(x*r9 - z*r10))/(r4*(r8 - z*r9 - x*r10)) + (((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*r17) - ((1 + ((2*depth + z)*r9)/r11)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + ((2*depth + z)*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + (-1 + 2*nu)*r10*((z/r8 - r9)/(r8 - z*r9 - x*r10) - ((2*depth + z)/r11 + r9)/(r11 + (2*depth + z)*r9 - x*r10)) + 2*((2*(1 - nu)*x*(2*nu + depth/r11)*(1 + (2*depth + z)/r11))/r14 + (2*depth*(1 - nu)*x*(2*depth + z))/(r6*(2*depth + z + r11)) + ((depth + z)*((2*depth*x*(2*depth + z))/r7 + (2*nu*x*(1 + (2*depth + z)/r11))/r14))/r11 - ((depth + z)*(2*depth + z)*(-((depth*x)/(r1 + r2 + r5)) - (2*nu*x)/(2*depth + z + r11) + (1 - 2*nu)*r12))/r6 + (-((depth*x)/(r1 + r2 + r5)) - (2*nu*x)/(2*depth + z + r11) + (1 - 2*nu)*r12)/r11 + (2*(1 - nu)*(depth/r11 + r9)*r10)/(r11 + (2*depth + z)*r9 - x*r10) - (2*(1 - nu)*(depth/r11 + r9)*((2*depth + z)/r11 + r9)*(x*r9 + (2*depth + z)*r10))/r17 - (2*depth*(1 - nu)*(2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) - 2*(1 - 2*nu)*(1 - nu)*r12*((1 + (2*depth + z)/r11)/(2*depth + z + r11) - (r9*((2*depth + z)/r11 + r9))/(r11 + (2*depth + z)*r9 - x*r10)) - ((depth + z)*(((2*depth + z + r11*r9)*r12*((((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9))/r17 - (1 + ((2*depth + z)*r9)/r11)/(r11 + (2*depth + z)*r9 - x*r10)))/r11 + ((1 + ((2*depth + z)*r9)/r11)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r11 - ((2*depth + z)*(2*depth + z + r11*r9)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r6 + (depth*(-(((2*depth + z)*r10)/(r1 + r2 + r5)) - ((2*depth + z + r11*r9)*r10)/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + (2*r5*(x*r9 + (2*depth + z)*r10))/r7 - (x*r9 + (2*depth + z)*r10)/(r1 + r2 + r5) + (((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*r17) - ((1 + ((2*depth + z)*r9)/r11)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + ((2*depth + z)*(2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/r11 - (depth*(2*depth + z)*(r10 - ((2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r1 + r2 + r5) - ((2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r6))/(r11 + (2*depth + z)*r9 - x*r10) + ((depth + z)*((2*depth + z)/r11 + r9)*(r9*r10 + ((2*depth + z + r11*r9)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r11 + (depth*(r10 - ((2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r1 + r2 + r5) - ((2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r11))/r17 - (r9*r10 + ((2*depth + z + r11*r9)*r12*(2*(1 - nu)*r9 - (2*depth + z + r11*r9)/(r11 + (2*depth + z)*r9 - x*r10)))/r11 + (depth*(r10 - ((2*depth + z)*(x*r9 + (2*depth + z)*r10))/(r1 + r2 + r5) - ((2*depth + z + r11*r9)*(x*r9 + (2*depth + z)*r10))/(r11*(r11 + (2*depth + z)*r9 - x*r10))))/r11)/(r11 + (2*depth + z)*r9 - x*r10))))/(8.*PI*(1 - nu)) + (burgers_vector.c[2]*(y*r10*(-(((z/r8 - r9)*(-z + r8*r9))/(r8*r13)) + (-1 + (z*r9)/r8)/(r8*(r8 - z*r9 - x*r10)) - (z*(-z + r8*r9))/(r4*(r8 - z*r9 - x*r10)) + (((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9))/(r11*r17) - (1 + ((2*depth + z)*r9)/r11)/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + ((2*depth + z)*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))) + 2*(1 - nu)*((y*r10)/(r15*(1 + r2/r15)) + (y*r10)/(r18*(1 + r2/r18)) + ((x*y*r8*r16)/r19 + (y*z*r10)/(r8*(r2*r9 + x*(x*r9 - z*r10))))/(1 + (r2*(r1 + r2 + r3)*r16)/r19) - (-((x*y*r11*r16)/r20) + (y*(2*depth + z)*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))))/(1 + (r2*(r1 + r2 + r5)*r16)/r20)) + 2*((y*(depth + z)*r10*((-2*depth*r5)/r7 + depth/(r1 + r2 + r5) - ((depth/r11 + r9)*((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9))/r17 + ((depth/r11 + r9)*(1 + ((2*depth + z)*r9)/r11))/(r11 + (2*depth + z)*r9 - x*r10) - (depth*(2*depth + z)*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - (y*(depth + z)*((2*depth + z)/r11 + r9)*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*r17) - (y*(depth + z)*(2*depth + z)*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + (y*r10*(1 + (depth*(2*depth + z))/(r1 + r2 + r5) + ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + 2*(1 - nu)*(-((y*(depth/r11 + r9)*((2*depth + z)/r11 + r9)*r10)/r17) - (depth*y*(2*depth + z)*r10)/(r6*(r11 + (2*depth + z)*r9 - x*r10)) - (y*r10)/(r18*(1 + r2/r18)) + (-((x*y*r11*r16)/r20) + (y*(2*depth + z)*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))))/(1 + (r2*(r1 + r2 + r5)*r16)/r20)))))/(8.*PI*(1 - nu)) + (burgers_vector.c[0]*(y*(-(z/r4) + (2*depth + z)/r6 - r9*(-(((z/r8 - r9)*(-z + r8*r9))/(r8*r13)) + (-1 + (z*r9)/r8)/(r8*(r8 - z*r9 - x*r10)) - (z*(-z + r8*r9))/(r4*(r8 - z*r9 - x*r10)) + (((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9))/(r11*r17) - (1 + ((2*depth + z)*r9)/r11)/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + ((2*depth + z)*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10)))) + 2*((y*(depth + z)*((-2*depth*(2*depth + z))/r7 - (2*nu*(1 + (2*depth + z)/r11))/r14))/r11 - (y*(depth + z)*(2*depth + z)*(depth/(r1 + r2 + r5) + (2*nu)/(2*depth + z + r11)))/r6 + (y*(depth/(r1 + r2 + r5) + (2*nu)/(2*depth + z + r11)))/r11 + (y*(depth + z)*r9*((2*depth*r5)/r7 - depth/(r1 + r2 + r5) + ((depth/r11 + r9)*((2*depth + z)/r11 + r9)*(2*depth + z + r11*r9))/r17 - ((depth/r11 + r9)*(1 + ((2*depth + z)*r9)/r11))/(r11 + (2*depth + z)*r9 - x*r10) + (depth*(2*depth + z)*(2*depth + z + r11*r9))/(r6*(r11 + (2*depth + z)*r9 - x*r10))))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) - (y*(depth + z)*r9*((2*depth + z)/r11 + r9)*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*r17) - (y*(depth + z)*(2*depth + z)*r9*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + (y*r9*(1 - 2*nu - (depth*(2*depth + z))/(r1 + r2 + r5) - ((depth/r11 + r9)*(2*depth + z + r11*r9))/(r11 + (2*depth + z)*r9 - x*r10)))/(r11*(r11 + (2*depth + z)*r9 - x*r10)) + 2*(1 - nu)*(-((y*(2*nu + depth/r11)*(1 + (2*depth + z)/r11))/r14) - (depth*y*(2*depth + z))/(r6*(2*depth + z + r11)) + (y*r9*(depth/r11 + r9)*((2*depth + z)/r11 + r9))/r17 + (depth*y*(2*depth + z)*r9)/(r6*(r11 + (2*depth + z)*r9 - x*r10)) + (1 - 2*nu)*r12*(-((y*r10)/(r18*(1 + r2/r18))) + (-((x*y*r11*r16)/r20) + (y*(2*depth + z)*r10)/(r11*(r2*r9 + x*(x*r9 + (2*depth + z)*r10))))/(1 + (r2*(r1 + r2 + r5)*r16)/r20))))))/(8.*PI*(1 - nu));
        return uzz;
    }
};

class DisplacmentGradientSystemReplace: public virtual DisplacmentGradient{
public:
    Rotation_Matrix_Z *rotate;
    Rotation_Matrix_Z_Vector *rotate_vector;
    DisplacmentGradient *distortion_gradients;
    double phi;
    DisplacmentGradientSystemReplace(Vector<double> burgers_vector, double depth, double nu, double phi, double kappa, DisplacmentGradient& distorsion_gradients)//mock
    {
        rotate = new Rotation_Matrix_Z(phi);
        rotate_vector = new Rotation_Matrix_Z_Vector(phi);
        distortion_gradients=&distorsion_gradients;
        this->phi=phi;
        this->depth=depth;
        this->nu=nu;
        this->kappa=kappa;
        rotate_vector->Basis(burgers_vector);
        this->burgers_vector=rotate_vector->vector;
        distortion_gradients->burgers_vector=this->burgers_vector;
    }
    double uxxcalc(double x, double y, double z) override
    {
        rotate->Basis(x,y,z);
        double uxx=distortion_gradients->uxxcalc(rotate->xc, rotate->yc, rotate->zc)*cos(phi)*cos(phi) - distortion_gradients->uxycalc(rotate->xc, rotate->yc, rotate->zc)*cos(phi)*sin(phi) - distortion_gradients->uyxcalc(rotate->xc, rotate->yc, rotate->zc)*cos(phi)*sin(phi) + distortion_gradients->uyycalc(rotate->xc, rotate->yc, rotate->zc)*sin(phi)*sin(phi);
        return uxx;
    }
    double uxycalc(double x, double y, double z) override
    {
        double uxy=distortion_gradients->uxycalc (rotate -> xc, rotate -> yc, rotate -> zc)*pow(cos(phi), 2) + distortion_gradients->uxxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos(phi)*sin(phi) - distortion_gradients->uyycalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos(phi)*sin(phi) -distortion_gradients->uyxcalc(rotate -> xc, rotate -> yc, rotate -> zc)*pow(sin(phi), 2);
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override
    {
        double uxz=distortion_gradients->uxzcalc(rotate -> xc, rotate -> yc, rotate -> zc)*cos(phi) - distortion_gradients->uyzcalc(rotate -> xc, rotate -> yc, rotate -> zc)*sin (phi);
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override
    {
        double uyx=distortion_gradients->uyxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*pow (cos (phi), 2) + distortion_gradients ->uxxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi)*sin (phi) - distortion_gradients ->uyycalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi)*sin (phi) - distortion_gradients ->uxycalc (rotate -> xc, rotate -> yc, rotate -> zc)*pow (sin (phi), 2);
        return uyx;
    }
    double uyycalc(double x, double y, double z) override
    {
        double uyy=distortion_gradients ->uyycalc (rotate -> xc, rotate -> yc, rotate -> zc)*pow (cos (phi), 2) + distortion_gradients ->uxycalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi)*sin (phi) + distortion_gradients ->uyxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi)*sin (phi) + distortion_gradients ->uxxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*pow (sin (phi), 2);
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override
    {
        double uyz=distortion_gradients ->uyzcalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi) +distortion_gradients ->uxzcalc (rotate -> xc, rotate -> yc, rotate -> zc)*sin (phi);
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override
    {
        double uzx=distortion_gradients ->uzxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi) -distortion_gradients ->uzycalc (rotate -> xc, rotate -> yc, rotate -> zc)*sin (phi);
        return uzx;
    }
    double uzycalc(double x, double y, double z) override
    {
        double uzy=distortion_gradients ->uzycalc (rotate -> xc, rotate -> yc, rotate -> zc)*cos (phi) +distortion_gradients ->uzxcalc (rotate -> xc, rotate -> yc, rotate -> zc)*sin (phi);
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override
    {
        double uzz=distortion_gradients ->uzzcalc (rotate -> xc, rotate -> yc, rotate -> zc);
        return uzz;
    }
};

class InitialisationGeometry{
public:
    Vector<double> x_vector,y_vector,z_vector, glide_plane_vector, b_vector;
    std::vector <double> kappa, phi;
    std::vector <Vector<double>> exit_point;
    double dislocation_depth;
    InitialisationGeometry(Vector<double> diffraction_vector,Vector<double> normal_vector,std::vector <Vector<double>> direction_vector_list, Vector<double> burgers_vector, int nubmer_of_dislocation, double dislocation_depth, std::vector <double> segment_lenght){
        x_vector=Vector_Normalization<double, double>(diffraction_vector);
        z_vector=Vector_Normalization<double, double>(Vector_Inverse(normal_vector));
        y_vector=Vector_Multiplication<double, double>(z_vector, x_vector);
        try{
        Glide_Plane(direction_vector_list);
        }
        catch (const char *str){
            std::cout << str;
        }
        burgers_vector.c[0]=burgers_vector.c[0]*nubmer_of_dislocation;
        burgers_vector.c[1]=burgers_vector.c[1]*nubmer_of_dislocation;
        burgers_vector.c[2]=burgers_vector.c[2]*nubmer_of_dislocation;
        Burgers_Vector_Into_System_Coordinates(burgers_vector);
        this->dislocation_depth=dislocation_depth;
        Angle_Between_Dislocation_Axis_And_Calculation_Axis(direction_vector_list);
    }
    void Glide_Plane(std::vector <Vector<double>> direction_vector_list){//untested
        switch (direction_vector_list.size()) {
            case 0:
                std::cerr<<"Can't calculate a glide plane, list is empty!"<<std::endl;
                break;
            case 1:
                std::cerr<<"Can't calculate a glide plane, list has only one vector!"<<std::endl;
                break;
            default:
                Glide_Plane_Calculate(direction_vector_list);
                break;
        }
    }
    void Glide_Plane_Calculate(std::vector <Vector<double>> direction_vector_list){//untested
        Vector<double> comparison_vector;
        for(int i=1;i<direction_vector_list.size();i++){
            glide_plane_vector=Vector_Multiplication<double, double>(direction_vector_list[i-1], direction_vector_list[i]);
            if (i!=1 and comparison_vector!=glide_plane_vector){
                throw "This vectors haven't overall glide plane";
            }
            comparison_vector=glide_plane_vector;
        }
    }
    void Exit_Point_Coordinate(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght){//untested
        Vector<double> exit_point_calculation;
        for(int i=0;i<direction_vector_list.size();i++)
        {
            exit_point_calculation.c[2]=-dislocation_depth;
            if(kappa[i]>=(PI_2-0.005*PI_2) && kappa[i]<=(PI_2+0.005*PI_2)){
                exit_point_calculation.c[0]=(segment_lenght[i])*cos(phi[i]);
                exit_point_calculation.c[1]=(segment_lenght[i])*sin(phi[i]);
            }
            else{
                exit_point_calculation.c[0]=(dislocation_depth/tan(kappa[i])+segment_lenght[i])*cos(phi[i]);
                exit_point_calculation.c[1]=(dislocation_depth/tan(kappa[i])+segment_lenght[i])*sin(phi[i]);
            }
            exit_point.push_back(exit_point_calculation);
        }
    }
    void Burgers_Vector_Into_System_Coordinates(Vector<double> b_initial){//not tested
        b_vector.c[0]=-((-(b_initial.c[2]*y_vector.c[1]*z_vector.c[0]) + b_initial.c[1]*y_vector.c[2]*z_vector.c[0] + b_initial.c[2]*y_vector.c[0]*z_vector.c[1] - b_initial.c[0]*y_vector.c[2]*z_vector.c[1] - b_initial.c[1]*y_vector.c[0]*z_vector.c[2] + b_initial.c[0]*y_vector.c[1]*z_vector.c[2])/ (x_vector.c[2]*y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[1]=-((b_initial.c[2]*x_vector.c[1]*z_vector.c[0] - b_initial.c[1]*x_vector.c[2]*z_vector.c[0] - b_initial.c[2]*x_vector.c[0]*z_vector.c[1] + b_initial.c[0]*x_vector.c[2]*z_vector.c[1] + b_initial.c[1]*x_vector.c[0]*z_vector.c[2] - b_initial.c[0]*x_vector.c[1]*z_vector.c[2])/(x_vector.c[2]* y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[2]=-((b_initial.c[2]*x_vector.c[1]*y_vector.c[0] - b_initial.c[1]*x_vector.c[2]*y_vector.c[0] - b_initial.c[2]*x_vector.c[0]*y_vector.c[1] + b_initial.c[0]*x_vector.c[2]*y_vector.c[1] + b_initial.c[1]*x_vector.c[0]*y_vector.c[2] - b_initial.c[0]*x_vector.c[1]* y_vector.c[2])/(-(x_vector.c[2]*y_vector.c[1]*z_vector.c[0]) + x_vector.c[1]*y_vector.c[2]*z_vector.c[0] + x_vector.c[2]*y_vector.c[0]*z_vector.c[1] - x_vector.c[0]*y_vector.c[2]*z_vector.c[1] - x_vector.c[1]*y_vector.c[0]*z_vector.c[2] + x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
    }
    void Angle_Between_Dislocation_Axis_And_Calculation_Axis(std::vector <Vector<double>> direction_vector_list){//untested
        double kappa_calculate,phi_calculate;
        for(int i=0;i<direction_vector_list.size();i++)
        {
            Vector <double> vector_normal_tayz=Vector_Multiplication<double, double>(direction_vector_list[i], z_vector);
            Vector <double> vector_normal_xz=Vector_Multiplication<double, double>(x_vector, z_vector);
            if(Vector_Absolute_Value(vector_normal_tayz)==0 || Vector_Absolute_Value(vector_normal_tayz)==0){
                phi_calculate=0;
            }
            else{
                phi_calculate=acos(Scalar_Multiplication(vector_normal_tayz, vector_normal_xz)/(Vector_Absolute_Value(vector_normal_xz)*Vector_Absolute_Value(vector_normal_tayz)));
            }
            kappa_calculate=acos(Scalar_Multiplication(z_vector, direction_vector_list[i])/(Vector_Absolute_Value(z_vector)*Vector_Absolute_Value(direction_vector_list[i])));
            if(kappa_calculate>PI_2){
                kappa_calculate=kappa_calculate-PI_2;
            }
            kappa.push_back(kappa_calculate);
            phi.push_back(phi_calculate);
        }
    }
};

int main(int argc, const char * argv[]) {
    //Input data
    //Diffraction vector:
    Vector<double> diffraction_vector;
    diffraction_vector.c[0]=1;
    diffraction_vector.c[1]=0;
    diffraction_vector.c[2]=0;
    //Normal to surface vector:
    Vector<double> normal_vector;
    normal_vector.c[0]=0;
    normal_vector.c[1]=0;
    normal_vector.c[2]=-1;
    //Vector of tay-vectors:
    std::vector <Vector<double>> tay_vector;
    Vector<double> tay_1;
    tay_1.c[0]=-1;
    tay_1.c[1]=0;
    tay_1.c[2]=-1;
    tay_vector.push_back(tay_1);
    Vector<double> tay_2;
    tay_2.c[0]=1;
    tay_2.c[1]=0;
    tay_2.c[2]=0;
    tay_vector.push_back(tay_2);
    //Burger's vectorof dislocation:
    Vector<double> burgers_vector;
    burgers_vector.c[0]=1;
    burgers_vector.c[1]=1;
    burgers_vector.c[2]=1;
    //Number of dislocation in one place:
    int number_of_dislocation=1;
    //Dislocation depth:
    double dislocation_depth=100;
    //Segment lengt (thi first and second position must be zeros! If the dislocation have more, then to segments, then add into third position a value of mid segment lenght)
    std::vector<double> segment_lenght;
    segment_lenght.push_back(0);
    segment_lenght.push_back(0);
    //Start of initialization
    InitialisationGeometry *init_obj = new InitialisationGeometry(diffraction_vector, normal_vector,tay_vector, burgers_vector, number_of_dislocation,dislocation_depth,segment_lenght);
    AngularDislocation *test = new AngularDislocation();
    test->burgers_vector=burgers_vector;
    test->depth=dislocation_depth;
    test->kappa=1.0;
    test->nu=0.4;
    std::cout<<test->uxxcalc(5, 5, 5)<<std::endl;
    std::cout<<test->uxycalc(5, 5, 5)<<std::endl;
    std::cout<<test->uxzcalc(5, 5, 5)<<std::endl;
    std::cout<<test->uyxcalc(5, 5, 5)<<std::endl;
    std::cout<<test->uyycalc(5, 5, 5)<<std::endl;
    std::cout<<test->uyzcalc(5, 5, 5)<<std::endl;
    std::cout<<test->uzxcalc(5, 5, 5)<<std::endl;
    std::cout<<test->uzycalc(5, 5, 5)<<std::endl;
    std::cout<<test->uzzcalc(5, 5, 5)<<std::endl;
    return 0;
}
