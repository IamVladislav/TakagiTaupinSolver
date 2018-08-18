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


class DisplacmentGradient{
public:
    virtual ~DisplacmentGradient(){}
    double uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz;
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

class TestDisp: public DisplacmentGradient{
public:
    double uxxcalc(double x, double y, double z) override
    {
        uxx=10;
        return uxx;
    }
    double uxycalc(double x, double y, double z) override
    {
        uxy=5;
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override
    {
        uxz=5;
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override
    {
        uyx=5;
        return uyx;
    }
    double uyycalc(double x, double y, double z) override
    {
        uyy=5;
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override
    {
        uyy=5;
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override
    {
        uzx=5;
        return uzx;
    }
    double uzycalc(double x, double y, double z) override
    {
        uzy=5;
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override
    {
        uzz=5;
        return uzz;
    }
};

class DisplacmentGradientSystemReplace: public virtual DisplacmentGradient{
public:
    Rotation_Matrix_Z *rotate;
    DisplacmentGradient *displacment;
    double l;//mock
    DisplacmentGradientSystemReplace(double bu, double h, double nu, DisplacmentGradient& obj)//mock
    {
        rotate = new Rotation_Matrix_Z(l);
        displacment=&obj;
        l=5;//mock
    }
    double uxxcalc(double x, double y, double z) override
    {
        rotate->Basis(x,y,z);
        uxx=displacment->uxxcalc(rotate->xc, rotate->yc, rotate->zc)*cos(l)*cos(l) - displacment->uxycalc(rotate->xc, rotate->yc, rotate->zc)*cos(l)*sin(l) - displacment->uyxcalc(rotate->xc, rotate->yc, rotate->zc)*cos(l)*sin(l) + displacment->uyycalc(rotate->xc, rotate->yc, rotate->zc)*sin(l)*sin(l);
        return uxx;
    }
    double uxycalc(double x, double y, double z) override //mock
    {
        uxy=displacment->uxxcalc(rotate->xc, rotate->yc, rotate->zc);
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override //mock
    {
        uxz=5;
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override //mock
    {
        uyx=5;
        return uyx;
    }
    double uyycalc(double x, double y, double z) override //mock
    {
        uyy=5;
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override //mock
    {
        uyz=5;
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override //mock
    {
        uzx=5;
        return uzx;
    }
    double uzycalc(double x, double y, double z) override //mock
    {
        uzy=5;
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override //mock
    {
        uzz=5;
        return uzz;
    }
};

class InitializationGeometry{
public:
    Vector<double> x_vector,y_vector,z_vector;
    InitializationGeometry(Vector<int> diffraction_vector,Vector<int> normal_vector,std::vector <Vector<int>> direction_vector_list, Vector<int> burgers_vector, int nubmer_of_dislocation, double dislocation_depth){
        y_vector=Vector_Multiplication(z_vector, x_vector);
    }
    void Third_Vector_Of_System(){//mock
        y_vector.c[0]=1;
    }
    void Glide_Plane(){//mock
        
    }
    void Exit_Point_Coordinate(){//mock
        
    }
};

int main(int argc, const char * argv[]) {
    TestDisp p;
    p.uxxcalc(1, 1, 1);
    //std::cout << p.uxx << std::endl;
    Rotation_Matrix_Z t(10);
    //std::cout<<t.Test()<<std::endl;
    DisplacmentGradientSystemReplace obj(1, 1, 1, p);
    std::cout<<obj.uxxcalc(1, 1, 1)<<std::endl;
    std::vector <Vector<int>> test;
    Vector<int> test1;
    test1.c[0]=1;
    test1.c[1]=2;
    test1.c[2]=3;
    test.push_back(test1);
    std::cout<<test[0].c[1]<<std::endl;
    return 0;
}
