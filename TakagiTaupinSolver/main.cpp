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
    Vector<double> b;
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
    TestDisp(){
        b.c[0]=3;
        b.c[1]=3;
        b.c[2]=3;
    }
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

class InitialisationGeometry{
public:
    Vector<double> x_vector,y_vector,z_vector, glide_plane_vector, b_vector;
    InitialisationGeometry(Vector<int> diffraction_vector,Vector<int> normal_vector,std::vector <Vector<int>> direction_vector_list, Vector<int> burgers_vector, int nubmer_of_dislocation, double dislocation_depth){
        x_vector=Vector_Normalization<int, double>(diffraction_vector);
        z_vector=Vector_Normalization<int, double>(Vector_Inverse(normal_vector));
        y_vector=Vector_Multiplication<double, double>(z_vector, x_vector);
        try{
        Glide_Plane(direction_vector_list);
        }
        catch (const char *str){
            std::cout << str;
        }
        Burgers_Vector_Into_System_Coordinates(burgers_vector);
    }
    void Glide_Plane(std::vector <Vector<int>> direction_vector_list){//mock
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
    void Glide_Plane_Calculate(std::vector <Vector<int>> direction_vector_list){//mock
        Vector<double> comparison_vector;
        for(int i=1;i<direction_vector_list.size();i++){
            glide_plane_vector=Vector_Multiplication<int, double>(direction_vector_list[i-1], direction_vector_list[i]);
            if (i!=1 and comparison_vector!=glide_plane_vector){
                throw "This vectors haven't overall glide plane";
            }
            comparison_vector=glide_plane_vector;
        }
    }
    void Exit_Point_Coordinate(){//mock
        
    }
    void Burgers_Vector_Into_System_Coordinates(Vector<int> b_initial){//not tested
        b_vector.c[0]=-((-(b_initial.c[2]*y_vector.c[1]*z_vector.c[0]) + b_initial.c[1]*y_vector.c[2]*z_vector.c[0] + b_initial.c[2]*y_vector.c[0]*z_vector.c[1] - b_initial.c[0]*y_vector.c[2]*z_vector.c[1] - b_initial.c[1]*y_vector.c[0]*z_vector.c[2] + b_initial.c[0]*y_vector.c[1]*z_vector.c[2])/ (x_vector.c[2]*y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[1]=-((b_initial.c[2]*x_vector.c[1]*z_vector.c[0] - b_initial.c[1]*x_vector.c[2]*z_vector.c[0] - b_initial.c[2]*x_vector.c[0]*z_vector.c[1] + b_initial.c[0]*x_vector.c[2]*z_vector.c[1] + b_initial.c[1]*x_vector.c[0]*z_vector.c[2] - b_initial.c[0]*x_vector.c[1]*z_vector.c[2])/(x_vector.c[2]* y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[2]=-((b_initial.c[2]*x_vector.c[1]*y_vector.c[0] - b_initial.c[1]*x_vector.c[2]*y_vector.c[0] - b_initial.c[2]*x_vector.c[0]*y_vector.c[1] + b_initial.c[0]*x_vector.c[2]*y_vector.c[1] + b_initial.c[1]*x_vector.c[0]*y_vector.c[2] - b_initial.c[0]*x_vector.c[1]* y_vector.c[2])/(-(x_vector.c[2]*y_vector.c[1]*z_vector.c[0]) + x_vector.c[1]*y_vector.c[2]*z_vector.c[0] + x_vector.c[2]*y_vector.c[0]*z_vector.c[1] - x_vector.c[0]*y_vector.c[2]*z_vector.c[1] - x_vector.c[1]*y_vector.c[0]*z_vector.c[2] + x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
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
    Vector<int> diffraction_vector;
    diffraction_vector.c[0]=1;
    diffraction_vector.c[1]=0;
    diffraction_vector.c[2]=0;
    Vector<int> normal_vector;
    normal_vector.c[0]=0;
    normal_vector.c[1]=0;
    normal_vector.c[2]=-1;
    Vector<int> burgers_vector;
    burgers_vector.c[0]=1;
    burgers_vector.c[1]=1;
    burgers_vector.c[2]=1;
    test.push_back(diffraction_vector);
    test.push_back(normal_vector);
    test.push_back(burgers_vector);
    InitialisationGeometry p_test(diffraction_vector, normal_vector, test, burgers_vector, 2, 2);
    //std::cout << p_test.glide_plane_vector << std::endl;
    std::cout << p_test.b_vector << std::endl;
    return 0;
}
