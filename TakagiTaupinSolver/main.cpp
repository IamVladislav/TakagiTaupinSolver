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

class TestDisp: public DisplacmentGradient{
public:
    TestDisp(){
    }
    double uxxcalc(double x, double y, double z) override
    {
        double uxx=10;
        return uxx;
    }
    double uxycalc(double x, double y, double z) override
    {
        double uxy=5;
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override
    {
        double uxz=5;
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override
    {
        double uyx=5;
        return uyx;
    }
    double uyycalc(double x, double y, double z) override
    {
        double uyy=5;
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override
    {
        double uyz=5;
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override
    {
        double uzx=5;
        return uzx;
    }
    double uzycalc(double x, double y, double z) override
    {
        double uzy=5;
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override
    {
        double uzz=5;
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
    double uxycalc(double x, double y, double z) override //mock
    {
        double uxy=distortion_gradients->uxxcalc(rotate->xc, rotate->yc, rotate->zc);
        return uxy;
    }
    double uxzcalc(double x, double y, double z) override //mock
    {
        double uxz=5;
        return uxz;
    }
    double uyxcalc(double x, double y, double z) override //mock
    {
        double uyx=5;
        return uyx;
    }
    double uyycalc(double x, double y, double z) override //mock
    {
        double uyy=5;
        return uyy;
    }
    double uyzcalc(double x, double y, double z) override //mock
    {
        double uyz=5;
        return uyz;
    }
    double uzxcalc(double x, double y, double z) override //mock
    {
        double uzx=5;
        return uzx;
    }
    double uzycalc(double x, double y, double z) override //mock
    {
        double uzy=5;
        return uzy;
    }
    double uzzcalc(double x, double y, double z) override //mock
    {
        double uzz=5;
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
        Burgers_Vector_Into_System_Coordinates(burgers_vector);
        this->dislocation_depth=dislocation_depth;
        Angle_Between_Dislocation_Axis_And_Calculation_Axis(direction_vector_list);
    }
    void Glide_Plane(std::vector <Vector<double>> direction_vector_list){//not tested
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
    void Glide_Plane_Calculate(std::vector <Vector<double>> direction_vector_list){//not tested
        Vector<double> comparison_vector;
        for(int i=1;i<direction_vector_list.size();i++){
            glide_plane_vector=Vector_Multiplication<double, double>(direction_vector_list[i-1], direction_vector_list[i]);
            if (i!=1 and comparison_vector!=glide_plane_vector){
                throw "This vectors haven't overall glide plane";
            }
            comparison_vector=glide_plane_vector;
        }
    }
    void Exit_Point_Coordinate(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght){//TODO: Add a if-else for 1/0 exception
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
    void Angle_Between_Dislocation_Axis_And_Calculation_Axis(std::vector <Vector<double>> direction_vector_list){//TODO: find case and solution, if tay and phi more than pi, what;s happend?
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
    TestDisp p;
    p.uxxcalc(1, 1, 1);
    //std::cout << p.uxx << std::endl;
    Rotation_Matrix_Z t(10);
    //std::cout<<t.Test()<<std::endl;
    //std::cout<<obj.uxxcalc(1, 1, 1)<<std::endl;
    std::vector <Vector<double>> test;
    std::vector <double> segment_lenght;
    Vector<double> diffraction_vector;
    diffraction_vector.c[0]=1;
    diffraction_vector.c[1]=0;
    diffraction_vector.c[2]=0;
    Vector<double> normal_vector;
    normal_vector.c[0]=0;
    normal_vector.c[1]=0;
    normal_vector.c[2]=-1;
    Vector<double> burgers_vector;
    burgers_vector.c[0]=0;
    burgers_vector.c[1]=1;
    burgers_vector.c[2]=1;
    test.push_back(diffraction_vector);
    test.push_back(normal_vector);
    test.push_back(burgers_vector);
    segment_lenght.push_back(0);
    segment_lenght.push_back(0);
    segment_lenght.push_back(0);
    InitialisationGeometry p_test(diffraction_vector, normal_vector, test, burgers_vector, 2, 2, segment_lenght);
    DisplacmentGradientSystemReplace obj(burgers_vector, 1, 1, 1, 1, p);
    //std::cout << p_test.glide_plane_vector << std::endl;
    //std::cout<<std::endl;
    //std::cout << p_test.phi[2] << std::endl;
    return 0;
}
