//
//  Distorsion_Source.hpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#ifndef Distorsion_Source_hpp
#define Distorsion_Source_hpp

const double PI_2=1.5708;
const double PI=3.14159;
const double INF=10001;

#include <stdio.h>
#include <cmath>
#include <vector>
#include "Vector_And_Operations.hpp"
#include "System_Repclacement_Library.hpp"

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

class BeamDislocation: public DisplacmentGradient{//untested
public:
    double theta;
    Parallel_Shfit *parallel_shift;
    BeamDislocation(double nu, double kappa, Vector<double> exit_point_coordinate, double segment_lenght, Vector<double> burgers_vector);
    double uxxcalc(double x, double y, double z);
    double uxycalc(double x, double y, double z);
    double uxzcalc(double x, double y, double z);
    double uyxcalc(double x, double y, double z);
    double uyycalc(double x, double y, double z);
    double uyzcalc(double x, double y, double z);
    double uzxcalc(double x, double y, double z);
    double uzycalc(double x, double y, double z);
    double uzzcalc(double x, double y, double z);
};

class AngularDislocation: public DisplacmentGradient{
public:
    Parallel_Shfit *parallel_shift;
    AngularDislocation(double nu, double kappa, double segment_lenght, Vector<double> burgers_vector);
    double uxxcalc(double x, double y, double z);
    double uxycalc(double x, double y, double z);
    double uxzcalc(double x, double y, double z);
    double uyxcalc(double x, double y, double z);
    double uyycalc(double x, double y, double z);
    double uyzcalc(double x, double y, double z);
    double uzxcalc(double x, double y, double z);
    double uzycalc(double x, double y, double z);
    double uzzcalc(double x, double y, double z);
};

class DisplacmentGradientSystemReplace: public DisplacmentGradient{
public:
    Rotation_Matrix_Z *rotate;
    Rotation_Matrix_Z_Vector *rotate_vector;
    DisplacmentGradient *distortion_gradients;
    Vector<double> exit_point_coordinate;
    double phi, segment_lenght;
    DisplacmentGradientSystemReplace(double nu, double depth, double phi, double kappa, Vector<double> exit_point_coordinate, double segment_lenght, Vector<double> burgers_vector, int type);
    double uxxcalc(double x, double y, double z);
    double uxycalc(double x, double y, double z);
    double uxzcalc(double x, double y, double z);
    double uyxcalc(double x, double y, double z);
    double uyycalc(double x, double y, double z);
    double uyzcalc(double x, double y, double z);
    double uzxcalc(double x, double y, double z);
    double uzycalc(double x, double y, double z);
    double uzzcalc(double x, double y, double z);
};

class Model_Of_Polygonal_Dislocation: public DisplacmentGradient{
public:
    std::vector <DisplacmentGradientSystemReplace> final_model;
    Model_Of_Polygonal_Dislocation(std::vector <DisplacmentGradientSystemReplace> final_model);
    double uxxcalc(double x, double y, double z);
    double uxycalc(double x, double y, double z);
    double uxzcalc(double x, double y, double z);
    double uyxcalc(double x, double y, double z);
    double uyycalc(double x, double y, double z);
    double uyzcalc(double x, double y, double z);
    double uzxcalc(double x, double y, double z);
    double uzycalc(double x, double y, double z);
    double uzzcalc(double x, double y, double z);
};

#endif /* Distorsion_Source_hpp */
