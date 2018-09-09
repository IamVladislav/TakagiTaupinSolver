//
//  main.cpp
//  UnitTest
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include <iostream>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "Vector_And_Operations.hpp"
#include "System_Repclacement_Library.hpp"
#include "/Users/vladislavmetel/Desktop/TakagiTaupinSolver/TakagiTaupinSolver/System_Repclacement_Library.cpp"
#include "/Users/vladislavmetel/Desktop/TakagiTaupinSolver/TakagiTaupinSolver/Distorsion_Source.cpp"
#include "/Users/vladislavmetel/Desktop/TakagiTaupinSolver/TakagiTaupinSolver/Initialization_And_Calculation_Geometry.cpp"

TEST (Vector_Operation_Test, Equals) {
    Vector<double> first_vector;
    first_vector.c[0]=1;
    first_vector.c[1]=1;
    first_vector.c[2]=1;
    Vector<double> second_vector;
    second_vector.c[0]=1;
    second_vector.c[1]=1;
    second_vector.c[2]=1;
    ASSERT_TRUE(first_vector==second_vector);
}

TEST (Vector_Operation_Test, NotEquals) {
    Vector<int> first_vector;
    first_vector.c[0]=1;
    first_vector.c[1]=1;
    first_vector.c[2]=1;
    Vector<int> second_vector;
    second_vector.c[0]=1;
    second_vector.c[1]=0;
    second_vector.c[2]=1;
    ASSERT_TRUE(first_vector!=second_vector);
}

TEST (Vector_Operation_Test, Vector_Multiplication) {
    Vector<double> first_vector;
    Vector<double> second_vector;
    Vector<double> result_vector_expected;
    Vector<double> result_vector_actuall;
    //Case 1:
    first_vector.c[0]=1;
    first_vector.c[1]=1;
    first_vector.c[2]=1;
    second_vector.c[0]=-1;
    second_vector.c[1]=0;
    second_vector.c[2]=-1;
    result_vector_expected.c[0]=-1;
    result_vector_expected.c[1]=0;
    result_vector_expected.c[2]=1;
    result_vector_actuall=Vector_Multiplication<double, double>(first_vector, second_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
    //Case 2:
    first_vector.c[0]=1;
    first_vector.c[1]=0;
    first_vector.c[2]=0;
    second_vector.c[0]=0;
    second_vector.c[1]=1;
    second_vector.c[2]=0;
    result_vector_expected.c[0]=0;
    result_vector_expected.c[1]=0;
    result_vector_expected.c[2]=1;
    result_vector_actuall=Vector_Multiplication<double, double>(first_vector, second_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
}

TEST (Vector_Operation_Test, Vector_Absolute_Value) {
    Vector<int> first_vector;
    double result_actuall, result_expected;
    //Case 1:
    first_vector.c[0]=1;
    first_vector.c[1]=1;
    first_vector.c[2]=1;
    result_expected=sqrt(first_vector.c[0]*first_vector.c[0] + first_vector.c[1]*first_vector.c[1] + first_vector.c[2]*first_vector.c[2]);
    result_actuall=Vector_Absolute_Value(first_vector);
    ASSERT_EQ(result_actuall,result_expected);
    //Case 2:
    first_vector.c[0]=-1;
    first_vector.c[1]=2;
    first_vector.c[2]=1;
    result_expected=sqrt(first_vector.c[0]*first_vector.c[0] + first_vector.c[1]*first_vector.c[1] + first_vector.c[2]*first_vector.c[2]);
    result_actuall=Vector_Absolute_Value(first_vector);
    ASSERT_EQ(result_actuall,result_expected);
    result_actuall=Vector_Absolute_Value(first_vector);
    ASSERT_EQ(result_actuall,result_expected);
}

TEST (Vector_Operation_Test, Vector_Normalization) {
    Vector<int> first_vector;
    Vector<double> result_vector_expected;
    Vector<double> result_vector_actuall;
    //Case 1:
    first_vector.c[0]=1;
    first_vector.c[1]=1;
    first_vector.c[2]=1;
    result_vector_expected.c[0]=first_vector.c[0]/Vector_Absolute_Value(first_vector);
    result_vector_expected.c[1]=first_vector.c[1]/Vector_Absolute_Value(first_vector);
    result_vector_expected.c[2]=first_vector.c[2]/Vector_Absolute_Value(first_vector);
    result_vector_actuall=Vector_Normalization<int, double>(first_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
    //Case 2:
    first_vector.c[0]=1;
    first_vector.c[1]=0;
    first_vector.c[2]=0;
    result_vector_expected.c[0]=first_vector.c[0]/Vector_Absolute_Value(first_vector);
    result_vector_expected.c[1]=first_vector.c[1]/Vector_Absolute_Value(first_vector);
    result_vector_expected.c[2]=first_vector.c[2]/Vector_Absolute_Value(first_vector);
    result_vector_actuall=Vector_Normalization<int, double>(first_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
}

TEST (Vector_Operation_Test, Vector_Inverse) {
    Vector<double> first_vector;
    Vector<double> result_vector_expected;
    Vector<double> result_vector_actuall;
    //Case 1:
    first_vector.c[0]=1.0;
    first_vector.c[1]=-1.55;
    first_vector.c[2]=2.3;
    result_vector_expected.c[0]=-first_vector.c[0];
    result_vector_expected.c[1]=-first_vector.c[1];
    result_vector_expected.c[2]=-first_vector.c[2];
    result_vector_actuall=Vector_Inverse(first_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
    //Case 2:
    first_vector.c[0]=1;
    first_vector.c[1]=0;
    first_vector.c[2]=0;
    result_vector_expected.c[0]=-first_vector.c[0];
    result_vector_expected.c[1]=-first_vector.c[1];
    result_vector_expected.c[2]=-first_vector.c[2];
    result_vector_actuall=Vector_Inverse(first_vector);
    ASSERT_EQ(result_vector_actuall,result_vector_expected);
}

TEST (Vector_Operation_Test, Scalar_Multiplication) {
    Vector<double> first_vector;
    Vector<double> second_vector;
    double result_actuall, result_expected;
    //Case 1:
    first_vector.c[0]=1.0;
    first_vector.c[1]=-1.55;
    first_vector.c[2]=2.3;
    second_vector.c[0]=1;
    second_vector.c[1]=0;
    second_vector.c[2]=-2;
    result_expected=first_vector.c[0]*second_vector.c[0] + first_vector.c[1]*second_vector.c[1] + first_vector.c[2]*second_vector.c[2];
    result_actuall=Scalar_Multiplication(first_vector, second_vector);
    ASSERT_EQ(result_actuall,result_expected);
    //Case 2:
    first_vector.c[0]=1;
    first_vector.c[1]=0;
    first_vector.c[2]=0;
    second_vector.c[0]=0;
    second_vector.c[1]=1;
    second_vector.c[2]=0;
    result_expected=first_vector.c[0]*second_vector.c[0] + first_vector.c[1]*second_vector.c[1] + first_vector.c[2]*second_vector.c[2];
    result_actuall=Scalar_Multiplication(first_vector, second_vector);
    ASSERT_EQ(result_actuall,result_expected);
}

TEST (System_Replacement, Rotation_Matrix_Z){
    double phi,x,y,z,x_expected,y_expected,z_expected,xinv_expected,yinv_expected,zinv_expected;
    Rotation_Matrix_Z *test_object;
    //Case 1:
    phi = 0.43;
    x=2;
    y=3;
    z=2;
    x_expected=x*cos(phi)+y*sin(phi);
    y_expected=-x*sin(phi)+y*cos(phi);
    z_expected=z;
    xinv_expected=x*cos(phi)-y*sin(phi);
    yinv_expected=x*sin(phi)+y*cos(phi);
    zinv_expected=z;
    test_object = new Rotation_Matrix_Z(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 2:
    phi = 1.80;
    x=1;
    y=30;
    z=200;
    x_expected=x*cos(phi)+y*sin(phi);
    y_expected=-x*sin(phi)+y*cos(phi);
    z_expected=z;
    xinv_expected=x*cos(phi)-y*sin(phi);
    yinv_expected=x*sin(phi)+y*cos(phi);
    zinv_expected=z;
    test_object = new Rotation_Matrix_Z(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 3:
    phi = 3.20;
    x=-1;
    y=30;
    z=200;
    x_expected=x*cos(phi)+y*sin(phi);
    y_expected=-x*sin(phi)+y*cos(phi);
    z_expected=z;
    xinv_expected=x*cos(phi)-y*sin(phi);
    yinv_expected=x*sin(phi)+y*cos(phi);
    zinv_expected=z;
    test_object = new Rotation_Matrix_Z(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 4:
    phi = 5;
    x=-1;
    y=-30;
    z=-2;
    x_expected=x*cos(phi)+y*sin(phi);
    y_expected=-x*sin(phi)+y*cos(phi);
    z_expected=z;
    xinv_expected=x*cos(phi)-y*sin(phi);
    yinv_expected=x*sin(phi)+y*cos(phi);
    zinv_expected=z;
    test_object = new Rotation_Matrix_Z(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
}

TEST (System_Replacement, Rotation_Matrix_Y){
    double phi,x,y,z,x_expected,y_expected,z_expected,xinv_expected,yinv_expected,zinv_expected;
    Rotation_Matrix_Y *test_object;
    //Case 1:
    phi = 0.43;
    x=2;
    y=3;
    z=2;
    x_expected=x*cos(phi)+z*sin(phi);
    y_expected=y;
    z_expected=-x*sin(phi)+z*cos(phi);
    xinv_expected=x*cos(phi)-z*sin(phi);
    yinv_expected=y;
    zinv_expected=x*sin(phi)+z*cos(phi);
    test_object = new Rotation_Matrix_Y(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 2:
    phi = 1.80;
    x=1;
    y=30;
    z=200;
    x_expected=x*cos(phi)+z*sin(phi);
    y_expected=y;
    z_expected=-x*sin(phi)+z*cos(phi);
    xinv_expected=x*cos(phi)-z*sin(phi);
    yinv_expected=y;
    zinv_expected=x*sin(phi)+z*cos(phi);
    test_object = new Rotation_Matrix_Y(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 3:
    phi = 3.20;
    x=-1;
    y=30;
    z=200;
    x_expected=x*cos(phi)+z*sin(phi);
    y_expected=y;
    z_expected=-x*sin(phi)+z*cos(phi);
    xinv_expected=x*cos(phi)-z*sin(phi);
    yinv_expected=y;
    zinv_expected=x*sin(phi)+z*cos(phi);
    test_object = new Rotation_Matrix_Y(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 4:
    phi = 5;
    x=-1;
    y=-30;
    z=-2;
    x_expected=x*cos(phi)+z*sin(phi);
    y_expected=y;
    z_expected=-x*sin(phi)+z*cos(phi);
    xinv_expected=x*cos(phi)-z*sin(phi);
    yinv_expected=y;
    zinv_expected=x*sin(phi)+z*cos(phi);
    test_object = new Rotation_Matrix_Y(phi);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
}

TEST (System_Replacement, Parallel_Shift){
    double x_shift,y_shift,z_shift,x,y,z,x_expected,y_expected,z_expected,xinv_expected,yinv_expected,zinv_expected;
    Parallel_Shfit *test_object;
    //Case 1:
    x_shift=10;
    y_shift=2;
    z_shift=3;
    x=2;
    y=3;
    z=2;
    x_expected=x-x_shift;
    y_expected=y-y_shift;
    z_expected=z-z_shift;
    xinv_expected=x+x_shift;;
    yinv_expected=y+y_shift;;
    zinv_expected=z+z_shift;;
    test_object = new Parallel_Shfit(x_shift,y_shift,z_shift);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
    //Case 2:
    x_shift=-10;
    y_shift=0;
    z_shift=0;
    x=-1000;
    y=8;
    z=124;
    x_expected=x-x_shift;
    y_expected=y-y_shift;
    z_expected=z-z_shift;
    xinv_expected=x+x_shift;;
    yinv_expected=y+y_shift;;
    zinv_expected=z+z_shift;;
    test_object = new Parallel_Shfit(x_shift,y_shift,z_shift);
    test_object->Basis(x, y, z);
    test_object->Basis_Inverse(x, y, z);
    ASSERT_NEAR(test_object->xc,x_expected,0.01);
    ASSERT_NEAR(test_object->yc,y_expected,0.01);
    ASSERT_NEAR(test_object->zc,z_expected,0.01);
    ASSERT_NEAR(test_object->xcinv,xinv_expected,0.01);
    ASSERT_NEAR(test_object->ycinv,yinv_expected,0.01);
    ASSERT_NEAR(test_object->zcinv,zinv_expected,0.01);
}

TEST (System_Replacement, Rotation_Matrix_Z_Vector){
    double phi;
    Vector<double> vector_expected, vectorinv_expected,vector_input;
    Rotation_Matrix_Z_Vector *test_object;
    //Case 1:
    phi = 0.43;
    vector_input.c[0]=3;
    vector_input.c[0]=-3;
    vector_input.c[0]=1;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[1]*sin(phi);
    vector_expected.c[1]=-vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vector_expected.c[2]=vector_input.c[2];
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[1]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vectorinv_expected.c[2]=vector_input.c[2];
    test_object = new Rotation_Matrix_Z_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 2:
    phi = 1.80;
    vector_input.c[0]=0;
    vector_input.c[0]=-3;
    vector_input.c[0]=1;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[1]*sin(phi);
    vector_expected.c[1]=-vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vector_expected.c[2]=vector_input.c[2];
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[1]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vectorinv_expected.c[2]=vector_input.c[2];
    test_object = new Rotation_Matrix_Z_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 3:
    phi = 3.20;
    vector_input.c[0]=3;
    vector_input.c[0]=-3;
    vector_input.c[0]=45;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[1]*sin(phi);
    vector_expected.c[1]=-vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vector_expected.c[2]=vector_input.c[2];
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[1]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vectorinv_expected.c[2]=vector_input.c[2];
    test_object = new Rotation_Matrix_Z_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 4:
    phi = 5;
    vector_input.c[0]=0;
    vector_input.c[0]=0;
    vector_input.c[0]=13;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[1]*sin(phi);
    vector_expected.c[1]=-vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vector_expected.c[2]=vector_input.c[2];
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[1]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[0]*sin(phi)+vector_input.c[1]*cos(phi);
    vectorinv_expected.c[2]=vector_input.c[2];
    test_object = new Rotation_Matrix_Z_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
}

TEST (System_Replacement, Rotation_Matrix_Y_Vector){
    double phi;
    Vector<double> vector_expected, vectorinv_expected,vector_input;
    Rotation_Matrix_Y_Vector *test_object;
    //Case 1:
    phi = 0.43;
    vector_input.c[0]=3;
    vector_input.c[0]=-3;
    vector_input.c[0]=1;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[2]*sin(phi);
    vector_expected.c[1]=vector_input.c[1];
    vector_expected.c[2]=-vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[2]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[1];
    vectorinv_expected.c[2]=vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    test_object = new Rotation_Matrix_Y_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 2:
    phi = 1.80;
    vector_input.c[0]=0;
    vector_input.c[0]=-3;
    vector_input.c[0]=1;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[2]*sin(phi);
    vector_expected.c[1]=vector_input.c[1];
    vector_expected.c[2]=-vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[2]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[1];
    vectorinv_expected.c[2]=vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    test_object = new Rotation_Matrix_Y_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 3:
    phi = 3.20;
    vector_input.c[0]=3;
    vector_input.c[0]=-3;
    vector_input.c[0]=45;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[2]*sin(phi);
    vector_expected.c[1]=vector_input.c[1];
    vector_expected.c[2]=-vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[2]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[1];
    vectorinv_expected.c[2]=vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    test_object = new Rotation_Matrix_Y_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
    //Case 4:
    phi = 5;
    vector_input.c[0]=0;
    vector_input.c[0]=0;
    vector_input.c[0]=13;
    vector_expected.c[0]=vector_input.c[0]*cos(phi)+vector_input.c[2]*sin(phi);
    vector_expected.c[1]=vector_input.c[1];
    vector_expected.c[2]=-vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    vectorinv_expected.c[0]=vector_input.c[0]*cos(phi)-vector_input.c[2]*sin(phi);
    vectorinv_expected.c[1]=vector_input.c[1];
    vectorinv_expected.c[2]=vector_input.c[0]*sin(phi)+vector_input.c[2]*cos(phi);
    test_object = new Rotation_Matrix_Y_Vector(phi);
    test_object->Basis(vector_input);
    test_object->Basis_Inverse(vector_input);
    ASSERT_NEAR(test_object->vector.c[0],vector_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector.c[1],vector_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector.c[2],vector_expected.c[2],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[0],vectorinv_expected.c[0],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[1],vectorinv_expected.c[1],0.01);
    ASSERT_NEAR(test_object->vector_inv.c[2],vectorinv_expected.c[2],0.01);
}

TEST (Distrotion_Field, Angular_Dislocation){
    Vector<double> burgers_vector, exit_point_coordinate;
    double x, y, z, angle, depth, nu, segment_lenght, uxx_expected, uxy_expected, uxz_expected, uyx_expected,uyy_expected,uyz_expected,uzx_expected,uzy_expected,uzz_expected;
    AngularDislocation *test_obj;
    double precision = 0.000001;
    //Case 1:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=3;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=INF;
    exit_point_coordinate.c[1]=INF;
    exit_point_coordinate.c[2]=-100;
    x=5;
    y=6;
    z=8;
    nu = 0.4;
    segment_lenght = 0;
    angle = PI_2-0.01*PI_2;
    depth = 100;
    uxx_expected = 0.0633306;
    uxy_expected = -0.0689766;
    uxz_expected = 0.00845708;
    uyx_expected = 0.0823332;
    uyy_expected = -0.0984343;
    uyz_expected = 0.0308324;
    uzx_expected = 0.0266809;
    uzy_expected = -0.0725122;
    uzz_expected = 0.0348514;
    test_obj = new AngularDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 2:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=-1;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=INF;
    exit_point_coordinate.c[1]=INF;
    exit_point_coordinate.c[2]=-300;
    x=124;
    y=-23;
    z=400;
    nu = 0.3;
    segment_lenght = 0;
    angle = PI_2/2;
    depth = 300;
    uxx_expected = -0.000532479;
    uxy_expected = -0.0085049;
    uxz_expected = -0.000302145;
    uyx_expected = 0.00127729;
    uyy_expected = 0.001341;
    uyz_expected = -0.000289173;
    uzx_expected = -0.000609022;
    uzy_expected = -0.0051017;
    uzz_expected = -0.0000574385;
    test_obj = new AngularDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 3:
    burgers_vector.c[0]=0;
    burgers_vector.c[1]=0;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=INF;
    exit_point_coordinate.c[1]=INF;
    exit_point_coordinate.c[2]=-150;
    x=0;
    y=5;
    z=300;
    nu = 0.2;
    segment_lenght = 100;
    angle = PI_2/3;
    depth = 150;
    uxx_expected = 0.000011972;
    uxy_expected = 0.000591855;
    uxz_expected = -5.809652113384596e-6;
    uyx_expected = -0.000390611;
    uyy_expected = -0.0000142997;
    uyz_expected = -0.000111124;
    uzx_expected = 0.000197452;
    uzy_expected = 0.00333298;
    uzz_expected = 0.0000133638;
    test_obj = new AngularDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
}

TEST (Distrotion_Field, Beam_Dislocation){
    Vector<double> burgers_vector, exit_point_coordinate;
    double x, y, z, angle, depth, nu, segment_lenght, uxx_expected, uxy_expected, uxz_expected, uyx_expected,uyy_expected,uyz_expected,uzx_expected,uzy_expected,uzz_expected;
    BeamDislocation *test_obj;
    double precision = 0.000001;
    //Case 1:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=3;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=0;
    exit_point_coordinate.c[1]=0;
    exit_point_coordinate.c[2]=-100;
    x=5;
    y=6;
    z=8;
    nu = 0.4;
    segment_lenght = 0;
    angle = 0.001;
    depth = 100;
    uxx_expected = -0.0722531;
    uxy_expected = 0.0700228;
    uxz_expected = -0.000103008;
    uyx_expected = -0.099078;
    uyy_expected = 0.0692551;
    uyz_expected = 0.000297339;
    uzx_expected = -0.0445411;
    uzy_expected = 0.0357161;
    uzz_expected = 0.0000774161;
    test_obj = new BeamDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 2:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=-1;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=-300;
    exit_point_coordinate.c[1]=0;
    exit_point_coordinate.c[2]=-300;
    x=124;
    y=-23;
    z=400;
    nu = 0.3;
    segment_lenght = 0;
    angle = PI_2/2;
    depth = 300;
    uxx_expected = 0.000304388;
    uxy_expected = -0.00232471;
    uxz_expected = -0.000287393;
    uyx_expected = 0.000272403;
    uyy_expected = 0.000177961;
    uyz_expected = -0.000272164;
    uzx_expected = -0.0000110171;
    uzy_expected = -0.0018624;
    uzz_expected = -0.0000278826;
    test_obj = new BeamDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 3:
    burgers_vector.c[0]=0;
    burgers_vector.c[1]=0;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=86.60;
    exit_point_coordinate.c[1]=0;
    exit_point_coordinate.c[2]=-150;
    x=0;
    y=5;
    z=300;
    nu = 0.2;
    segment_lenght = 100;
    angle = PI_2/3;
    depth = 150;
    uxx_expected = 0.000138294;
    uxy_expected = 0.00158921;
    uxz_expected = -0.0000750038;
    uyx_expected = -0.000798065;
    uyy_expected = -0.000253243;
    uyz_expected = 0.000333745;
    uzx_expected = -0.000591728;
    uzy_expected = -0.00778931;
    uzz_expected = 0.000331918;
    test_obj = new BeamDislocation(nu, angle, depth, exit_point_coordinate, segment_lenght, burgers_vector);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
}

TEST (Model_Rotation, Angular_Dislocation){
    Vector<double> burgers_vector, exit_point_coordinate;
    double x, y, z, depth, nu, segment_lenght, uxx_expected, uxy_expected, uxz_expected, uyx_expected,uyy_expected,uyz_expected,uzx_expected,uzy_expected,uzz_expected,phi,kappa;
    DisplacmentGradientSystemReplace *test_obj;
    double precision = 0.000001;
    //Case 1:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=3;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=-50;
    exit_point_coordinate.c[1]=-100;
    exit_point_coordinate.c[2]=-100;
    x=5;
    y=6;
    z=8;
    nu = 0.4;
    segment_lenght = 0;
    depth = 100;
    phi = 1.10714871779409;
    kappa = 0.8410686705679302;
    uxx_expected = 0.21984977;
    uxy_expected = 0.01426786;
    uxz_expected = -0.15207801;
    uyx_expected = 0.19006567;
    uyy_expected = 0.00429901;
    uyz_expected = -0.11927170;
    uzx_expected = 0.22915184;
    uzy_expected = 0.07138488;
    uzz_expected = -0.19552116 ;
    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 1);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 2:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=-1;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=-150;
    exit_point_coordinate.c[1]=-173.205;
    exit_point_coordinate.c[2]=-300;
    x=124;
    y=-23;
    z=400;
    nu = 0.3;
    segment_lenght = 0;
    depth = 300;
    phi = 0.8570719478501309;
    kappa = 0.652251026319475;
    uxx_expected = -0.00001399;
    uxy_expected = -0.00628310;
    uxz_expected = -0.00034305;
    uyx_expected = 0.00006900;
    uyy_expected = 0.00175178;
    uyz_expected = 0.00012733;
    uzx_expected = 0.00143529;
    uzy_expected = -0.00431896;
    uzz_expected = -0.00058953;
    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 1);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //    //Case 3:
    //    burgers_vector.c[0]=0;
    //    burgers_vector.c[1]=0;
    //    burgers_vector.c[2]=3;
    //    exit_point_coordinate.c[0]=-600;
    //    exit_point_coordinate.c[1]=-259.8076211353315;
    //    exit_point_coordinate.c[2]=-150;
    //    x=0;
    //    y=5;
    //    z=300;
    //    nu = 0.2;
    //    segment_lenght = 100;
    //    depth = 150;
    //    phi = 0.4086378550975923;
    //    kappa = 1.345282920896765;
    //    uxx_expected = -0.00005336;
    //    uxy_expected = 0.00064878;
    //    uxz_expected = 0.00003711;
    //    uyx_expected = -0.00053581;
    //    uyy_expected = 0.00005031;
    //    uyz_expected = -0.00011076;
    //    uzx_expected = -0.00155506;
    //    uzy_expected = 0.00424228;
    //    uzz_expected = 0.00000749;
    //    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 1);
    //    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    //    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    //    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    //    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    //    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    //    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    //    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    //    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    //    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
}

TEST (Model_Rotation, Beam_Dislocation){
    Vector<double> burgers_vector, exit_point_coordinate;
    double x, y, z, depth, nu, segment_lenght, uxx_expected, uxy_expected, uxz_expected, uyx_expected,uyy_expected,uyz_expected,uzx_expected,uzy_expected,uzz_expected,phi,kappa;
    DisplacmentGradientSystemReplace *test_obj;
    double precision = 0.000001;
    //Case 1:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=3;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=-50;
    exit_point_coordinate.c[1]=-100;
    exit_point_coordinate.c[2]=-100;
    x=5;
    y=6;
    z=8;
    nu = 0.4;
    segment_lenght = 0;
    depth = 100;
    phi = 1.10714871779409;
    kappa = 0.8410686705679302;
    uxx_expected = 0.157768;
    uxy_expected = 0.0786772;
    uxz_expected = -0.157188;
    uyx_expected = 0.106111;
    uyy_expected = 0.0668987;
    uyz_expected = -0.120004;
    uzx_expected = 0.208833;
    uzy_expected = 0.093635;
    uzz_expected = -0.198133;
    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 2);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //Case 2:
    burgers_vector.c[0]=3;
    burgers_vector.c[1]=-1;
    burgers_vector.c[2]=3;
    exit_point_coordinate.c[0]=-150;
    exit_point_coordinate.c[1]=-173.205;
    exit_point_coordinate.c[2]=-300;
    x=124;
    y=-23;
    z=400;
    nu = 0.3;
    segment_lenght = 0;
    depth = 300;
    phi = 0.8570719478501309;
    kappa = 0.652251026319475;
    uxx_expected = 0.000851692;
    uxy_expected = -7.47757e-05;
    uxz_expected = -0.000341787;
    uyx_expected = -0.0009311;
    uyy_expected = 0.000563886;
    uyz_expected = 0.000140909;
    uzx_expected = 0.00213133;
    uzy_expected = -0.00100848;
    uzz_expected = -0.00054643;
    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 2);
    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
    //    //Case 3:
    //    burgers_vector.c[0]=0;
    //    burgers_vector.c[1]=0;
    //    burgers_vector.c[2]=3;
    //    exit_point_coordinate.c[0]=-600;
    //    exit_point_coordinate.c[1]=-259.8076211353315;
    //    exit_point_coordinate.c[2]=-150;
    //    x=0;
    //    y=5;
    //    z=300;
    //    nu = 0.2;
    //    segment_lenght = 100;
    //    depth = 150;
    //    phi = 0.4086378550975923;
    //    kappa = 1.345282920896765;
    //    uxx_expected = -0.00005336;
    //    uxy_expected = 0.00064878;
    //    uxz_expected = 0.00003711;
    //    uyx_expected = -0.00053581;
    //    uyy_expected = 0.00005031;
    //    uyz_expected = -0.00011076;
    //    uzx_expected = -0.00155506;
    //    uzy_expected = 0.00424228;
    //    uzz_expected = 0.00000749;
    //    test_obj = new DisplacmentGradientSystemReplace(nu, depth, phi, kappa, exit_point_coordinate, segment_lenght, burgers_vector, 1);
    //    ASSERT_NEAR(test_obj->uxxcalc(x, y, z), uxx_expected, precision);
    //    ASSERT_NEAR(test_obj->uxycalc(x, y, z), uxy_expected, precision);
    //    ASSERT_NEAR(test_obj->uxzcalc(x, y, z), uxz_expected, precision);
    //    ASSERT_NEAR(test_obj->uyxcalc(x, y, z), uyx_expected, precision);
    //    ASSERT_NEAR(test_obj->uyycalc(x, y, z), uyy_expected, precision);
    //    ASSERT_NEAR(test_obj->uyzcalc(x, y, z), uyz_expected, precision);
    //    ASSERT_NEAR(test_obj->uzxcalc(x, y, z), uzx_expected, precision);
    //    ASSERT_NEAR(test_obj->uzycalc(x, y, z), uzy_expected, precision);
    //    ASSERT_NEAR(test_obj->uzzcalc(x, y, z), uzz_expected, precision);
}

TEST (Initilalization_And_Calculation_Geometry, Glide_Plane_Calculation){
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
    //Nu
    double nu=0.4;
    InitialisationGeometry *test_obj = new InitialisationGeometry(diffraction_vector, normal_vector, tay_vector, segment_lenght, dislocation_depth, nu, burgers_vector, number_of_dislocation);
    //Someday i replace upper code on object mock, becouse this initialization really don't matter, but i need the object for testing.
    //Case 1:
    std::vector <Vector<double>> tay_vector_test_1;
    tay_1.c[0]=1;
    tay_1.c[1]=0;
    tay_1.c[2]=0;
    tay_vector_test_1.push_back(tay_1);
    tay_2.c[0]=0;
    tay_2.c[1]=1;
    tay_2.c[2]=0;
    tay_vector_test_1.push_back(tay_2);
    test_obj->Glide_Plane_Calculate(tay_vector_test_1);
    ASSERT_EQ(test_obj->glide_plane_vector.c[0], 0);
    ASSERT_EQ(test_obj->glide_plane_vector.c[1], 0);
    ASSERT_EQ(test_obj->glide_plane_vector.c[2], 1);
    //Case 2:
    std::vector <Vector<double>> tay_vector_test_2;
    tay_1.c[0]=1;
    tay_1.c[1]=1;
    tay_1.c[2]=1;
    tay_vector_test_2.push_back(tay_1);
    tay_2.c[0]=-1;
    tay_2.c[1]=-1;
    tay_2.c[2]=-1;
    tay_vector_test_2.push_back(tay_2);
    test_obj->Glide_Plane_Calculate(tay_vector_test_2);
    ASSERT_EQ(test_obj->glide_plane_vector.c[0], 0);
    ASSERT_EQ(test_obj->glide_plane_vector.c[1], 0);
    ASSERT_EQ(test_obj->glide_plane_vector.c[2], 0);
    //Case 3:(Fake input data)
    std::vector <Vector<double>> tay_vector_test_3;
    tay_1.c[0]=0;
    tay_1.c[1]=-1;
    tay_1.c[2]=1;
    tay_vector_test_3.push_back(tay_1);
    tay_2.c[0]=-1;
    tay_2.c[1]=-1;
    tay_2.c[2]=0;
    tay_vector_test_3.push_back(tay_2);
    Vector<double> tay_3;
    tay_3.c[0]=1;
    tay_3.c[1]=0;
    tay_3.c[2]=1;
    tay_vector_test_3.push_back(tay_3);
    test_obj->Glide_Plane_Calculate(tay_vector_test_3);
    ASSERT_EQ(test_obj->glide_plane_vector.c[0], 1);
    ASSERT_EQ(test_obj->glide_plane_vector.c[1], -1);
    ASSERT_EQ(test_obj->glide_plane_vector.c[2], -1);
    //Case 4:
    std::vector <Vector<double>> tay_vector_test_4;
    tay_1.c[0]=1;
    tay_1.c[1]=0;
    tay_1.c[2]=0;
    tay_vector_test_4.push_back(tay_1);
    tay_2.c[0]=1;
    tay_2.c[1]=1;
    tay_2.c[2]=1;
    tay_vector_test_4.push_back(tay_2);
    tay_3.c[0]=10;
    tay_3.c[1]=1;
    tay_3.c[2]=1;
    tay_vector_test_4.push_back(tay_3);
    try{
        test_obj->Glide_Plane_Calculate(tay_vector_test_4);
    }
    catch (const char *str){
        SUCCEED();
    }
}

int main(int argc, char * argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
