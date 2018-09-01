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

int main(int argc, char * argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
