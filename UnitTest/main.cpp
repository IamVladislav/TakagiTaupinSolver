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

int main(int argc, char * argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
