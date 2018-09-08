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
#include "System_Repclacement_Library.hpp"
#include "Distorsion_Source.hpp"
#include "Initialization_And_Calculation_Geometry.hpp"

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
    //Nu
    double nu=0.4;
    //Start of initialization
//    InitialisationGeometry *init_obj = new InitialisationGeometry(diffraction_vector, normal_vector,tay_vector, burgers_vector, number_of_dislocation,dislocation_depth,segment_lenght);
//    AngularDislocation *test = new AngularDislocation();
//    test->burgers_vector=burgers_vector;
//    test->depth=dislocation_depth;
//    test->kappa=1.0;
//    test->nu=0.4;
//    std::cout<<test->uxxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uxycalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uxzcalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uyxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uyycalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uyzcalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uzxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uzycalc(5, 5, 5)<<std::endl;
//    std::cout<<test->uzzcalc(5, 5, 5)<<std::endl;
    //
//    BeamDislocation *test1 = new BeamDislocation();
//    test1->burgers_vector=burgers_vector;
//    test1->depth=dislocation_depth;
//    test1->theta=1.0;
//    test1->nu=0.4;
//    std::cout<<test1->uxxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uxycalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uxzcalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uyxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uyycalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uyzcalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uzxcalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uzycalc(5, 5, 5)<<std::endl;
//    std::cout<<test1->uzzcalc(5, 5, 5)<<std::endl;
//    Vector<double> exit_point_coordinates;
//    exit_point_coordinates.c[0]=-1;
//    exit_point_coordinates.c[1]=-1;
//    exit_point_coordinates.c[2]=-1;
//    DisplacmentGradientSystemReplace *test2 = new DisplacmentGradientSystemReplace(0.4, dislocation_depth, 1.0, 1.0, exit_point_coordinates,burgers_vector, 2);
    InitialisationGeometry *init_obj = new InitialisationGeometry(diffraction_vector, normal_vector, tay_vector, segment_lenght, dislocation_depth, nu, burgers_vector, number_of_dislocation);
    //Test
    Model_Of_Polygonal_Dislocation *test;
    test=init_obj->ModelOutput();
    std::cout<<test->uxxcalc(10, 10, 10)<<std::endl;
    return 0;
}
