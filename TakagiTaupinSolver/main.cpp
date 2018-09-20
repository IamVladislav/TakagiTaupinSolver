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
#include <fstream>

int main(int argc, const char * argv[]) {
    //Diffraction vector:
    Vector<double> diffraction_vector;
    diffraction_vector.c[0]=2;
    diffraction_vector.c[1]=-2;
    diffraction_vector.c[2]=0;
    //Normal to surface vector:
    Vector<double> normal_vector;
    normal_vector.c[0]=1;
    normal_vector.c[1]=1;
    normal_vector.c[2]=1;
    //Vector of tay-vectors:
    std::vector <Vector<double>> tay_vector;
    Vector<double> tay_1;
    tay_1.c[0]=-1;
    tay_1.c[1]=-1;
    tay_1.c[2]=0;
    tay_vector.push_back(tay_1);
    Vector<double> tay_2;
    tay_2.c[0]=-1;
    tay_2.c[1]=-1;
    tay_2.c[2]=0;
    tay_vector.push_back(tay_2);
    Vector<double> tay_3;
    tay_3.c[0]=0;
    tay_3.c[1]=-1;
    tay_3.c[2]=1;
    tay_vector.push_back(tay_3);
    Vector<double> tay_4;
    tay_4.c[0]=0;
    tay_4.c[1]=-1;
    tay_4.c[2]=1;
    tay_vector.push_back(tay_4);
    Vector<double> tay_5;
    tay_5.c[0]=1;
    tay_5.c[1]=0;
    tay_5.c[2]=1;
    tay_vector.push_back(tay_5);
    Vector<double> tay_6;
    tay_6.c[0]=1;
    tay_6.c[1]=0;
    tay_6.c[2]=1;
    tay_vector.push_back(tay_6);
    //Burger's vectorof dislocation:
    Vector<double> burgers_vector;
    burgers_vector.c[0]=0;
    burgers_vector.c[1]=-1;
    burgers_vector.c[2]=1;
    //Number of dislocation in one place:
    int number_of_dislocation=1;
    //Dislocation depth:
    std::vector <double> dislocation_depth;
    dislocation_depth.push_back(15);
    dislocation_depth.push_back(15);
    dislocation_depth.push_back(15);
    dislocation_depth.push_back(15);
    dislocation_depth.push_back(15);
    dislocation_depth.push_back(15);
    //Segment lengt
    std::vector<double> segment_lenght;
    segment_lenght.push_back(0);
    segment_lenght.push_back(0);
    segment_lenght.push_back(0);
    segment_lenght.push_back(-40);
    segment_lenght.push_back(-40);
    segment_lenght.push_back(-40);
    //Nu
    double nu=0.4;
    //Type_vector
    std::vector<int> type_list;
    type_list.push_back(1);
    type_list.push_back(2);
    type_list.push_back(1);
    type_list.push_back(1);
    type_list.push_back(1);
    type_list.push_back(2);
    InitialisationGeometry *init_obj = new InitialisationGeometry(diffraction_vector, normal_vector, tay_vector, segment_lenght, type_list, dislocation_depth, nu, burgers_vector, number_of_dislocation);
    Model_Of_Polygonal_Dislocation *model_obj;
    model_obj=init_obj->ModelOutput();
    
    //Calculation section
    int i_max,k_max;
    double x_min, x_max, y_min, y_max, x, y, z, x_step, y_step;
    i_max=1000;
    k_max=1000;
    x_min=-20.01;
    x_max=59.99;
    y_min=-30.01;
    y_max=59.99;
    x_step=(x_max-x_min)/i_max;
    y_step=(y_max-y_min)/k_max;
    x=x_min;
    y=y_min;
    std::ofstream fout;
    fout.open("Distribution_of_uxz");
    for(int i=0;i<i_max;i++){
        x=x+x_step;
        for(int k=0;k<k_max;k++){
            y=y+y_step;
            fout << model_obj->uyxcalc(x, y, 0) << std::endl;
        }
        y=y_min;
    }
    fout.close();
    return 0;
}
