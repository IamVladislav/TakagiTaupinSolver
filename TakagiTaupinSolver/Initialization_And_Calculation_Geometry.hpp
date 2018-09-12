//
//  Initialization_And_Calculation_Geometry.hpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 08/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#ifndef Initialization_And_Calculation_Geometry_hpp
#define Initialization_And_Calculation_Geometry_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
#include "Vector_And_Operations.hpp"
#include "System_Repclacement_Library.hpp"
#include "Distorsion_Source.hpp"

class InitialisationGeometry{
public:
    Vector<double> x_vector,y_vector,z_vector, glide_plane_vector, b_vector;
    std::vector <double> kappa, phi, dislocation_depth;
    std::vector <Vector<double>> exit_point;
    std::vector <DisplacmentGradientSystemReplace> final_model;
    double nu;
    InitialisationGeometry(Vector<double> diffraction_vector,Vector<double> normal_vector,std::vector <Vector<double>> direction_vector_list,std::vector <double> segment_lenght, std::vector <int> type_list, std::vector <double> dislocation_depth, double nu, Vector<double> burgers_vector, int nubmer_of_dislocation);
    void Model_Creation(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenghtm, std::vector <int> type_list);
    Model_Of_Polygonal_Dislocation* ModelOutput();
    void Glide_Plane(std::vector <Vector<double>> direction_vector_list);
    void Glide_Plane_Calculate(std::vector <Vector<double>> direction_vector_list);
    void Exit_Point_Coordinate(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght);
    void Burgers_Vector_Into_System_Coordinates(Vector<double> b_initial);
    void Angle_Between_Dislocation_Axis_And_Calculation_Axis(std::vector <Vector<double>> direction_vector_list);
    template <typename type> bool is_Right_Normal(Vector<type> first_vector, Vector<type> second_vector);
};

#endif /* Initialization_And_Calculation_Geometry_hpp */
