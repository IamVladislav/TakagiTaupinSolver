//
//  Initialization_And_Calculation_Geometry.cpp
//  TakagiTaupinSolver
//
//  Created by Vladislav Metel on 08/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include "Initialization_And_Calculation_Geometry.hpp"

InitialisationGeometry::InitialisationGeometry(Vector<double> diffraction_vector,Vector<double> normal_vector,std::vector <Vector<double>> direction_vector_list,std::vector <double> segment_lenght, std::vector <int> type_list, std::vector <double> dislocation_depth, double nu, Vector<double> burgers_vector, int nubmer_of_dislocation){
    this->nu=nu;
    this->dislocation_depth=dislocation_depth;
    x_vector=Vector_Normalization<double, double>(diffraction_vector);
    z_vector=Vector_Normalization<double, double>(Vector_Inverse(normal_vector));
    y_vector=Vector_Multiplication<double, double>(z_vector, x_vector);
    try{
        Glide_Plane(direction_vector_list);
    }
    catch (const char *str){
        std::cout << str;
    }
    burgers_vector.c[0]=burgers_vector.c[0]*nubmer_of_dislocation;
    burgers_vector.c[1]=burgers_vector.c[1]*nubmer_of_dislocation;
    burgers_vector.c[2]=burgers_vector.c[2]*nubmer_of_dislocation;
    Burgers_Vector_Into_System_Coordinates(burgers_vector);
    this->dislocation_depth=dislocation_depth;
    Angle_Between_Dislocation_Axis_And_Calculation_Axis(direction_vector_list);
    Exit_Point_Coordinate(direction_vector_list, segment_lenght);
    Segment_Lenght(segment_lenght);
    Model_Creation(direction_vector_list, type_list);
}
void InitialisationGeometry::Model_Creation(std::vector <Vector<double>> direction_vector_list, std::vector <int> type_list){
    DisplacmentGradientSystemReplace object_1(this->nu, this->dislocation_depth[0], this->phi[0], this->kappa[0], this->exit_point[0], segment_lenght[0], Vector_Inverse(this->b_vector), type_list[0]);
        final_model.push_back(object_1);
    DisplacmentGradientSystemReplace object_2(this->nu, this->dislocation_depth[1], this->phi[1], this->kappa[1], this->exit_point[1], segment_lenght[1], this->b_vector, type_list[1]);
    final_model.push_back(object_2);
    DisplacmentGradientSystemReplace object_3(this->nu, this->dislocation_depth[2], this->phi[2], this->kappa[2], this->exit_point[2], segment_lenght[2], this->b_vector, type_list[2]);
    final_model.push_back(object_3);
}
Model_Of_Polygonal_Dislocation* InitialisationGeometry::ModelOutput(){//mock?
    Model_Of_Polygonal_Dislocation* model = new Model_Of_Polygonal_Dislocation(final_model);
    return model;
}
void InitialisationGeometry::Glide_Plane(std::vector <Vector<double>> direction_vector_list){
    switch (direction_vector_list.size()) {
        case 0:
            std::cerr<<"Can't calculate a glide plane, list is empty!"<<std::endl;
            break;
        case 1:
            std::cerr<<"Can't calculate a glide plane, list has only one vector!"<<std::endl;
            break;
        default:
            Glide_Plane_Calculate(direction_vector_list);
            if(!is_Right_Normal(glide_plane_vector, Vector_Inverse(z_vector))){
                glide_plane_vector=Vector_Inverse(glide_plane_vector);
            }
            break;
    }
}
void InitialisationGeometry::Glide_Plane_Calculate(std::vector <Vector<double>> direction_vector_list){
    Vector<double> comparison_vector;
    Vector<double> zero_vector;
    zero_vector.c[0]=0;
    zero_vector.c[1]=0;
    zero_vector.c[2]=0;
    for(int i=1;i<direction_vector_list.size();i++){
        glide_plane_vector=Vector_Multiplication<double, double>(direction_vector_list[i-1], direction_vector_list[i]);
        if (i!=1 and (comparison_vector!=glide_plane_vector && Vector_Inverse(comparison_vector)!=glide_plane_vector) && glide_plane_vector==zero_vector){
            throw "This vectors haven't overall glide plane";
        }
        comparison_vector=glide_plane_vector;
    }
}
void InitialisationGeometry::Exit_Point_Coordinate(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght){
    Vector<double> exit_point_calculation;
    for(int i=0;i<direction_vector_list.size();i++)
    {
        exit_point_calculation.c[2]=-dislocation_depth[i];
        if(kappa[i]>=(-0.01*PI_2) && kappa[i]<=(+0.01*PI_2)){
            exit_point_calculation.c[0]=0;
            exit_point_calculation.c[1]=0;
        }
        else{
            if(kappa[i]>=(PI_2-0.005*PI_2) && kappa[i]<=(PI_2+0.005*PI_2)){
                exit_point_calculation.c[0]=INF;
                exit_point_calculation.c[1]=INF;
            }
            else{
                exit_point_calculation.c[0]=-dislocation_depth[i]*tan(kappa[i])*cos(phi[i]);
                exit_point_calculation.c[1]=-dislocation_depth[i]*tan(kappa[i])*sin(phi[i]);
            }
        }
        exit_point.push_back(exit_point_calculation);
    }
}
void InitialisationGeometry::Burgers_Vector_Into_System_Coordinates(Vector<double> b_initial){
    b_vector.c[0]=-((-(b_initial.c[2]*y_vector.c[1]*z_vector.c[0]) + b_initial.c[1]*y_vector.c[2]*z_vector.c[0] + b_initial.c[2]*y_vector.c[0]*z_vector.c[1] - b_initial.c[0]*y_vector.c[2]*z_vector.c[1] - b_initial.c[1]*y_vector.c[0]*z_vector.c[2] + b_initial.c[0]*y_vector.c[1]*z_vector.c[2])/ (x_vector.c[2]*y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
    b_vector.c[1]=-((b_initial.c[2]*x_vector.c[1]*z_vector.c[0] - b_initial.c[1]*x_vector.c[2]*z_vector.c[0] - b_initial.c[2]*x_vector.c[0]*z_vector.c[1] + b_initial.c[0]*x_vector.c[2]*z_vector.c[1] + b_initial.c[1]*x_vector.c[0]*z_vector.c[2] - b_initial.c[0]*x_vector.c[1]*z_vector.c[2])/(x_vector.c[2]* y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
    b_vector.c[2]=-((b_initial.c[2]*x_vector.c[1]*y_vector.c[0] - b_initial.c[1]*x_vector.c[2]*y_vector.c[0] - b_initial.c[2]*x_vector.c[0]*y_vector.c[1] + b_initial.c[0]*x_vector.c[2]*y_vector.c[1] + b_initial.c[1]*x_vector.c[0]*y_vector.c[2] - b_initial.c[0]*x_vector.c[1]* y_vector.c[2])/(-(x_vector.c[2]*y_vector.c[1]*z_vector.c[0]) + x_vector.c[1]*y_vector.c[2]*z_vector.c[0] + x_vector.c[2]*y_vector.c[0]*z_vector.c[1] - x_vector.c[0]*y_vector.c[2]*z_vector.c[1] - x_vector.c[1]*y_vector.c[0]*z_vector.c[2] + x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
}
void InitialisationGeometry::Angle_Between_Dislocation_Axis_And_Calculation_Axis(std::vector <Vector<double>> direction_vector_list){
    double kappa_calculate,phi_calculate;
    for(int i=0;i<direction_vector_list.size();i++)
    {
        if(!is_Right_Normal(direction_vector_list[i], z_vector)){
            direction_vector_list[i]=Vector_Inverse(direction_vector_list[i]);
        }
        Vector <double> vector_normal_tayz=Vector_Multiplication<double, double>(direction_vector_list[i], z_vector);
        Vector <double> vector_normal_xz=Vector_Multiplication<double, double>(x_vector, z_vector);
        if(Vector_Absolute_Value(vector_normal_tayz)==0 || Vector_Absolute_Value(vector_normal_tayz)==0){
            phi_calculate=0;
        }
        else{
            phi_calculate=acos(Scalar_Multiplication(vector_normal_tayz, vector_normal_xz)/(Vector_Absolute_Value(vector_normal_xz)*Vector_Absolute_Value(vector_normal_tayz)));
            double angle1, angle2;
            angle1=acos(Scalar_Multiplication(z_vector, direction_vector_list[i])/(Vector_Absolute_Value(z_vector)*Vector_Absolute_Value(direction_vector_list[i])));
            angle2=acos(Scalar_Multiplication(y_vector, direction_vector_list[i])/(Vector_Absolute_Value(y_vector)*Vector_Absolute_Value(direction_vector_list[i])));
            if(((angle1 > 0 && angle1 < PI_2) && (angle2 > PI_2 && angle2 < PI)) || ((angle2 > 0 && angle2 < PI_2) && (angle1 > PI_2 && angle1 < PI))){
                phi_calculate=2*PI-phi_calculate;
            }
        }
        kappa_calculate=acos(Scalar_Multiplication(z_vector, direction_vector_list[i])/(Vector_Absolute_Value(z_vector)*Vector_Absolute_Value(direction_vector_list[i])));
        if(kappa_calculate>PI_2 && kappa_calculate<PI){
            kappa_calculate=kappa_calculate-PI_2;
        }
        if(kappa_calculate>=PI){
            kappa_calculate=kappa_calculate-PI;
        }
        kappa.push_back(kappa_calculate);
        phi.push_back(phi_calculate);
    }
}
template <typename type> bool InitialisationGeometry::is_Right_Normal(Vector<type> first_vector, Vector<type> second_vector){
    double angle = acos(Scalar_Multiplication(first_vector, second_vector)/(Vector_Absolute_Value(first_vector)*Vector_Absolute_Value(second_vector)));
    if (angle >= 0 && angle <= PI_2){
        return true;
    }
    else{
        return false;
    }
}
void InitialisationGeometry::Segment_Lenght(std::vector <double> segment_lenght){
    Vector<double> zeros;
    zeros.c[0]=0;
    zeros.c[1]=0;
    zeros.c[2]=0;
    for(int i=0;i<segment_lenght.size();i++){
        this->segment_lenght.push_back(zeros);
    }
    std::vector <int> parallel_segment_number;
    std::vector <int> matcher;
    for(int i=0;i<segment_lenght.size();i++){
        if (exit_point[i].c[0] == INF && exit_point[i].c[1] == INF){
            parallel_segment_number.push_back(i);
        }
    }
    for(int i=0;i<parallel_segment_number.size();i++){
        this->segment_lenght[parallel_segment_number[i]].c[0]=-segment_lenght[parallel_segment_number[i]]*cos(phi[parallel_segment_number[i]]);
        this->segment_lenght[parallel_segment_number[i]].c[1]=-segment_lenght[parallel_segment_number[i]]*sin(phi[parallel_segment_number[i]]);
    }
    for(int i=0;i<segment_lenght.size();i++){
        for(int p=0;p<parallel_segment_number.size();p++){
            if(segment_lenght[parallel_segment_number[p]]==segment_lenght[i]){
                this->segment_lenght[i].c[0]=this->segment_lenght[parallel_segment_number[p]].c[0];
                this->segment_lenght[i].c[1]=this->segment_lenght[parallel_segment_number[p]].c[1];
            }
        }
    }
}
