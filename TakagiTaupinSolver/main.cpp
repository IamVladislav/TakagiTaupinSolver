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

class InitialisationGeometry{
public:
    Vector<double> x_vector,y_vector,z_vector, glide_plane_vector, b_vector;
    std::vector <double> kappa, phi;
    std::vector <Vector<double>> exit_point;
    std::vector <DisplacmentGradientSystemReplace> final_model;
    double dislocation_depth, nu;
    InitialisationGeometry(Vector<double> diffraction_vector,Vector<double> normal_vector,std::vector <Vector<double>> direction_vector_list,std::vector <double> segment_lenght,  double dislocation_depth, double nu, Vector<double> burgers_vector, int nubmer_of_dislocation){
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
        Model_Creation(direction_vector_list, segment_lenght);
    }
    void Model_Creation(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght){//mock
        for(int i=0;i<direction_vector_list.size();i++){
            if(exit_point[i].c[0]>=(INF-0.01*INF) && exit_point[i].c[0]<=(INF+0.01*INF)){
                DisplacmentGradientSystemReplace angular_dislocation(this->nu, this->dislocation_depth, this->phi[i], this->kappa[i], this->exit_point[i], segment_lenght[i],this->b_vector, 1);
                final_model.push_back(angular_dislocation);
                std::cout<<"Parallel segment created!"<<std::endl;
            }
            else{
                DisplacmentGradientSystemReplace angular_dislocation(this->nu, this->dislocation_depth, this->phi[i], this->kappa[i], this->exit_point[i], segment_lenght[i], this->b_vector, 1);
                final_model.push_back(angular_dislocation);
                DisplacmentGradientSystemReplace beam_dislocation(this->nu, this->dislocation_depth, this->phi[i], this->kappa[i], this->exit_point[i], segment_lenght[i], this->b_vector, 2);
                final_model.push_back(beam_dislocation);
                std::cout<<"Exit segment created!"<<std::endl;
            }
        }
    }
    Model_Of_Polygonal_Dislocation* ModelOutput(){//mock
        Model_Of_Polygonal_Dislocation* model = new Model_Of_Polygonal_Dislocation(final_model);
        return model;
    }
    void Glide_Plane(std::vector <Vector<double>> direction_vector_list){//untested
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
    void Glide_Plane_Calculate(std::vector <Vector<double>> direction_vector_list){//untested
        Vector<double> comparison_vector;
        for(int i=1;i<direction_vector_list.size();i++){
            glide_plane_vector=Vector_Multiplication<double, double>(direction_vector_list[i-1], direction_vector_list[i]);
            if (i!=1 and comparison_vector!=glide_plane_vector){
                throw "This vectors haven't overall glide plane";
            }
            comparison_vector=glide_plane_vector;
        }
    }
    void Exit_Point_Coordinate(std::vector <Vector<double>> direction_vector_list, std::vector <double> segment_lenght){//untested
        Vector<double> exit_point_calculation;
        for(int i=0;i<direction_vector_list.size();i++)
        {
            exit_point_calculation.c[2]=-dislocation_depth;
            if(kappa[i]>=(-0.01*PI_2) && kappa[i]<=(+0.01*PI_2)){
                exit_point_calculation.c[0]=(segment_lenght[i])*cos(phi[i]);
                exit_point_calculation.c[1]=(segment_lenght[i])*sin(phi[i]);
            }
            else{
                if(kappa[i]>=(PI_2-0.005*PI_2) && kappa[i]<=(PI_2+0.005*PI_2)){
                    exit_point_calculation.c[0]=INF;
                    exit_point_calculation.c[1]=INF;
                }
                else{
                    exit_point_calculation.c[0]=(dislocation_depth/tan(kappa[i])+segment_lenght[i])*cos(phi[i]);
                    exit_point_calculation.c[1]=(dislocation_depth/tan(kappa[i])+segment_lenght[i])*sin(phi[i]);
                }
            }
            exit_point.push_back(exit_point_calculation);
        }
    }
    void Burgers_Vector_Into_System_Coordinates(Vector<double> b_initial){//not tested
        b_vector.c[0]=-((-(b_initial.c[2]*y_vector.c[1]*z_vector.c[0]) + b_initial.c[1]*y_vector.c[2]*z_vector.c[0] + b_initial.c[2]*y_vector.c[0]*z_vector.c[1] - b_initial.c[0]*y_vector.c[2]*z_vector.c[1] - b_initial.c[1]*y_vector.c[0]*z_vector.c[2] + b_initial.c[0]*y_vector.c[1]*z_vector.c[2])/ (x_vector.c[2]*y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[1]=-((b_initial.c[2]*x_vector.c[1]*z_vector.c[0] - b_initial.c[1]*x_vector.c[2]*z_vector.c[0] - b_initial.c[2]*x_vector.c[0]*z_vector.c[1] + b_initial.c[0]*x_vector.c[2]*z_vector.c[1] + b_initial.c[1]*x_vector.c[0]*z_vector.c[2] - b_initial.c[0]*x_vector.c[1]*z_vector.c[2])/(x_vector.c[2]* y_vector.c[1]*z_vector.c[0] - x_vector.c[1]*y_vector.c[2]*z_vector.c[0] - x_vector.c[2]*y_vector.c[0]*z_vector.c[1] + x_vector.c[0]*y_vector.c[2]*z_vector.c[1] + x_vector.c[1]*y_vector.c[0]*z_vector.c[2] - x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
        b_vector.c[2]=-((b_initial.c[2]*x_vector.c[1]*y_vector.c[0] - b_initial.c[1]*x_vector.c[2]*y_vector.c[0] - b_initial.c[2]*x_vector.c[0]*y_vector.c[1] + b_initial.c[0]*x_vector.c[2]*y_vector.c[1] + b_initial.c[1]*x_vector.c[0]*y_vector.c[2] - b_initial.c[0]*x_vector.c[1]* y_vector.c[2])/(-(x_vector.c[2]*y_vector.c[1]*z_vector.c[0]) + x_vector.c[1]*y_vector.c[2]*z_vector.c[0] + x_vector.c[2]*y_vector.c[0]*z_vector.c[1] - x_vector.c[0]*y_vector.c[2]*z_vector.c[1] - x_vector.c[1]*y_vector.c[0]*z_vector.c[2] + x_vector.c[0]*y_vector.c[1]*z_vector.c[2]));
    }
    void Angle_Between_Dislocation_Axis_And_Calculation_Axis(std::vector <Vector<double>> direction_vector_list){//untested
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
