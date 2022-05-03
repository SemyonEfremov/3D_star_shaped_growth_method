#include <iostream>
#include <array>
#include <vector>
#include <queue>
#include <list>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <omp.h>
#include <stdio.h>
#include "parameters.h"
#include "MetaDataWrite.h"
#include "DistanceModel.h"
#include "GridComputationMethods.h"
#include "CellGrowthMethod.h"



int main()
{
    double rotational_angle_phi = -pi / 4.0, rotational_angle_theta = -acos(1.0 / sqrt(3.0));
    clock_t starting_point, ending_point;
    time_t start, end;
    starting_point = clock();
    time(&start);
    //regular_tetrahedron_vertices = vertices_initialization(periodic_distance_type, periodic_distance_symmetry);
    //DistancePrintToFile(max_length_value, max_dimension_value, 200, 100, regular_distance_factor_const, semiregular_distance_factor_const, seed_distribution_type_const, distance_type_const, periodic_distance_type_const, periodic_distance_symmetry_const, folder, "example_distance");
    //PrintMultiDistInterpolationToFile(50, 25, 25);
    PrintPolyhedron(50, periodic_distance_type_const, periodic_distance_symmetry_const);
    for (int index_rot_theta = 0; index_rot_theta < rotation_discretization; index_rot_theta++)
    {
        for (int index_rot_phi = 0; index_rot_phi < rotation_discretization; index_rot_phi++)
        {
            rotational_angle_theta = 0.0;//(double)index_rot_theta * 6.28 / (rotation_discretization - 1);
            rotational_angle_phi = 0.0;//(double)index_rot_phi * 6.28 / (rotation_discretization - 1);
            InitializeGrowth(0.0, 0.0, max_length_value, max_dimension_value, edge_radius_value, regular_distance_factor_const, semiregular_distance_factor_const, seed_distribution_type_const, distance_type_const, periodic_distance_type_const, periodic_distance_symmetry_const);
        }
    }
    //InitializeRandomExploration(max_length_value, max_dimension_value, edge_radius_value, 50, seed_distribution_type_const, distance_type_const);
    ending_point = clock();
    time(&end);

    std::cout << "The evaluation time is " << (double)(ending_point - starting_point) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
    std::cout << "The evaluation time is " << (double)(difftime(end, start)) << " seconds." << std::endl;



    return 0;
}
