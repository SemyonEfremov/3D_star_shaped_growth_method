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



namespace PostProcessing
{
    std::array<int, 3> IndexFlatToTriple(int flat_index, int num_x, int num_y, int num_z)
    {
        std::array<int, 3> result;
        result[0] = (int)(floor(floor(flat_index / num_z) / num_y));
        result[1] = (int)((int)floor(flat_index / num_z) % num_y);
        result[2] = (int)(flat_index % num_z);

        return(result);
    }

    double SurfaceForStructure(std::array<double, 3>& point, std::array<double, 3>& central_point, const std::string& surface_type)
    {
        double result;
        double min_length = std::min(central_point[0], std::min(central_point[1], central_point[2]));

        if (surface_type == std::string("sphere"))
        {
            result = (point[0] - central_point[0]) * (point[0] - central_point[0]) + (point[1] - central_point[1]) * (point[1] - central_point[1]) + (point[2] - central_point[2]) * (point[2] - central_point[2]) - 0.8 * min_length * min_length;
        }

        return(result);
    }

    void DetectiongEdge(int flat_index, std::vector<unsigned int>& volume_data, std::vector<double>& x_coordinates, std::vector<double>& y_coordinates, std::vector<double>& z_coordinates, std::array<double, 3>& central_point, std::array<int, 3>& dimensions, double min_radius, int thickness)
    {
        bool InsideDomainIndicator = true;
        std::array<double, 3> dist_comp;
        std::array<int, 3> index = IndexFlatToTriple(flat_index, dimensions[0], dimensions[1], dimensions[2]), index_local, index_neighbor;
        
        for (int flat_index_local = 1; flat_index_local < (thickness * thickness * thickness); flat_index_local++)
        {
            if (volume_data[flat_index] != 2)
            {
                InsideDomainIndicator = true;
                index_local = IndexFlatToTriple(flat_index_local, thickness, thickness, thickness);
                index_neighbor[0] = index[0] + index_local[0];
                index_neighbor[1] = index[1] + index_local[1];
                index_neighbor[2] = index[2] + index_local[2];
                if (index_neighbor[0] < dimensions[0])
                {
                    dist_comp[0] = x_coordinates[index_neighbor[0]] - x_coordinates[index[0]];
                }
                else
                {
                    dist_comp[0] = x_coordinates[index_neighbor[0] - dimensions[0]] + x_coordinates[dimensions[0] - 1] - x_coordinates[index[0]];
                    InsideDomainIndicator = false;
                }
                if (index_neighbor[1] < dimensions[1])
                {
                    dist_comp[1] = y_coordinates[index_neighbor[1]] - y_coordinates[index[1]];
                }
                else
                {
                    dist_comp[1] = y_coordinates[index_neighbor[1] - dimensions[1]] + y_coordinates[dimensions[1] - 1] - y_coordinates[index[1]];
                    InsideDomainIndicator = false;
                }
                if (index_neighbor[2] < dimensions[2])
                {
                    dist_comp[2] = z_coordinates[index_neighbor[2]] - z_coordinates[index[2]];
                }
                else
                {
                    dist_comp[2] = z_coordinates[index_neighbor[2] - dimensions[2]] + z_coordinates[dimensions[2] - 1] - z_coordinates[index[2]];
                    InsideDomainIndicator = false;
                }
                if (SquaredEucledianDistance(dist_comp[0], dist_comp[1], dist_comp[2], 0.0, 0.0, 0.0) <= min_radius)
                {
                    if (InsideDomainIndicator == false)
                    {
                        volume_data[flat_index] = 2;
                    }
                    else
                    {
                        if (volume_data[index_neighbor[0] * dimensions[1] * dimensions[2] + index_neighbor[1] * dimensions[2] + index_neighbor[2]] == 0)
                        {
                            volume_data[flat_index] = 2;
                        }
                    }
                }
                if (volume_data[flat_index] != 2)
                {
                    InsideDomainIndicator = true;
                    index_neighbor[0] = index[0] - index_local[0];
                    index_neighbor[1] = index[1] - index_local[1];
                    index_neighbor[2] = index[2] - index_local[2];
                    if (index_neighbor[0] >= 0)
                    {
                        dist_comp[0] = x_coordinates[index_neighbor[0]] - x_coordinates[index[0]];
                    }
                    else
                    {
                        dist_comp[0] = x_coordinates[dimensions[0] - 1] - x_coordinates[dimensions[0] + index_neighbor[0] - 1] + x_coordinates[index[0]];
                        InsideDomainIndicator = false;
                    }
                    if (index_neighbor[1] >= 0)
                    {
                        dist_comp[1] = y_coordinates[index_neighbor[1]] - y_coordinates[index[1]];
                    }
                    else
                    {
                        dist_comp[1] = y_coordinates[dimensions[1] - 1] - y_coordinates[dimensions[1] + index_neighbor[1] - 1] + y_coordinates[index[1]];
                        InsideDomainIndicator = false;
                    }
                    if (index_neighbor[2] >= 0)
                    {
                        dist_comp[2] = z_coordinates[index_neighbor[2]] - z_coordinates[index[2]];
                    }
                    else
                    {
                        dist_comp[2] = z_coordinates[dimensions[2] - 1] - z_coordinates[dimensions[2] + index_neighbor[2] - 1] + z_coordinates[index[2]];
                        InsideDomainIndicator = false;
                    }
                    if (SquaredEucledianDistance(dist_comp[0], dist_comp[1], dist_comp[2], 0.0, 0.0, 0.0) <= min_radius)
                    {
                        if (InsideDomainIndicator == false)
                        {
                            volume_data[flat_index] = 2;
                        }
                        else
                        {
                            if (volume_data[index_neighbor[0] * dimensions[1] * dimensions[2] + index_neighbor[1] * dimensions[2] + index_neighbor[2]] == 0)
                            {
                                volume_data[flat_index] = 2;
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<unsigned int> ReadVolumeData(std::vector<double>& x_coordinates, std::vector<double>& y_coordinates, std::vector<double>& z_coordinates, unsigned int thickness, const std::string& surface_type)
    {
        std::vector<unsigned int> result;
        std::array<double, 3> point, central_point = { 0.5 * (x_coordinates[x_coordinates.size() - 1] + x_coordinates[0]), 0.5 * (y_coordinates[y_coordinates.size() - 1] + y_coordinates[0]) , 0.5 * (z_coordinates[z_coordinates.size() - 1] + z_coordinates[0]) };
        std::cout << "central point = { " << central_point[0] << ", " << central_point[1] << ", " << central_point[2] << " }" << std::endl;
        double min_radius = 1.01 * std::max(EucledianDistance(x_coordinates[0], y_coordinates[0], z_coordinates[0], x_coordinates[thickness], y_coordinates[0], z_coordinates[0]), std::max(EucledianDistance(x_coordinates[0], y_coordinates[0], z_coordinates[0], x_coordinates[0], y_coordinates[thickness], z_coordinates[0]), EucledianDistance(x_coordinates[0], y_coordinates[0], z_coordinates[0], x_coordinates[0], y_coordinates[0], z_coordinates[thickness])));
        min_radius = min_radius * min_radius;
        std::array<int, 3> dimensions = { x_coordinates.size(), y_coordinates.size(), z_coordinates.size() };
        std::cout << "2" << std::endl;
        for (int i = 0; i < dimensions[0]; i++)
        {
            for (int j = 0; j < dimensions[1]; j++)
            {
                for (int k = 0; k < dimensions[2]; k++)
                {
                    point = { x_coordinates[i], y_coordinates[j], z_coordinates[k] };
                    if (SurfaceForStructure(point, central_point, surface_type) <= 0.0) { result.push_back(1); }
                    else { result.push_back(0); }
                }
            }
        }
        std::cout << "1" << std::endl;

        for (int flat_index = 0; flat_index < result.size(); flat_index++)
        {
            if (result[flat_index] == 1)
            {
                DetectiongEdge(flat_index, result, x_coordinates, y_coordinates, z_coordinates, central_point, dimensions, min_radius, thickness);
            }
        }

        return(result);
    }

    std::vector<unsigned int> StructureIntoVolume(std::list<std::vector<unsigned int>>& structure, std::array<double, 3>& size, std::array<int, 3>& dimensions, unsigned int thickness, const std::string& surface_type)
    {
        std::cout << std::endl << "Fitting the structure into the shape..." << std::endl;

        std::vector<unsigned int> structure_closed = structure.front();
        std::vector<unsigned int> structure_open = structure.back();
        std::vector<unsigned int> result;

        std::vector<double> x_grid = GridCentersComputation(GridComputation(size[0], dimensions[0] + 1));
        std::vector<double> y_grid = GridCentersComputation(GridComputation(size[1], dimensions[1] + 1));
        std::vector<double> z_grid = GridCentersComputation(GridComputation(size[2], dimensions[2] + 1));
        std::cout << std::endl << "Computing volumetric data..." << std::endl;
        std::vector<unsigned int> volume = ReadVolumeData(x_grid, y_grid, z_grid, thickness, surface_type);

        std::cout << std::endl << "Computing the resulting internal structure..." << std::endl;

        std::cout << structure_open.size() << " " << volume.size() << " " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << " " << x_grid.size() << " " << y_grid.size() << " " << z_grid.size() << std::endl;

//#pragma omp parallel for
        for (int flat_index = 0; flat_index < (dimensions[0] * dimensions[1] * dimensions[2]); flat_index++)
        {
            if (volume[flat_index] != 0)
            {
                if (volume[flat_index] == 1) { result.push_back(structure_open[flat_index]); }
                else { result.push_back(structure_closed[flat_index]); }
            }
            else { result.push_back(0); }
        }

        return(result);
    }

    void SkeletonWriteToFile(std::vector<unsigned int>& identificator, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers)
    {
        int index_x = 0, index_y = 0, index_z = 0, size_x = x_centers.size(), size_y = y_centers.size(), size_z = z_centers.size();

        std::ofstream myfile_grid;
        myfile_grid.open("sceleton_points.csv");
        myfile_grid << "X, Y, Z\n";

        for (int index = 0; index < identificator.size(); index++)
        {
            if (identificator[index] == 1)
            {
                index_x = (int)(floor(floor(index / size_z) / size_y));
                index_y = (int)((int)(floor(index / size_z)) % size_y);
                index_z = (int)(index % size_z);

                myfile_grid << x_centers[index_x] << ", " << y_centers[index_y] << ", " << z_centers[index_z] << "\n";
            }
        }
    }
}