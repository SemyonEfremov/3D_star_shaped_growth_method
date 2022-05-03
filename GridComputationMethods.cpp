#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include "parameters.h"






std::vector<unsigned int> StructureExpansion(std::vector<unsigned int>& structure, unsigned int expansion_ratio_x, unsigned int expansion_ratio_y, unsigned int expansion_ratio_z, unsigned int nx, unsigned int ny, unsigned int nz)
{
    std::vector<unsigned int> result;

    for (unsigned int iteration_x = 0; iteration_x < expansion_ratio_x; iteration_x++)
    {
        for (unsigned int i = 0; i < nx; i++)
        {
            for (unsigned int iteration_y = 0; iteration_y < expansion_ratio_y; iteration_y++)
            {
                for (unsigned int j = 0; j < ny; j++)
                {
                    for (unsigned int iteration_z = 0; iteration_z < expansion_ratio_z; iteration_z++)
                    {
                        for (unsigned int k = 0; k < nz; k++)
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * i + nz * j + k)]);
                        }
                    }
                }
            }
        }
    }

    return(result);
}

std::vector<unsigned int> StructureEdging(std::vector<unsigned int>& structure, unsigned int nx, unsigned int ny, unsigned int nz)
{
    std::vector<unsigned int> result;

    for (unsigned int i = 0; i < nx; i++)
    {
        if (i == 0)
        {
            for (unsigned int j = 0; j < (ny + 2); j++)
            {
                for (unsigned int k = 0; k < (nz + 2); k++)
                {
                    result.push_back(0);
                }
            }
        }
        for (unsigned int j = 0; j < ny; j++)
        {
            if (j == 0)
            {
                for (unsigned int k = 0; k < (nz + 2); k++)
                {
                    result.push_back(0);
                }
            }
            for (unsigned int k = 0; k < nz; k++)
            {
                if (k == 0) { result.push_back(0); }
                result.push_back(structure[(unsigned int)(ny * nz * i + nz * j + k)]);
                if (k == (nz - 1)) { result.push_back(0); }
            }
            if (j == (ny - 1))
            {
                for (unsigned int k = 0; k < (nz + 2); k++)
                {
                    result.push_back(0);
                }
            }
        }
        if (i == (nx - 1))
        {
            for (unsigned int j = 0; j < (ny + 2); j++)
            {
                for (unsigned int k = 0; k < (nz + 2); k++)
                {
                    result.push_back(0);
                }
            }
        }
    }

    return(result);
}

std::vector<unsigned int> StructureEdgingManual(std::vector<unsigned int>& structure, unsigned int extra_voxels_x, unsigned int extra_voxels_y, unsigned int extra_voxels_z, unsigned int nx, unsigned int ny, unsigned int nz)
{
    std::vector<unsigned int> result;

    for (unsigned int i = 0; i < nx; i++)
    {
        if (i == 0)
        {
            for (unsigned int l = 0; l < extra_voxels_x; l++)
            {
                for (int j = -(int)extra_voxels_y; j < (int)(ny + extra_voxels_y); j++)
                {
                    for (int k = -(int)extra_voxels_z; k < (int)(nz + extra_voxels_z); k++)
                    {
                        if (((j >= 0) && (j < ny)) && ((k >= 0) && (k < nz)))
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * (nx - 1 - (extra_voxels_x - 1 - l)) + nz * j + k)]);
                        }
                        else
                        {
                            if (!((j >= 0) && (j < ny)) && ((k >= 0) && (k < nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * (nx - 1 - (extra_voxels_x - 1 - l)) + nz * fabs(ny - fabs(j)) + k)]);
                            }
                            if (((j >= 0) && (j < ny)) && !((k >= 0) && (k < nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * (nx - 1 - (extra_voxels_x - 1 - l)) + nz * j + fabs(nz - fabs(k)))]);
                            }
                            if (!((j >= 0) && (j < ny)) && !((k >= 0) && (k < nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * (nx - 1 - (extra_voxels_x - 1 - l)) + nz * fabs(ny - fabs(j)) + fabs(nz - fabs(k)))]);
                            }
                        }
                    }
                }
            }
            //std::cout << "x passed " << std::endl;
        }
        for (unsigned int j = 0; j < ny; j++)
        {
            if (j == 0)
            {
                for (unsigned int l = 0; l < extra_voxels_y; l++)
                {
                    for (int k = -(int)extra_voxels_z; k < (int)(nz + extra_voxels_z); k++)
                    {
                        if ((k >= 0) && (k < nz))
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * i + nz * (ny - 1 - (extra_voxels_y - 1 - l)) + k)]);
                        }
                        else
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * i + nz * (ny - 1 - (extra_voxels_y - 1 - l)) + fabs(nz - fabs(k)))]);
                        }
                    }
                }
                //std::cout << "y passed" << std::endl;
            }
            for (unsigned int k = 0; k < nz; k++)
            {
                if (k == 0)
                {
                    for (unsigned int l = 0; l < extra_voxels_z; l++)
                    {
                        result.push_back(structure[(unsigned int)(ny * nz * i + nz * j + (nz - 1 - (extra_voxels_z - 1 - l)))]);
                    }
                    //std::cout << "z passed" << std::endl;
                }
                result.push_back(structure[(unsigned int)(ny * nz * i + nz * j + k)]);
                if (k == (nz - 1))
                {
                    for (unsigned int l = 0; l < extra_voxels_z; l++)
                    {
                        result.push_back(structure[(unsigned int)(ny * nz * i + nz * j + l)]);
                    }
                    //std::cout << "z passed" << std::endl;
                }
            }
            if (j == (ny - 1))
            {
                for (unsigned int l = 0; l < extra_voxels_y; l++)
                {
                    for (int k = -(int)extra_voxels_z; k < (int)(nz + extra_voxels_z); k++)
                    {
                        if ((k >= 0) && (k < nz))
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * i + nz * l + k)]);
                        }
                        else
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * i + nz * l + fabs(nz - fabs(k)))]);
                        }
                    }
                }
                //std::cout << "y passed" << std::endl;
            }
        }
        if (i == (nx - 1))
        {
            for (unsigned int l = 0; l < extra_voxels_x; l++)
            {
                for (int j = -(int)extra_voxels_y; j < (int)(ny + extra_voxels_y); j++)
                {
                    for (int k = -(int)extra_voxels_z; k < (int)(nz + extra_voxels_z); k++)
                    {
                        if (((j >= 0) && (j < ny)) && ((k >= 0) && (k < nz)))
                        {
                            result.push_back(structure[(unsigned int)(ny * nz * l + nz * j + k)]);
                        }
                        else
                        {
                            if (((j < 0) || (j >= ny)) && ((k >= 0) && (k < nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * l + nz * fabs(ny - fabs(j)) + k)]);
                            }
                            if (((j >= 0) && (j < ny)) && ((k < 0) || (k >= nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * l + nz * j + fabs(nz - fabs(k)))]);
                            }
                            if (((j < 0) || (j >= ny)) && ((k < 0) || (k >= nz)))
                            {
                                result.push_back(structure[(unsigned int)(ny * nz * l + nz * fabs(ny - fabs(j)) + fabs(nz - fabs(k)))]);
                            }
                        }
                    }
                }
            }
            //std::cout << "x passed" << std::endl;
        }
    }

    return(result);
}

void CellGridExpansion(std::vector<std::array<double, 3>>& seeds_grid_sample, std::array<double, 3>& sample_length, std::array<int, 3>& sample_size)
{
    std::vector<std::array<double, 3>> temp_seeds_grid_sample;
    std::array<double, 3> seed;

    double min_expansion = std::min(expansion_rate_y * sample_length[1], std::min(expansion_rate_x * sample_length[0], expansion_rate_z * sample_length[2]));
    double expansion_ratio_x = (expansion_rate_x * sample_length[0]) / min_expansion;
    double expansion_ratio_y = (expansion_rate_y * sample_length[1]) / min_expansion;
    double expansion_ratio_z = (expansion_rate_z * sample_length[2]) / min_expansion;
    sample_length[0] = expansion_ratio_x;
    sample_length[1] = expansion_ratio_y;
    sample_length[2] = expansion_ratio_z;
    sample_size[0] = (int)(expansion_rate_x * (sample_size[0] - 1) + 1);
    sample_size[1] = (int)(expansion_rate_y * (sample_size[1] - 1) + 1);
    sample_size[2] = (int)(expansion_rate_z * (sample_size[2] - 1) + 1);
    
    for (const auto& grid_element : seeds_grid_sample)
    {
        for (unsigned int i = 0; i < expansion_rate_x; i++)
        {
            seed[0] = (grid_element[0] + i * (min_expansion * sample_length[0]) / expansion_rate_x) / min_expansion;
            for (unsigned int j = 0; j < expansion_rate_y; j++)
            {
                seed[1] = (grid_element[1] + j * (min_expansion * sample_length[1]) / expansion_rate_y) / min_expansion;
                for (unsigned int k = 0; k < expansion_rate_z; k++)
                {
                    seed[2] = (grid_element[2] + k * (min_expansion * sample_length[2]) / expansion_rate_z) / min_expansion;

                    temp_seeds_grid_sample.push_back(seed);
                }
            }
        }
    }

    seeds_grid_sample = temp_seeds_grid_sample;

    std::ofstream myfile_grid;
    myfile_grid.open("example_extended_grid.csv");
    myfile_grid << "X, Y, Z\n";
    for (const auto& point : seeds_grid_sample)
    {
        myfile_grid << point[0] << ", " << point[1] << ", " << point[2] << "\n";
    }
    myfile_grid.close();
}

void CellGridInitialization(std::vector<std::array<double, 3>>& seeds_grid_sample, std::array<double, 3>& sample_length, std::array<int, 3>& sample_size, double min_length, int min_dimension, const std::string& grid_type, const std::string& distance_type)
{
    std::array<double, 3> temp_point_coordinates;

    double l = 0.0, l1 = 0.0;

    if (grid_type == std::string("triclinic"))
    {
        l = 4.0 * min_length / 3.0;

        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.4 * l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l / 3.0, 0.8 * l / 3.0, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { (1.0 / 3.0 + 0.5) * l, (0.8 / 3.0 + 0.4) * l, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 2.0 * l / 3.0, 1.6 * l / 3.0, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { (2.0 / 3.0 - 0.5) * l, (1.6 / 3.0 + 0.4) * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(l * ((double)min_dimension - 1) + 1), (int)(0.8 * l * ((double)min_dimension - 1) + 1), min_dimension };
        sample_length = { l, 0.8 * l, 3.0 * l / 4.0 };
    }
    
    if (grid_type == std::string("cubic_primitive"))
    {
        temp_point_coordinates = { 0.5 * min_length, 0.5 * min_length, 0.5 * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { min_length, min_length, min_length };
    }

    if (grid_type == std::string("hexagonal"))
    {
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = {min_length, sqrt(3.0) * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(2.0 * ((double)min_dimension - 1) + 1), (int)(2.0 * sqrt(3.0) * ((double)min_dimension - 1) + 1), min_dimension };
        //sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { 2.0 * min_length, 2.0 * sqrt(3.0) * min_length, min_length };
    }

    if (grid_type == std::string("trigonal"))
    {
        double k = 4.0;

        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 0.5 * sqrt(3.0) * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, min_length / (2.0 * sqrt(3.0)), k * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 2.0 * min_length / sqrt(3.0), k * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0,  min_length / sqrt(3.0), k * 2.0 * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 5.0 * sqrt(3.0) * min_length / 6.0, k * 2.0 * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, (int)(sqrt(3.0) * ((double)min_dimension - 1) + 1), (int)(k * ((double)min_dimension - 1) + 1) };
        sample_length = { min_length, sqrt(3.0) * min_length, k * min_length };
    }

    if (grid_type == std::string("tetragonal_primitive"))
    {
        l = 3.0;

        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { min_length, min_length, l * min_length };
    }
    
    if (grid_type == std::string("tetragonal_body-centered"))
    {
        l = 3.0;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 0.5 * min_length, 0.5 * l * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { min_length, min_length, l * min_length };
    }

    if (grid_type == std::string("orthorhombic_primitive"))
    {
        l = 4.0;
        l1 = 2.5;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(l1 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { l1 * min_length, min_length, l * min_length };
    }

    if (grid_type == std::string("orthorhombic_base-centered"))
    {
        l = 4.0;
        l1 = 2.5;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l1 * min_length, 0.5 * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(l1 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { l1 * min_length, min_length, l * min_length };
    }

    if (grid_type == std::string("orthorhombic_body-centered"))
    {
        l = 4.0;
        l1 = 2.5;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l1 * min_length, 0.5 * min_length, 0.5 * l * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(l1 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { l1 * min_length, min_length, l * min_length };
    }

    if (grid_type == std::string("orthorhombic_face-centered"))
    {
        l = 4.0;
        l1 = 2.5;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * min_length, 0.5 * l * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l1 * min_length, 0.0, 0.5 * l * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l1 * min_length, 0.5 * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(l1 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = { l1 * min_length, min_length, l * min_length };
    }

    if (grid_type == std::string("monoclinic_primitive"))
    {
        l = 3.0 / 2.0;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 18.0, 0.0, l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 9.0, 0.0, 2.0 * l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(5.0 * l / 6.0 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = {5.0 * l * min_length / 6.0, min_length, l * min_length };
    }

    if (grid_type == std::string("monoclinic_base-centered"))
    {
        l = 3.0 / 2.0;
        
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 18.0, 0.0, l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 9.0, 0.0, 2.0 * l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 12.0, 2.0 * l * min_length / 6.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 25.0 * l * min_length / 36.0, 2.0 * l * min_length / 6.0, l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 5.0 * l * min_length / 36.0, 2.0 * l * min_length / 6.0, 2.0 * l * min_length / 3.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)(5.0 * l / 6.0 * ((double)min_dimension - 1) + 1), min_dimension, (int)(l * ((double)min_dimension - 1) + 1) };
        sample_length = {5.0 * l * min_length / 6.0, min_length, l * min_length };
    }

    if (grid_type == std::string("bitruncated_cubic"))
    {
        l = min_length / (2.0 * sqrt(2.0));
        temp_point_coordinates = { l / sqrt(2.0), 0.5 * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.5 * min_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, l / sqrt(2.0), 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 3.0 * l / sqrt(2.0), 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 0.0, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * min_length, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * min_length, 0.0, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * min_length, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l / sqrt(2.0), 0.0, 0.5 * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.0, 0.5 * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l / sqrt(2.0), 0.5 * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 3.0 * l / sqrt(2.0), 0.5 * min_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { min_length, min_length, min_length };
    }

    if (grid_type == std::string("cubic_body-centered"))
    {
        l = min_length;
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { min_length, min_length, min_length };
    }

    if (grid_type == std::string("cubic_face-centered"))
    {
        l = min_length;
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.0, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { min_length, min_length, min_length };
    }

    if (grid_type == std::string("diamond_cubic"))
    {
        l = min_length;
        temp_point_coordinates = { 0.75 * l, 0.25 * l, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * l, 0.75 * l, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * l, 0.25 * l, 0.75 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * l, 0.75 * l, 0.75 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0 * l, 0.0 * l, 1.0 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.0 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0 * l, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.0 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { min_dimension, min_dimension, min_dimension };
        sample_length = { min_length, min_length, min_length };
    }

    if (distance_type == "multiple")
    {
        std::cout << "...seeds grid expansion...\n" << std::endl;
        CellGridExpansion(seeds_grid_sample, sample_length, sample_size);
    }
}

std::vector<double> GridComputation(double length, int size)
{
    std::vector<double> grid;

    double step = length / (size - 1);

    for (int i = 0; i < size; i++)
    {
        grid.push_back(i * step);
    }

    std::cout << "grid dimension step is " << grid[1] - grid[0] << std::endl;

    return grid;
}

std::vector<double> GridCentersComputation(std::vector<double>& x)
{
    std::vector<double> centers;

    for (int i = 0; i < x.size() - 1; i++)
    {
        centers.push_back((x[i + 1] + x[i]) / 2.0);
    }

    return centers;
}

std::vector<std::array<double, 3>> PrintGridToFile(double max_length, int max_dimension, const std::string& grid_type)
{
    std::vector<std::array<double, 3>> seeds_grid_sample;
    std::array<double, 3> temp_point_coordinates, sample_length;
    std::array<int, 3> sample_size;

    double l = 0.0;

    if (grid_type == std::string("cubic"))
    {
        temp_point_coordinates = { 0.25 * max_length, 0.25 * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.75 * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.25 * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.75 * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.25 * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.75 * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.25 * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.75 * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, max_length, max_length };
    }

    if (grid_type == std::string("triangular_prismatic"))
    {
        temp_point_coordinates = { 0.0, 0.0, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 0.0, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.25 * sqrt(3.0) * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.25 * sqrt(3.0) * max_length, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.0, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5, 0.0, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.25 * sqrt(3.0) * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.25 * sqrt(3.0) * max_length, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, (int)(0.5 * sqrt(3.0) * ((double)max_dimension - 1) + 1), max_dimension };
        //sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, 0.5 * sqrt(3.0) * max_length, max_length };
    }

    if (grid_type == std::string("hexagonal_prismatic"))
    {
        l = max_length / (1.0 + sqrt(3.0));

        temp_point_coordinates = { 0.5 * l, 0.0, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 1.5 * l, 0.0, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * sqrt(3.0) * l, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * (2.0 + sqrt(3.0)) * l, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l,  0.5 * sqrt(3.0) * l, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.5 * (2.0 + sqrt(3.0)) * l, 0.25 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.0, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 1.5 * l, 0.0, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * sqrt(3.0) * l, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * (2.0 + sqrt(3.0)) * l, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l,  0.5 * sqrt(3.0) * l, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.5 * (2.0 + sqrt(3.0)) * l, 0.75 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { (int)((2.0 * l) * ((double)max_dimension - 1) + 1), (int)(((1.0 + sqrt(3.0)) * l) * ((double)max_dimension - 1) + 1), max_dimension };
        sample_length = { (2.0 * l) * max_length, ((1.0 + sqrt(3.0)) * l) * max_length, max_length };
    }

    if (grid_type == std::string("alternated_cubic"))
    {
        temp_point_coordinates = { 0.0, 0.0, 0.5 * max_length / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 0.0, 0.5 * max_length / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0, 0.5 * max_length, 0.5 * max_length / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 0.5 * max_length, 0.5 * max_length / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.25 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.25 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * max_length, 0.75 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * max_length, 0.75 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, (int)(1.0 / sqrt(2.0) * ((double)max_dimension - 1) + 1) };
        sample_length = { max_length, max_length, 1.0 / sqrt(2.0) * max_length };
    }

    if (grid_type == std::string("bitruncated_cubic"))
    {
        l = max_length / (2.0 * sqrt(2.0));
        temp_point_coordinates = { l / sqrt(2.0), 0.5 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.5 * max_length, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, l / sqrt(2.0), 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 3.0 * l / sqrt(2.0), 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 0.0, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * max_length, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 0.0, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * max_length, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l / sqrt(2.0), 0.0, 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.0, 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l / sqrt(2.0), 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 3.0 * l / sqrt(2.0), 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l / sqrt(2.0), 0.5 * max_length, max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.5 * max_length, max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, l / sqrt(2.0), max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, 3.0 * l / sqrt(2.0), max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, max_length, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { max_length, 0.5 * max_length, l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * max_length, max_length, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { max_length, 0.5 * max_length, 3.0 * l / sqrt(2.0) };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l / sqrt(2.0), max_length, 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 3.0 * l / sqrt(2.0), 0.0, 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { max_length, l / sqrt(2.0), 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { max_length, 3.0 * l / sqrt(2.0), 0.5 * max_length };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, max_length, max_length };
    }

    if (grid_type == std::string("body-centered_cubic"))
    {
        l = max_length;
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.0, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.0, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, l, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, max_length, max_length };
    }

    if (grid_type == std::string("face-centered_cubic"))
    {
        l = max_length;
        temp_point_coordinates = { 0.0, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.0, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, 0.0, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.0, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0, l, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.0, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, l, 0.0 };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { l, l, l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, max_length, max_length };
    }

    if (grid_type == std::string("diamond_cubic"))
    {
        l = max_length;
        temp_point_coordinates = { 0.75 * l, 0.25 * l, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * l, 0.75 * l, 0.25 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.25 * l, 0.25 * l, 0.75 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.75 * l, 0.75 * l, 0.75 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0 * l, 0.0 * l, 1.0 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.0 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.0 * l, 0.5 * l, 0.5 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        temp_point_coordinates = { 0.5 * l, 0.5 * l, 0.0 * l };
        seeds_grid_sample.push_back(temp_point_coordinates);

        sample_size = { max_dimension, max_dimension, max_dimension };
        sample_length = { max_length, max_length, max_length };
    }

    std::ofstream myfile_grid;
    myfile_grid.open("example_grid.csv");
    myfile_grid << "X, Y, Z\n";
    for (const auto& point: seeds_grid_sample)
    {
        myfile_grid << point[0] << ", " << point[1] << ", " << point[2] << "\n";
    }
    myfile_grid.close();
    
    return(seeds_grid_sample);
}