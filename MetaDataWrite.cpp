#include <list>
#include <vector>
#include <string>
#include <ios>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <omp.h>
#include "parameters.h"

std::string NameInitialisation(double phi, double theta, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, std::string seed_distribution, std::string distance_type, std::string periodic_distance_type, std::string periodic_distance_symmetry)
{
    std::string result("");
    std::string temp("");
    std::ostringstream file_id_stream;

    file_id_stream << phi << "_" << theta << "_" << edge_radius << "_" << regular_distance_factor << "_" << semiregular_distance_factor;
    result = file_id_stream.str();
    file_id_stream.str("");
    file_id_stream.clear();

    for (unsigned int index = 0; index < seed_distribution_type_list.size() ; index++)
    {
        if (seed_distribution == seed_distribution_type_list[index])
        {
            file_id_stream << index;
            temp = file_id_stream.str();
            result = result + "_" + temp;
            file_id_stream.str("");
            file_id_stream.clear();
        }
    }
    for (unsigned int index = 0; index < distance_type_list.size(); index++)
    {
        if (distance_type == distance_type_list[index])
        {
            file_id_stream << index;
            temp = file_id_stream.str();
            result = result + "_" + temp;
            file_id_stream.str("");
            file_id_stream.clear();
        }
    }
    for (unsigned int index = 0; index < periodic_distance_type_list.size(); index++)
    {
        if (periodic_distance_type == periodic_distance_type_list[index])
        {
            file_id_stream << index;
            temp = file_id_stream.str();
            result = result + "_" + temp;
            file_id_stream.str("");
            file_id_stream.clear();
        }
    }
    for (unsigned int index = 0; index < periodic_distance_symmetry_list.size(); index++)
    {
        if (periodic_distance_symmetry == periodic_distance_symmetry_list[index])
        {
            file_id_stream << index;
            temp = file_id_stream.str();
            result = result + "_" + temp;
            file_id_stream.str("");
            file_id_stream.clear();
        }
    }
    
    return(result);
}

std::vector<unsigned int> RowToColumnMajor(std::vector<unsigned int>& row_data, int size_x, int size_y, int size_z)
{
    std::cout << "switching the solution to the column major format..." << std::endl << std::endl;
    std::vector<unsigned int> column_data;
    int i, j, k; // indexes correspond to x, y, z, respectively
    int row_maj_index;

    //#pragma omp parallel for
    for (int col_maj_index = 0; col_maj_index < size_x * size_y * size_z; col_maj_index++)
    {
        i = (int)(col_maj_index % size_x);
        j = (int)((int)(col_maj_index / size_x) % size_y);
        k = (int)((int)(col_maj_index / size_x) / size_y);
        row_maj_index = size_y * size_z * i + size_z * j + k;

        column_data.push_back(row_data[row_maj_index]);
    }

    return(column_data);
}

void PrintToFileLongInt(std::list<std::vector<long int>>& cells_distribution, std::string& name)
{
    std::vector<long int> solution_1 = cells_distribution.front();
    std::vector<long int> solution_2 = cells_distribution.back();

    std::string header_1("closed_cell_structure_");
    std::string header_2("open_cell_structure_");
    std::string format(".txt");

    std::string resulting_name_1 = header_1 + name + format;
    std::string resulting_name_2 = header_2 + name + format;

    std::ofstream out_cells_id;
    out_cells_id.open(resulting_name_1);
    if (out_cells_id.is_open())
    {
        for (auto element : solution_1)
        {
            out_cells_id << element << std::endl;
        }
    }

    out_cells_id.close();

    std::ofstream out_seeds_id;
    out_seeds_id.open(resulting_name_2);
    if (out_seeds_id.is_open())
    {
        for (auto element : solution_2)
        {
            out_seeds_id << element << std::endl;
        }
    }

    out_seeds_id.close();

    std::cout << resulting_name_1 << " and " << resulting_name_2 << " were successfully created." << std::endl;
}

void PrintToFileUnsignedInt(std::list<std::vector<unsigned int>>& cells_distribution, const std::string& folder, std::string& name)
{
    std::vector<unsigned int> solution_1 = cells_distribution.front();
    std::vector<unsigned int> solution_2 = cells_distribution.back();

    std::string header_1("cl_");
    std::string header_2("op_");
    std::string format(".txt");
    header_1 = folder + std::string("/homogenization/") + header_1;
    header_2 = folder + std::string("/homogenization/") + header_2;

    std::string resulting_name_1 = header_1 + name + format;
    std::string resulting_name_2 = header_2 + name + format;

    std::ofstream out_cells_id;
    out_cells_id.open(resulting_name_1);
    if (out_cells_id.is_open())
    {
        for (auto element : solution_1)
        {
            out_cells_id << element << std::endl;
        }
    }

    out_cells_id.close();

    std::ofstream out_seeds_id;
    out_seeds_id.open(resulting_name_2);
    if (out_seeds_id.is_open())
    {
        for (auto element : solution_2)
        {
            out_seeds_id << element << std::endl;
        }
    }

    out_seeds_id.close();

    std::cout << resulting_name_1 << " and " << resulting_name_2 << " were successfully created." << std::endl;
}

void write_ITK_metaimage(const std::vector<unsigned int>& volume, const std::string& folder, const std::string& name, const unsigned int dim_x, const unsigned int dim_y, const unsigned int dim_z)
{
    /**
    Writes a ITK metaimage file, which can be viewed by Paraview.
    Important: Assume the binary data will be stored as INT.
    See http://www.itk.org/Wiki/ITK/MetaIO/Documentation

    Generates a raw file containing the volume, and an associated mhd
    metaimage file.

    Assumes volume is column-major order (x increasing fastest).
    */

    std::ofstream volume_file(folder + "/" + name + ".raw", std::ios::out | std::ios::binary);
    volume_file.write((char*)volume.data(), dim_x * dim_y * dim_z * sizeof(int));
    volume_file.close();
    std::ofstream mhd_file;
    mhd_file.open(folder + "/" + name + ".mhd");
    std::string dim_x_str = std::to_string(dim_x);
    std::string dim_y_str = std::to_string(dim_y);
    std::string dim_z_str = std::to_string(dim_z);
    mhd_file << "ObjectType = Image\nNDims = 3\nDimSize = " + dim_x_str + " " + dim_y_str + " " + dim_z_str + " " + "\nElementType = MET_INT\nElementDataFile = " + name + ".raw";
    mhd_file.close();
}

void WriteMetaData(std::list<std::vector<unsigned int>>& structure, const std::string& folder, std::string& name, int size_x, int size_y, int size_z)
{
    std::cout << "writing metadata..." << std::endl << std::endl;

    std::vector<unsigned int> row_maj_structure = structure.front();
    std::vector<unsigned int> col_maj_structure = RowToColumnMajor(row_maj_structure, size_x, size_y, size_z);

    write_ITK_metaimage(col_maj_structure, folder, "cl_" + name, size_x, size_y, size_z);

    row_maj_structure = structure.back();
    col_maj_structure = RowToColumnMajor(row_maj_structure, size_x, size_y, size_z);

    write_ITK_metaimage(col_maj_structure, folder, "op_" + name, size_x, size_y, size_z);
}

std::array<double, 3> ComputeNormal(std::array<double, 3> point_0, std::array<double, 3> point_1, std::array<double, 3> point_2)
{
    std::array<double, 3> result;
    result[0] = (point_1[1] - point_0[1]) * (point_2[2] - point_0[2]) - (point_2[1] - point_0[1]) * (point_1[2] - point_0[2]);
    result[1] = (point_2[0] - point_0[0]) * (point_1[2] - point_0[2]) - (point_1[0] - point_0[0]) * (point_2[2] - point_0[2]);
    result[2] = (point_1[0] - point_0[0]) * (point_2[1] - point_0[1]) - (point_2[0] - point_0[0]) * (point_1[1] - point_0[1]);

    return(result);
}

void SaveToSTL(std::vector<std::array<double, 3>>& vertices, std::vector<std::array<unsigned int, 3>>& triangles, const std::string& folder, const std::string& name)
{
    std::array<double, 3> point_0, point_1, point_2, normal_local;
    
    std::string resulting_name = folder + std::string("/") + name + std::string(".stl");

    std::ofstream surface_out;
    surface_out.open(resulting_name);

    surface_out << "solid STL generated by CustomCode" << std::endl;
    for (auto const& temp_triangle : triangles)
    {
        point_0 = vertices[(int)temp_triangle[0]];
        point_1 = vertices[(int)temp_triangle[1]];
        point_2 = vertices[(int)temp_triangle[2]];
        normal_local = ComputeNormal(point_0, point_1, point_2);
        surface_out << "facet normal" << "  " << normal_local[0] << "   " << normal_local[1] << "   " << normal_local[2] << std::endl;
        surface_out << "    outer loop" << std::endl;
        surface_out << "        vertex" << "    " << point_0[0] << "    " << point_0[1] << "    " << point_0[2] << std::endl;
        surface_out << "        vertex" << "    " << point_1[0] << "    " << point_1[1] << "    " << point_1[2] << std::endl;
        surface_out << "        vertex" << "    " << point_2[0] << "    " << point_2[1] << "    " << point_2[2] << std::endl;
        surface_out << "    endloop" << std::endl;
        surface_out << "endfacet" << std::endl;
    }

    surface_out.close();
}