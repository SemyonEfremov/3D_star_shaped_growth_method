#pragma once
#include <vector>
#include <list>

std::string NameInitialisation(double phi, double theta, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, std::string seed_distribution, std::string distance_type, std::string periodic_distance_type, std::string periodic_distance_symmetry);
std::vector<unsigned int> RowToColumnMajor(std::vector<unsigned int>& row_data, int size_x, int size_y, int size_z);
void PrintToFileLongInt(std::list<std::vector<long int>>& cells_distribution, std::string& name);
void PrintToFileUnsignedInt(std::list<std::vector<unsigned int>>& cells_distribution, const std::string& folder, std::string& name);
void write_ITK_metaimage(const std::vector<unsigned int>& volume, const std::string& folder, const std::string& name, const unsigned int dim_x, const unsigned int dim_y, const unsigned int dim_z);
void WriteMetaData(std::list<std::vector<unsigned int>>& structure, const std::string& folder, std::string& name, int size_x, int size_y, int size_z);
void SaveToSTL(std::vector<std::array<double, 3>>& vertices, std::vector<std::array<unsigned int, 3>>& triangles, const std::string& folder, const std::string& name);