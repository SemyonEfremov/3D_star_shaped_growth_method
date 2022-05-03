#pragma once



std::vector<unsigned int> StructureExpansion(std::vector<unsigned int>& structure, unsigned int expansion_ratio_x, unsigned int expansion_ratio_y, unsigned int expansion_ratio_z, unsigned int nx, unsigned int ny, unsigned int nz);
std::vector<unsigned int> StructureEdging(std::vector<unsigned int>& structure, unsigned int nx, unsigned int ny, unsigned int nz);
std::vector<unsigned int> StructureEdgingManual(std::vector<unsigned int>& structure, unsigned int extra_voxels_x, unsigned int extra_voxels_y, unsigned int extra_voxels_z, unsigned int nx, unsigned int ny, unsigned int nz);
void CellGridExpansion(std::vector<std::array<double, 3>>& seeds_grid_sample, std::array<double, 3>& sample_length, std::array<int, 3>& sample_size);
void CellGridInitialization(std::vector<std::array<double, 3>>& seeds_grid_sample, std::array<double, 3>& sample_length, std::array<int, 3>& sample_size, double max_length, int max_dimension, const std::string& grid_type, const std::string& distance_type);
std::vector<double> GridComputation(double length, int size);
std::vector<double> GridCentersComputation(std::vector<double>& x);
std::vector<std::array<double, 3>> PrintGridToFile(double max_length, int max_dimension, const std::string& grid_type);