#pragma once



std::list<std::vector<unsigned int>> SolutionPostProcessing(std::list<std::vector<long int>>& unprocessed_data);
void CheckContiniouty(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, int current_flat_index, int size_x, int size_y, int size_z, long int num_iteration);
int TestBlockNum(int x_block, int y_block, int z_block);
void TestIndexation(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block);
bool OpenCellCheck(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z);
void OpenCellFunction(int flat_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, int radius, int size_x, int size_y, int size_z);
void OpenStructCopmutation(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, int radius, int size_x, int size_y, int size_z);
int CloseStructCheckInternal(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id);
int CloseStructCheckBoundary(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id);
int DetectingNeighbourhood(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id);
std::priority_queue<std::array<double, 9>> ListInitialization(std::vector<std::array<double, 3>>& seeds_centers, std::vector<double>& centers_x, std::vector<double>& centers_y, std::vector<double>& centers_z, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double phi, double theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry);
long int GetCellIdentifier(std::array<int, 3>& block_identifier, int seed_identifier_1, int seeds_number);
std::list<std::vector<long int>> MethodGrowth(std::vector<std::array<double, 3>>& seeds, double phi, double theta, double length_x, double length_y, double length_z, int size_x, int size_y, int size_z, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry);
void InitializeGrowth(double phi, double theta, double max_length, int max_dimension, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, const std::string& seed_distribution, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry);
void InitializeRandomExploration(double max_length, int max_dimension, int edge_radius, unsigned int ensemble_size, const std::string& seed_distribution, const std::string& distance_type);