#pragma once



std::vector<std::array<double, 3>> vertices_initialization(std::string polyhedron_family, std::string polyhedron_type);
double SquaredEucledianDistance(double x0, double y0, double z0, double x1, double y1, double z1);
double EucledianDistance(double x0, double y0, double z0, double x1, double y1, double z1);
double PeriodicDistanceFunction(double theta, double phi, double d_theta);
double PeriodicDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double Distance(double x1, double y1, double z1, double x2, double y2, double z2, double phi, double theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry);
void DistancePrintToFile(double phi, double theta, double max_length, int max_dimension, unsigned int num_phi, unsigned int num_theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& grid_type, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry, const std::string& save_folder, const std::string& name);
void PrintMultiDistInterpolationToFile(unsigned int num_phi, unsigned int num_theta, unsigned int discretization_number, const std::string& save_folder, const std::string& name);
void PrintPolyhedron(unsigned int discretization_number, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry);