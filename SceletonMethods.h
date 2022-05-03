#pragma once
#include <vector>

namespace PostProcessing
{
	std::array<int, 3> IndexFlatToTriple(int flat_index, int num_x, int num_y, int num_z);
	std::vector<unsigned int> StructureIntoVolume(std::list<std::vector<unsigned int>>& structure, std::array<double, 3>& size, std::array<int, 3>& dimensions, unsigned int thickness, const std::string& surface_type);
	void SkeletonWriteToFile(std::vector<unsigned int>& identificator, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers);
}