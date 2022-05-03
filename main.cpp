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

const int edge_radius = 2;

const float pi = 3.141592653589793238462643;

int DetectingNeighbourhood(std::vector<long int> data, std::vector<double> set_x, std::vector<double> set_y, std::vector<double> set_z, std::array<int, 3> current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id)
{
    int result = 0;
    int index = 0;
    int flat_current_index = (size_z - 1)*(size_y - 1)*current_index[0] + (size_z - 1)*current_index[1] + current_index[2];
    int index_x, index_y, index_z;
    long int temp_cell_id = 0;
    bool change_block = false;

    if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
    {
        for (int i = - radius; i <= radius; i++)
        {
            for (int j =  - radius; j <= radius; j++)
            {
                for (int k = - radius; k <= radius; k++)
                {
                    index = (size_z - 1)*(size_y - 1)*(current_index[0] + i) + (size_z - 1)*(current_index[1] + j) + (current_index[2] + k);
                    if ((data[index] != - 1) && (data[index] != - 2) && (data[index] != - 3))
                    {
                        if (data[index] != current_cell_id)
                        {
                            result = 1;
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (int i = - radius; i <= radius; i++)
        {
            if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
            {
                index_x = abs(size_x - 1 - abs(current_index[0] + i));
                change_block = true;
            }
            else
            {
                index_x = current_index[0] + i;
            }
            for (int j =  - radius; j <= radius; j++)
            {
                if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
                {
                    index_y = abs(size_y - 1 - abs(current_index[1] + j));
                    change_block = true;
                }
                else
                {
                    index_y = current_index[1] + j;
                }
                for (int k = - radius; k <= radius; k++)
                {
                    if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
                    {
                        index_z = abs(size_z - 1 - abs(current_index[2] + k));
                        change_block = true;
                    }
                    else
                    {
                        index_z = current_index[2] + k;
                    }
                    index = (size_z - 1)*(size_y - 1)*index_x + (size_z - 1)*index_y + index_z;
                    if ((data[index] != - 1) && (data[index] != - 2) && (data[index] != - 3))
                    {
                        if (change_block == true)
                        {
                            result = 1;
                            temp_cell_id = data[index];
                            change_block = false;
                        }
                        else
                        {
                            if (data[index] != data[flat_current_index])
                            {
                                result = 1;
                                temp_cell_id = data[index];
                                change_block = false;
                            }
                        }
                    }
                }
            }
        }
    }

    return(result);
}

void PrintToFile(std::list<std::vector<long int>> cells_distribution, std::string name)
{
    std::vector<long int> solution_1 = cells_distribution.front();
    std::vector<long int> solution_2 = cells_distribution.back();

    std::string header_1("voronoy_with_cells_");
    std::string header_2("voronoy_with_seeds_");
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

    std:: cout << resulting_name_1 << " " << resulting_name_2 << std::endl;
}

double Distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double x_dist = x2 - x1;
    double y_dist = y2 - y1;
    double z_dist = z2 - z1;
    double phi = 0.5*pi + atan2(y_dist, x_dist);
    double theta = 0.0;
    if (z_dist != 0.0)
    {
        theta = atan(sqrt(x_dist*x_dist + y_dist*y_dist) / z_dist);
    }
    else
    {
        theta = 0.5*pi;
    }

    //return (0.25*sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist) / (4.0 + sin(6.0*theta)*cos(2.0*phi)));
    //return (sqrt(cos(0.5*pi)*cos(0.5*pi)*x_dist*x_dist / 10.0 + sin(0.5*pi)*sin(0.5*pi)*y_dist*y_dist / 10.0 + z_dist*z_dist));
    return (0.1*sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist) / (0.3 + 0.7*fabs(sin(1.0*theta)*sin(1.0*phi))));
}

std::vector<double> GridComputation(double length, int size)
{
    std::vector<double> grid;

    double step = length / (size - 1);

    for (int i = 0; i < size; i++)
    {
        grid.push_back(i*step);
    }

    return grid;
}

std::vector<double> GridCentersComputation(std::vector<double> x)
{
    std::vector<double> centers;

    for (int i = 0; i < x.size() - 1; i++)
    {
        centers.push_back((x[i + 1] + x[i]) / 2.0);
    }

    return centers;
}

std::priority_queue<std::array<double, 9>> ListInitialization(std::vector<std::array<double, 3>> seeds_centers, std::vector<double> centers_x, std::vector<double> centers_y, std::vector<double> centers_z, std::vector<double> x, std::vector<double> y, std::vector<double> z)
{
    int seeds_number = seeds_centers.size();
    int x_number = x.size();
    int y_number = y.size();
    int z_number = z.size();
    std::priority_queue<std::array<double,9>> Queue;
    std::array<double, 9> temp_for_push = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<int, 3> corresponding_center = {0, 0, 0};
    int seed_id = - 1;
    int voxel_counter = 0;
    bool error = false;
    double dist = 0.0;
    //std::cout << Queue.size() << std::endl;

    if (seeds_number == 0)
    {
        std::cout << "ERROR: missing seeds for growing." << std::endl;
        error = true;
    }

    if (error != true)
    {
        for (auto seed : seeds_centers)
        {
            seed_id = seed_id + 1;
            for (int i = 0; i < x_number - 1; i++)
            {
                if ((seed[0] >= x[i]) && (seed[0] <= x[i + 1]))
                {
                    corresponding_center[0] = i;
                    for (int j = 0; j < y_number - 1; j++)
                    {
                        if ((seed[1] >= y[j]) && (seed[1] <= y[j + 1]))
                        {
                            corresponding_center[1] = j;
                            for (int k = 0; k < z_number - 1; k++)
                            {
                                if ((seed[2] >= z[k]) && (seed[2] <= z[k + 1]))
                                {
                                    corresponding_center[2] = k;
                                }
                            }
                        }
                    }
                }
            }
            dist = Distance(seed[0], seed[1], seed[2], centers_x[corresponding_center[0]], centers_y[corresponding_center[1]], centers_z[corresponding_center[2]]);
            temp_for_push = {- dist, (double)voxel_counter, (double)corresponding_center[0], (double)corresponding_center[1], (double)corresponding_center[2], (double)seed_id, 0.0, 0.0, 0.0};
            Queue.push(temp_for_push);
            //std::cout << Queue.size() << std::endl;
            voxel_counter = voxel_counter + 1;
        }
    }

    return Queue;
}

long int GetCellIdentifier(std::array<int, 3> block_identifier, int seed_identifier_1, int seeds_number)
{
    long int cell_identifier = 0;
    int max_block = 15;

    if ((block_identifier[0] == 0) && (block_identifier[1] == 0) && (block_identifier[2] == 0))
    {
        cell_identifier = seed_identifier_1;
    }
    else
    {
        double seed_identifier = seed_identifier_1 + 1;
        if ((block_identifier[0] >= 0) && (block_identifier[1] >= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = (max_block*(max_block - 1)*block_identifier[0] + (max_block - 1)*block_identifier[1] + block_identifier[2] + seeds_number)*seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] <= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block*(max_block - 1)*block_identifier[0] + (max_block - 1)*block_identifier[1] + block_identifier[2] - seeds_number)*seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] >= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] <= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = - (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] <= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seeds_number*seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] >= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = - (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seeds_number*seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] >= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seeds_number*seeds_number*seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] <= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = - (max_block*(max_block - 1)*abs(block_identifier[0]) + (max_block - 1)*abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number)*seeds_number*seeds_number*seeds_number*seed_identifier;
        }
    }

    return cell_identifier;
}

std::list<std::vector<long int>> MethodGrowth(std::vector<std::array<double, 3>> seeds, double length_x, double length_y, double length_z, int size_x, int size_y, int size_z)
{
    std::vector<double> grid_x = GridComputation(length_x, size_x);
    std::vector<double> grid_y = GridComputation(length_y, size_y);
    std::vector<double> grid_z = GridComputation(length_z, size_z);
    std::vector<double> grid_centers_x = GridCentersComputation(grid_x);
    std::vector<double> grid_centers_y = GridCentersComputation(grid_y);
    std::vector<double> grid_centers_z = GridCentersComputation(grid_z);
    std::vector<long int> result_full, result_seeds;
    std::list<std::vector<long int>> final, temporary;

    for (int i = 0; i < grid_centers_x.size(); i++)
    {
        for (int j = 0; j < grid_centers_y.size(); j++)
        {
            for (int k = 0; k < grid_centers_z.size(); k++)
            {
                result_full.push_back(-1);
                result_seeds.push_back(-1);
            }
        }
    }

    std::priority_queue<std::array<double, 9>> Queue = ListInitialization(seeds, grid_centers_x, grid_centers_y, grid_centers_z, grid_x, grid_y, grid_z);

    std::array<int, 3> left_neighbour_block = {0, 0, 0};
    std::array<int, 3> right_neighbour_block = {0, 0, 0};
    std::array<int, 3> front_neighbour_block = {0, 0, 0};
    std::array<int, 3> back_neighbour_block = {0, 0, 0};
    std::array<int, 3> up_neighbour_block = {0, 0, 0};
    std::array<int, 3> down_neighbour_block = {0, 0, 0};

    std::array<int, 3> left_neighbour_center = {0, 0, 0};
    std::array<int, 3> right_neighbour_center = {0, 0, 0};
    std::array<int, 3> front_neighbour_center = {0, 0, 0};
    std::array<int, 3> back_neighbour_center = {0, 0, 0};
    std::array<int, 3> up_neighbour_center = {0, 0, 0};
    std::array<int, 3> down_neighbour_center = {0, 0, 0};

    std::array<double, 9> current_info, temp_to_push;
    std::array<double, 3> temp_point, seed = {0.0, 0.0, 0.0};
    std::array<int, 3> current_cell;
    std::array<int, 3> current_block;
    double dist = 0.0;
    int current_seed_id, current_index, index, neighbourhood, current_neighbourhood, iterator = 0;
    long int voxel_id = Queue.size();
    int initial_id = voxel_id;

    double x_last = grid_x.back();
    double y_last = grid_y.back();
    double z_last = grid_z.back();

    while(Queue.size() != 0)
    {
        current_info = Queue.top();
        Queue.pop();
        current_cell[0] = (int)current_info[2];
        current_cell[1] = (int)current_info[3];
        current_cell[2] = (int)current_info[4];
        current_index = (size_z - 1)*(size_y - 1)*current_cell[0] + (size_z - 1)*current_cell[1] + current_cell[2]; //look here again
        current_seed_id = (int)current_info[5];
        current_block[0] = (int)current_info[6];
        current_block[1] = (int)current_info[7];
        current_block[2] = (int)current_info[8];
        seed = seeds[current_seed_id];
        //std::cout << "result = " << current_index << std::endl;

        //current_neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, current_cell, edge_radius, size_x, size_y, size_z, result_full[current_index]);
        if (result_full[current_index] == - 1)
        {
            result_full[current_index] = GetCellIdentifier(current_block, current_seed_id, seeds.size());
            result_seeds[current_index] = current_seed_id;

            left_neighbour_block = current_block;
            left_neighbour_center[1] = current_cell[1];
            left_neighbour_center[2] = current_cell[2];
            if (current_cell[0] == 0)
            {
                left_neighbour_center[0] = size_x - 2;
                left_neighbour_block[0] = left_neighbour_block[0] - 1;
            }
            else
            {
                left_neighbour_center[0] = current_cell[0] - 1;
            }
            index = (size_z - 1)*(size_y - 1)*left_neighbour_center[0] + (size_z - 1)*left_neighbour_center[1] + left_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, left_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[left_neighbour_center[0]] + left_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[left_neighbour_center[1]] + left_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[left_neighbour_center[2]] + left_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)left_neighbour_center[0], (double)left_neighbour_center[1], (double)left_neighbour_center[2], (double)current_seed_id, (double)left_neighbour_block[0], (double)left_neighbour_block[1], (double)left_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }

            right_neighbour_block = current_block;
            right_neighbour_center[1] = current_cell[1];
            right_neighbour_center[2] = current_cell[2];
            if (current_cell[0] == (size_x - 2))
            {
                right_neighbour_center[0] = 0;
                right_neighbour_block[0] = right_neighbour_block[0] + 1;
            }
            else
            {
                right_neighbour_center[0] = current_cell[0] + 1;
            }
            index = (size_z - 1)*(size_y - 1)*right_neighbour_center[0] + (size_z - 1)*right_neighbour_center[1] + right_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, right_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[right_neighbour_center[0]] + right_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[right_neighbour_center[1]] + right_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[right_neighbour_center[2]] + right_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)right_neighbour_center[0], (double)right_neighbour_center[1], (double)right_neighbour_center[2], (double)current_seed_id, (double)right_neighbour_block[0], (double)right_neighbour_block[1], (double)right_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }

            front_neighbour_block = current_block;
            front_neighbour_center[0] = current_cell[0];
            front_neighbour_center[2] = current_cell[2];
            if (current_cell[1] == 0)
            {
                front_neighbour_center[1] = size_y - 2;
                front_neighbour_block[1] = front_neighbour_block[1] - 1;
            }
            else
            {
                front_neighbour_center[1] = current_cell[1] - 1;
            }
            index = (size_z - 1)*(size_y - 1)*front_neighbour_center[0] + (size_z - 1)*front_neighbour_center[1] + front_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, front_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[front_neighbour_center[0]] + front_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[front_neighbour_center[1]] + front_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[front_neighbour_center[2]] + front_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)front_neighbour_center[0], (double)front_neighbour_center[1], (double)front_neighbour_center[2], (double)current_seed_id, (double)front_neighbour_block[0], (double)front_neighbour_block[1], (double)front_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }

            back_neighbour_block = current_block;
            back_neighbour_center[0] = current_cell[0];
            back_neighbour_center[2] = current_cell[2];
            if (current_cell[1] == (size_y - 2))
            {
                back_neighbour_center[1] = 0;
                back_neighbour_block[1] = back_neighbour_block[1] + 1;
            }
            else
            {
                back_neighbour_center[1] = current_cell[1] + 1;
            }
            index = (size_z - 1)*(size_y - 1)*back_neighbour_center[0] + (size_z - 1)*back_neighbour_center[1] + back_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, back_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[back_neighbour_center[0]] + back_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[back_neighbour_center[1]] + back_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[back_neighbour_center[2]] + back_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)back_neighbour_center[0], (double)back_neighbour_center[1], (double)back_neighbour_center[2], (double)current_seed_id, (double)back_neighbour_block[0], (double)back_neighbour_block[1], (double)back_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }

            up_neighbour_block = current_block;
            up_neighbour_center[0] = current_cell[0];
            up_neighbour_center[1] = current_cell[1];
            if (current_cell[2] == (size_z - 2))
            {
                up_neighbour_center[2] = 0;
                up_neighbour_block[2] = up_neighbour_block[2] + 1;
            }
            else
            {
                up_neighbour_center[2] = current_cell[2] + 1;
            }
            index = (size_z - 1)*(size_y - 1)*up_neighbour_center[0] + (size_z - 1)*up_neighbour_center[1] + up_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, up_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[up_neighbour_center[0]] + up_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[up_neighbour_center[1]] + up_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[up_neighbour_center[2]] + up_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)up_neighbour_center[0], (double)up_neighbour_center[1], (double)up_neighbour_center[2], (double)current_seed_id, (double)up_neighbour_block[0], (double)up_neighbour_block[1], (double)up_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }

            down_neighbour_block = current_block;
            down_neighbour_center[0] = current_cell[0];
            down_neighbour_center[1] = current_cell[1];
            if (current_cell[2] == 0)
            {
                down_neighbour_center[2] = size_z - 2;
                down_neighbour_block[2] = down_neighbour_block[2] - 1;
            }
            else
            {
                down_neighbour_center[2] = current_cell[2] - 1;
            }
            index = (size_z - 1)*(size_y - 1)*down_neighbour_center[0] + (size_z - 1)*down_neighbour_center[1] + down_neighbour_center[2];
            neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, down_neighbour_center, edge_radius, size_x, size_y, size_z, result_full[current_index]);
            if ((result_full[index] == - 1) && (neighbourhood == 0))
            {
                temp_point[0] = grid_centers_x[down_neighbour_center[0]] + down_neighbour_block[0]*x_last;
                temp_point[1] = grid_centers_y[down_neighbour_center[1]] + down_neighbour_block[1]*y_last;
                temp_point[2] = grid_centers_z[down_neighbour_center[2]] + down_neighbour_block[2]*z_last;
                dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2]);
                temp_to_push = {- dist, (double)voxel_id, (double)down_neighbour_center[0], (double)down_neighbour_center[1], (double)down_neighbour_center[2], (double)current_seed_id, (double)down_neighbour_block[0], (double)down_neighbour_block[1], (double)down_neighbour_block[2]};
                Queue.push(temp_to_push);
                voxel_id = voxel_id + 1;
            }
            else
            {
                if (result_full[index] == - 1)
                {
                    if (neighbourhood == 1)
                    {
                        result_full[index] = - 2;
                    }
                    else
                    {
                        result_full[index] = - 3;
                    }
                }
            }
            /*
            iterator = iterator + 1;
            if ((iterator >= 440320) && (iterator <= 460800))
            {
                if ((iterator % int((460800 - 409600) / 20.0)) == 0)
                {
                    temporary.push_back(result_full);
                    temporary.push_back(result_seeds);
                    std::string file_name;
                    std::stringstream int_to_string;
                    int_to_string << iterator;
                    int_to_string >> file_name;
                    PrintToFile(temporary, file_name);
                    temporary.clear();
                }
            }*/
            //std::cout << Queue.size() << std::endl;
        }
    }

    final.push_back(result_full);
    final.push_back(result_seeds);

    return final;
}

int main()
{
    std::vector<std::array<double, 3>> seeds;
    std::array<double, 3> temp = {0.25, 0.75, 0.25};
    seeds.push_back(temp);
    /*temp = {0.25, 0.25, 0.75};
    seeds.push_back(temp);
    temp = {0.75, 0.25, 0.25};
    seeds.push_back(temp);
    temp = {0.25, 0.75, 0.25};
    seeds.push_back(temp);
    temp = {0.25, 0.75, 0.75};
    seeds.push_back(temp);
    temp = {0.75, 0.25, 0.75};
    seeds.push_back(temp);
    temp = {0.75, 0.75, 0.25};
    seeds.push_back(temp);
    temp = {0.75, 0.75, 0.75};
    seeds.push_back(temp);*/
    //std::cout << seeds.size() << std::endl;

    clock_t starting_point, ending_point;
    starting_point = clock();
    std::list<std::vector<long int>> k = MethodGrowth(seeds, 1.0, 1.0, 1.0, 31, 31, 31);
    ending_point = clock();

    std::cout << "Computational time = " << (double)(ending_point - starting_point) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;

    //std::string file_id("lalochka");
    //PrintToFile(k, file_id);

    return 0;
}
