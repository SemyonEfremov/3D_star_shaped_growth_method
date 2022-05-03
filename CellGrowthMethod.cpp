#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <queue>
#include <list>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include "parameters.h"
#include "GridComputationMethods.h"
#include "DistanceModel.h"
#include "MetaDataWrite.h"
#include "SceletonMethods.h"





std::list<std::vector<unsigned int>> SolutionPostProcessing(std::list<std::vector<long int>>& unprocessed_data)
{
    std::cout << "Starting post processing of the solution..." << std::endl << std::endl;

    std::list<std::vector<unsigned int>> result;
    std::vector<unsigned int> proc_closed_struct, proc_open_struct;
    std::vector<long int> unproc_closed_struct, unproc_open_struct;
    
    unproc_closed_struct = unprocessed_data.front();
    unproc_open_struct = unprocessed_data.back();

    //#pragma omp parallel for
    for (int index = 0; index < unproc_closed_struct.size(); index++)
    {
        if ((unproc_closed_struct[index] < 0) && (unproc_closed_struct[index] >= -3))
        {
            proc_closed_struct.push_back(1);
        }
        else
        {
            proc_closed_struct.push_back(0);
        }
        if (unproc_open_struct[index] == -3)
        {
            proc_open_struct.push_back(1);
        }
        else
        {
            proc_open_struct.push_back(0);
        }
    }

    result.push_back(proc_closed_struct);
    result.push_back(proc_open_struct);

    return(result);
}

void CheckContiniouty(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, int current_flat_index, int size_x, int size_y, int size_z, long int num_iteration)
{
    std::array<int, 3> current_index;
    current_index[0] = (int)(floor(floor(current_flat_index / (size_z - 1)) / (size_y - 1)));
    current_index[1] = (int)((int)(floor(current_flat_index / (size_z - 1))) % (size_y - 1));
    current_index[2] = (int)(current_flat_index % (size_z - 1));
    int index_1, index_2, block_x_1 = 0, block_y_1 = 0, block_z_1 = 0, block_x_2 = 0, block_y_2 = 0, block_z_2 = 0;
    bool signalizer = false;

    if (current_index[2] == 0)
    {
        index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] + 1);
        index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + size_z - 2;
        block_z_1 = -1;
    }
    else
    {
        if (current_index[2] == size_z - 2)
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + 0;
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] - 1);
            block_z_2 = 1;
        }
        else
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] + 1);
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] - 1);
        }
    }
    if ((data[current_flat_index] == data[index_1]) && (x_block[current_flat_index] + block_x_2 == x_block[index_1]) && (y_block[current_flat_index] + block_y_2 == y_block[index_1]) && (z_block[current_flat_index] + block_z_2 == z_block[index_1]))
    {
        signalizer = true;
    }
    if ((data[current_flat_index] == data[index_2]) && (x_block[current_flat_index] + block_x_1 == x_block[index_2]) && (y_block[current_flat_index] + block_y_1 == y_block[index_2]) && (z_block[current_flat_index] + block_z_1 == z_block[index_2]))
    {
        signalizer = true;
    }

    block_x_1 = 0;
    block_y_1 = 0;
    block_z_1 = 0;
    block_x_2 = 0;
    block_y_2 = 0;
    block_z_2 = 0;

    if (current_index[1] == 0)
    {
        index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] + 1) + current_index[2];
        index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (size_y - 2) + current_index[2];
        block_y_1 = -1;
    }
    else
    {
        if (current_index[1] == size_y - 2)
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * 0 + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] - 1) + current_index[2];
            block_y_2 = 1;
        }
        else
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] + 1) + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] - 1) + current_index[2];
        }
    }
    if ((data[current_flat_index] == data[index_1]) && (x_block[current_flat_index] + block_x_2 == x_block[index_1]) && (y_block[current_flat_index] + block_y_2 == y_block[index_1]) && (z_block[current_flat_index] + block_z_2 == z_block[index_1]))
    {
        signalizer = true;
    }
    if ((data[current_flat_index] == data[index_2]) && (x_block[current_flat_index] + block_x_1 == x_block[index_2]) && (y_block[current_flat_index] + block_y_1 == y_block[index_2]) && (z_block[current_flat_index] + block_z_1 == z_block[index_2]))
    {
        signalizer = true;
    }

    block_x_1 = 0;
    block_y_1 = 0;
    block_z_1 = 0;
    block_x_2 = 0;
    block_y_2 = 0;
    block_z_2 = 0;

    if (current_index[0] == 0)
    {
        index_1 = (size_z - 1) * (size_y - 1) * (current_index[0] + 1) + (size_z - 1) * current_index[1] + current_index[2];
        index_2 = (size_z - 1) * (size_y - 1) * (size_x - 2) + (size_z - 1) * current_index[1] + current_index[2];
        block_x_1 = -1;
    }
    else
    {
        if (current_index[0] == size_x - 2)
        {
            index_1 = (size_z - 1) * (size_y - 1) * 0 + (size_z - 1) * current_index[1] + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * (current_index[0] - 1) + (size_z - 1) * current_index[1] + current_index[2];
            block_x_2 = 1;
        }
        else
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
        }
    }
    if ((data[current_flat_index] == data[index_1]) && (x_block[current_flat_index] + block_x_2 == x_block[index_1]) && (y_block[current_flat_index] + block_y_2 == y_block[index_1]) && (z_block[current_flat_index] + block_z_2 == z_block[index_1]))
    {
        signalizer = true;
    }
    if ((data[current_flat_index] == data[index_2]) && (x_block[current_flat_index] + block_x_1 == x_block[index_2]) && (y_block[current_flat_index] + block_y_1 == y_block[index_2]) && (z_block[current_flat_index] + block_z_1 == z_block[index_2]))
    {
        signalizer = true;
    }

    block_x_1 = 0;
    block_y_1 = 0;
    block_z_1 = 0;
    block_x_2 = 0;
    block_y_2 = 0;
    block_z_2 = 0;

    if (signalizer == false)
    {
        std::cout << "At " << num_iteration << "-th iteration a disconnected part detected with the index (" << current_index[0] << ", " << current_index[1] << ", " << current_index[2] << "), <" << data[current_flat_index] << ", " << x_block[current_flat_index] << ", " << y_block[current_flat_index] << ", " << z_block[current_flat_index] << ">." << std::endl;
        if (current_index[2] == 0)
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] + 1);
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + size_z - 2;
            block_z_1 = -1;
        }
        else
        {
            if (current_index[2] == size_z - 2)
            {
                index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + 0;
                index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] - 1);
                block_z_2 = 1;
            }
            else
            {
                index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] + 1);
                index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + (current_index[2] - 1);
            }
        }
        std::cout << "z-neighbor <" << data[index_1] << ", " << x_block[index_1] << ", " << y_block[index_1] << ", " << z_block[index_1] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_2 << ", " << y_block[current_flat_index] + block_y_2 << ", " << z_block[current_flat_index] + block_z_2 << ">, (" << (int)(floor(floor(index_1 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_1 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_1 % (size_z - 1)) << ")" << std::endl;

        std::cout << "z-neighbor <" << data[index_2] << ", " << x_block[index_2] << ", " << y_block[index_2] << ", " << z_block[index_2] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_1 << ", " << y_block[current_flat_index] + block_y_1 << ", " << z_block[current_flat_index] + block_z_1 << ">, (" << (int)(floor(floor(index_2 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_2 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_2 % (size_z - 1)) << ")" << std::endl;

        block_x_1 = 0;
        block_y_1 = 0;
        block_z_1 = 0;
        block_x_2 = 0;
        block_y_2 = 0;
        block_z_2 = 0;

        if (current_index[1] == 0)
        {
            index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] + 1) + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (size_y - 2) + current_index[2];
            block_y_1 = -1;
        }
        else
        {
            if (current_index[1] == size_y - 2)
            {
                index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * 0 + current_index[2];
                index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] - 1) + current_index[2];
                block_y_2 = 1;
            }
            else
            {
                index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] + 1) + current_index[2];
                index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * (current_index[1] - 1) + current_index[2];
            }
        }
        std::cout << "y-neighbor <" << data[index_1] << ", " << x_block[index_1] << ", " << y_block[index_1] << ", " << z_block[index_1] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_2 << ", " << y_block[current_flat_index] + block_y_2 << ", " << z_block[current_flat_index] + block_z_2 << ">, (" << (int)(floor(floor(index_1 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_1 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_1 % (size_z - 1)) << ")" << std::endl;

        std::cout << "y-neighbor <" << data[index_2] << ", " << x_block[index_2] << ", " << y_block[index_2] << ", " << z_block[index_2] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_1 << ", " << y_block[current_flat_index] + block_y_1 << ", " << z_block[current_flat_index] + block_z_1 << ">, (" << (int)(floor(floor(index_2 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_2 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_2 % (size_z - 1)) << ")" << std::endl;

        block_x_1 = 0;
        block_y_1 = 0;
        block_z_1 = 0;
        block_x_2 = 0;
        block_y_2 = 0;
        block_z_2 = 0;

        if (current_index[0] == 0)
        {
            index_1 = (size_z - 1) * (size_y - 1) * (current_index[0] + 1) + (size_z - 1) * current_index[1] + current_index[2];
            index_2 = (size_z - 1) * (size_y - 1) * (size_x - 2) + (size_z - 1) * current_index[1] + current_index[2];
            block_x_1 = -1;
        }
        else
        {
            if (current_index[0] == size_x - 2)
            {
                index_1 = (size_z - 1) * (size_y - 1) * 0 + (size_z - 1) * current_index[1] + current_index[2];
                index_2 = (size_z - 1) * (size_y - 1) * (current_index[0] - 1) + (size_z - 1) * current_index[1] + current_index[2];
                block_x_2 = 1;
            }
            else
            {
                index_1 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
                index_2 = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
            }
        }
        std::cout << "x-neighbor <" << data[index_1] << ", " << x_block[index_1] << ", " << y_block[index_1] << ", " << z_block[index_1] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_2 << ", " << y_block[current_flat_index] + block_y_2 << ", " << z_block[current_flat_index] + block_z_2 << ">, (" << (int)(floor(floor(index_1 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_1 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_1 % (size_z - 1)) << ")" << std::endl;

        std::cout << "x-neighbor <" << data[index_2] << ", " << x_block[index_2] << ", " << y_block[index_2] << ", " << z_block[index_2] << ">, <" << data[current_flat_index] << ", " << x_block[current_flat_index] + block_x_1 << ", " << y_block[current_flat_index] + block_y_1 << ", " << z_block[current_flat_index] + block_z_1 << ">, (" << (int)(floor(floor(index_2 / (size_z - 1)) / (size_y - 1))) << ", " << (int)((int)(floor(index_2 / (size_z - 1))) % (size_y - 1)) << ", " << (int)(index_2 % (size_z - 1)) << ")" << std::endl;

        block_x_1 = 0;
        block_y_1 = 0;
        block_z_1 = 0;
        block_x_2 = 0;
        block_y_2 = 0;
        block_z_2 = 0;
    }
}

int TestBlockNum(int x_block, int y_block, int z_block)
{
    int result = 0;
    if (x_block == -1)
    {
        if (y_block == -1)
        {
            if (z_block == -1)
            {
                result = 10;
            }
            else
            {
                if (z_block == 1)
                {
                    result = 20;
                }
                else
                {
                    result = 30;
                }
            }
        }
        else
        {
            if (y_block == 1)
            {
                if (z_block == -1)
                {
                    result = 40;
                }
                else
                {
                    if (z_block == 1)
                    {
                        result = 50;
                    }
                    else
                    {
                        result = 60;
                    }
                }
            }
            else
            {
                if (z_block == -1)
                {
                    result = 70;
                }
                else
                {
                    if (z_block == 1)
                    {
                        result = 80;
                    }
                    else
                    {
                        result = 90;
                    }
                }
            }
        }
    }
    else
    {
        if (x_block == 1)
        {
            if (y_block == -1)
            {
                if (z_block == -1)
                {
                    result = 100;
                }
                else
                {
                    if (z_block == 1)
                    {
                        result = 110;
                    }
                    else
                    {
                        result = 120;
                    }
                }
            }
            else
            {
                if (y_block == 1)
                {
                    if (z_block == -1)
                    {
                        result = 130;
                    }
                    else
                    {
                        if (z_block == 1)
                        {
                            result = 140;
                        }
                        else
                        {
                            result = 150;
                        }
                    }
                }
                else
                {
                    if (z_block == -1)
                    {
                        result = 160;
                    }
                    else
                    {
                        if (z_block == 1)
                        {
                            result = 170;
                        }
                        else
                        {
                            result = 180;
                        }
                    }
                }
            }
        }
        else
        {
            if (y_block == -1)
            {
                if (z_block == -1)
                {
                    result = 190;
                }
                else
                {
                    if (z_block == 1)
                    {
                        result = 200;
                    }
                    else
                    {
                        result = 210;
                    }
                }
            }
            else
            {
                if (y_block == 1)
                {
                    if (z_block == -1)
                    {
                        result = 220;
                    }
                    else
                    {
                        if (z_block == 1)
                        {
                            result = 230;
                        }
                        else
                        {
                            result = 240;
                        }
                    }
                }
                else
                {
                    if (z_block == -1)
                    {
                        result = 250;
                    }
                    else
                    {
                        if (z_block == 1)
                        {
                            result = 260;
                        }
                        else
                        {
                            result = 0;
                        }
                    }
                }
            }
        }
    }
    return (result);
}

void TestIndexation(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block)
{
    int temp_ind = 0;
    for (int i = 0; i < data.size(); i++)
    {
        temp_ind = TestBlockNum(x_block[i], y_block[i], z_block[i]);
        data[i] = data[i] + temp_ind;
    }
}

bool OpenCellCheck(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z)
{
    bool result = false, not_to_push = false;
    int index = 0, i, j, k;
    int flat_current_index = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
    int index_x, index_y, index_z, cur_x_block, cur_y_block, cur_z_block;
    long int temp_cell_id = -1, temp_cell_id_nochange = -1;
    std::array<bool, 3> change_block = { false, false, false };
    int change_block_x = 0, change_block_y = 0, change_block_z = 0, temp_change_block_x = 0, temp_change_block_y = 0, temp_change_block_z = 0, temp_cell_block_x = 0, temp_cell_block_y = 0, temp_cell_block_z = 0;
    int temp_cell_block_nochange_x = 0, temp_cell_block_nochange_y = 0, temp_cell_block_nochange_z = 0;

    cur_x_block = x_block[flat_current_index];
    cur_y_block = y_block[flat_current_index];
    cur_z_block = z_block[flat_current_index];

    std::vector<std::array<int, 4>> surrounding_neighborhood;
    std::array<int, 4> temp_compare, current_data_to_compare;
    std::vector<std::array<int, 3>> surrounding_index;
    std::array<int, 3> temp_index;

    if (surrounding_neighborhood.size() != 0)
    {
        std::cout << "ERROR, the vector IS NOT empty!" << std::endl;
    }

    if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
    {
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
            j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
            k = (int)(par_ind % (2 * radius + 1)) - radius;

            index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
            current_data_to_compare[0] = (int)data[index];
            current_data_to_compare[1] = x_block[index];
            current_data_to_compare[2] = y_block[index];
            current_data_to_compare[3] = z_block[index];

            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3))
            {
                for (int l = 0; l < surrounding_neighborhood.size(); l++)
                {
                    temp_compare = surrounding_neighborhood[l];
                    if ((current_data_to_compare[0] == temp_compare[0]) && (current_data_to_compare[1] == temp_compare[1]) && (current_data_to_compare[2] == temp_compare[2]) && (current_data_to_compare[3] == temp_compare[3]))
                    {
                        not_to_push = true;
                    }
                }
                if (not_to_push == false)
                {
                    surrounding_neighborhood.push_back(current_data_to_compare);
                }
                not_to_push = false;
            }
        }
    }
    else
    {
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
            if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
            {
                index_x = abs(size_x - 1 - abs(current_index[0] + i));
                change_block_x = (int)copysign(1.0, current_index[0] + i);
                change_block[0] = true;
            }
            else
            {
                index_x = current_index[0] + i;
            }
            j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
            if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
            {
                index_y = abs(size_y - 1 - abs(current_index[1] + j));
                change_block_y = (int)copysign(1.0, current_index[1] + j);
                change_block[1] = true;
            }
            else
            {
                index_y = current_index[1] + j;
            }
            k = (int)(par_ind % (2 * radius + 1)) - radius;
            if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
            {
                index_z = abs(size_z - 1 - abs(current_index[2] + k));
                change_block_z = (int)copysign(1.0, current_index[2] + k);
                change_block[2] = true;
            }
            else
            {
                index_z = current_index[2] + k;
            }

            index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
            if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3))
                {
                    current_data_to_compare[0] = (int)data[index];
                    current_data_to_compare[1] = x_block[index] - change_block_x;
                    current_data_to_compare[2] = y_block[index] - change_block_y;
                    current_data_to_compare[3] = z_block[index] - change_block_z;
                    for (int l = 0; l < surrounding_neighborhood.size(); l++)
                    {
                        temp_compare = surrounding_neighborhood[l];
                        if ((current_data_to_compare[0] == temp_compare[0]) && (current_data_to_compare[1] == temp_compare[1]) && (current_data_to_compare[2] == temp_compare[2]) && (current_data_to_compare[3] == temp_compare[3]))
                        {
                            not_to_push = true;
                        }
                    }
                    if (not_to_push == false)
                    {
                        surrounding_neighborhood.push_back(current_data_to_compare);
                    }
                    not_to_push = false;
                }
            }
            else
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3))
                {
                    current_data_to_compare[0] = (int)data[index];
                    current_data_to_compare[1] = x_block[index];
                    current_data_to_compare[2] = y_block[index];
                    current_data_to_compare[3] = z_block[index];
                    for (int l = 0; l < surrounding_neighborhood.size(); l++)
                    {
                        temp_compare = surrounding_neighborhood[l];
                        if ((current_data_to_compare[0] == temp_compare[0]) && (current_data_to_compare[1] == temp_compare[1]) && (current_data_to_compare[2] == temp_compare[2]) && (current_data_to_compare[3] == temp_compare[3]))
                        {
                            not_to_push = true;
                        }
                    }
                    if (not_to_push == false)
                    {
                        surrounding_neighborhood.push_back(current_data_to_compare);
                    }
                    not_to_push = false;
                }
            }

            change_block_z = 0;
            change_block[2] = false;
            change_block_y = 0;
            change_block[1] = false;
            change_block_x = 0;
            change_block[0] = false;
        }
    }

    if (surrounding_neighborhood.size() > 2)
    {
        result = true;
    }

    return result;
}

void OpenCellFunction(int flat_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, int radius, int size_x, int size_y, int size_z)
{
    std::array<int, 3> current_index;
    bool indicator = false;

    current_index[0] = (int)(floor(floor(flat_ind / (size_z - 1)) / (size_y - 1)));
    current_index[1] = (int)((int)(floor(flat_ind / (size_z - 1))) % (size_y - 1));
    current_index[2] = (int)(flat_ind % (size_z - 1));
    indicator = OpenCellCheck(data, x_block, y_block, z_block, current_index, radius, size_x, size_y, size_z);
    if (indicator == true)
    {
        data[flat_ind] = -3;
    }
}

std::array<int, 5> OpenCellCheckEucledianInternal(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, double minimum_neighborhood_radius)
{
    std::array<int, 5> result = { 0, 0, 0, 0, 0 };

    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    int k = (int)(par_ind % (2 * radius + 1)) - radius;

    if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[current_index[0] + i], y_centers[current_index[1] + j], z_centers[current_index[2] + k]) <= minimum_neighborhood_radius)
    {
        int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
        if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
        {
            result = { 1, data[index], x_block[index], y_block[index], z_block[index] };
        }
    }

    return(result);
}

std::array<int, 5> OpenCellCheckEucledianExternal(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, double minimum_neighborhood_radius, double length_x, double length_y, double length_z)
{
    std::array<int, 5> result = { 0, 0, 0, 0, 0 };
    int index_x = 0, index_y = 0, index_z = 0, change_block_x = 0, change_block_y = 0, change_block_z = 0;
    std::array<bool, 3> change_block = { false, false, false };
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
    {
        index_x = abs(size_x - 1 - abs(current_index[0] + i));
        change_block_x = (int)copysign(1.0, current_index[0] + i);
        change_block[0] = true;
    }
    else
    {
        index_x = current_index[0] + i;
    }
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
    {
        index_y = abs(size_y - 1 - abs(current_index[1] + j));
        change_block_y = (int)copysign(1.0, current_index[1] + j);
        change_block[1] = true;
    }
    else
    {
        index_y = current_index[1] + j;
    }
    int k = (int)(par_ind % (2 * radius + 1)) - radius;
    if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
    {
        index_z = abs(size_z - 1 - abs(current_index[2] + k));
        change_block_z = (int)copysign(1.0, current_index[2] + k);
        change_block[2] = true;
    }
    else
    {
        index_z = current_index[2] + k;
    }

    if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[index_x] + length_x * change_block_x, y_centers[index_y] + length_y * change_block_y, z_centers[index_z] + length_z * change_block_z) <= minimum_neighborhood_radius)
    {
        int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
        if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                result = { 1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z };
            }
        }
        else
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                    result = { 1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z };
            }
        }
    }

    return(result);
}

bool OpenCellCheckEucledian(int flat_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, double radius_eucledian, double length_x, double length_y, double length_z, int radius, int size_x, int size_y, int size_z, int neighborhood_elements_number)
{
    bool result = false;
    std::array<int, 3> current_index, current_block;
    current_index[0] = (int)(floor(floor(flat_ind / (size_z - 1)) / (size_y - 1)));
    current_index[1] = (int)((int)(floor(flat_ind / (size_z - 1))) % (size_y - 1));
    current_index[2] = (int)(flat_ind % (size_z - 1));
    current_block = { x_block[current_index[0]], y_block[current_index[1]], z_block[current_index[2]] };
    std::vector<std::array<int, 5>> neighbor_cells;
    std::array<int, 5> temp_cell = { 0, 0, 0, 0, 0 }, temp_cell_1 = { 0, 0, 0, 0, 0 }, temp_cell_2 = { 0, 0, 0, 0, 0 }, temp_cell_3 = { 0, 0, 0, 0, 0 };

    if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
    {
        for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
        {
            temp_cell = OpenCellCheckEucledianInternal(par_ind, data, x_block, y_block, z_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, radius_eucledian);
            neighbor_cells.push_back(temp_cell);
        }
    }
    else
    {
        for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
        {
            temp_cell = OpenCellCheckEucledianExternal(par_ind, data, x_block, y_block, z_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, radius_eucledian, length_x, length_y, length_z);
            neighbor_cells.push_back(temp_cell);

        }
    }
    for (const auto& element : neighbor_cells)
    {
        if (element[0] == 1)
        {
            if (temp_cell_1[0] == 0)
            {
                temp_cell_1 = element;
            }
            else
            {
                if (temp_cell_2[0] == 0)
                {
                    if (element != temp_cell_1)
                    {
                        temp_cell_2 = element;
                    }
                }
                else
                {
                    if (temp_cell_3[0] == 0)
                    {
                        if ((element != temp_cell_1) && (element != temp_cell_2))
                        {
                            temp_cell_3 = element;
                            break;
                        }
                    }
                    else { break; }
                }
            }
        }
    }
    if ((temp_cell_1[0] * temp_cell_2[0] * temp_cell_3[0]) != 0) { result = true; }
    return(result);
}

void OpenCellFunctionStep1(int flat_ind, std::vector<unsigned int>& identificator, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, double radius_eucledian, double length_x, double length_y, double length_z, int radius, int size_x, int size_y, int size_z, int neighborhood_elements_number)
{
    std::array<int, 3> current_index;
    bool indicator = false;

    current_index[0] = (int)(floor(floor(flat_ind / (size_z - 1)) / (size_y - 1)));
    current_index[1] = (int)((int)(floor(flat_ind / (size_z - 1))) % (size_y - 1));
    current_index[2] = (int)(flat_ind % (size_z - 1));
    indicator = OpenCellCheckEucledian(flat_ind, data, x_block, y_block, z_block, x_centers, y_centers, z_centers, radius_eucledian, length_x, length_y, length_z, radius, size_x, size_y, size_z, neighborhood_elements_number);
    if (indicator == true)
    {
        identificator[flat_ind] = 1;
        data[flat_ind] = -3;
    }
}

void OpenCellFunctionStep2(int flat_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, double radius_eucledian, double length_x, double length_y, double length_z, int radius, int size_x, int size_y, int size_z)
{
    std::array<int, 3> current_index, neighbor_index, change_block_direction = { 0, 0, 0 };
    std::array<bool, 3> change_block = { false, false, false };
    int index = 0, i = 0, j = 0, k = 0;
    current_index[0] = (int)(floor(floor(flat_ind / (size_z - 1)) / (size_y - 1)));
    current_index[1] = (int)((int)(floor(flat_ind / (size_z - 1))) % (size_y - 1));
    current_index[2] = (int)(flat_ind % (size_z - 1));
    if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
    {
//#pragma omp parallel for
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
            j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
            k = (int)(par_ind % (2 * radius + 1)) - radius;
            if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[current_index[0] + i], y_centers[current_index[1] + j], z_centers[current_index[2] + k]) <= radius_eucledian)
            {
                index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
                if ((data[index] == -2) || (data[index] == -1)) { data[index] = -3; }
                
            }
        }
    }
    else
    {
//#pragma omp parallel for
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
            if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
            {
                neighbor_index[0] = abs(size_x - 1 - abs(current_index[0] + i));
                change_block_direction[0] = (int)copysign(1.0, current_index[0] + i);
                change_block[0] = true;
            }
            else
            {
                neighbor_index[0] = current_index[0] + i;
            }
            j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
            if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
            {
                neighbor_index[1] = abs(size_y - 1 - abs(current_index[1] + j));
                change_block_direction[1] = (int)copysign(1.0, current_index[1] + j);
                change_block[1] = true;
            }
            else
            {
                neighbor_index[1] = current_index[1] + j;
            }
            k = (int)(par_ind % (2 * radius + 1)) - radius;
            if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
            {
                neighbor_index[2] = abs(size_z - 1 - abs(current_index[2] + k));
                change_block_direction[2] = (int)copysign(1.0, current_index[2] + k);
                change_block[2] = true;
            }
            else
            {
                neighbor_index[2] = current_index[2] + k;
            }
            if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[neighbor_index[0]] + length_x * change_block_direction[0], y_centers[neighbor_index[1]] + length_y * change_block_direction[1], z_centers[neighbor_index[2]] + length_z * change_block_direction[2]) <= radius_eucledian)
            {
                index = (size_z - 1) * (size_y - 1) * neighbor_index[0] + (size_z - 1) * neighbor_index[1] + neighbor_index[2];
                if (data[index] == -2) { data[index] = -3; }
            }
        }
    }
}

void OpenStructCopmutation(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, double length_x, double length_y, double length_z, int radius, int size_x, int size_y, int size_z, const std::string& neighborhood_check_type)
{
    if (neighborhood_check_type == "by_index")
    {
        #pragma omp parallel for
        for (int flat_ind = 0; flat_ind < data.size(); flat_ind++)
        {
            if (data[flat_ind] == -1)
            {
                data[flat_ind] = -2;
            }
            if (data[flat_ind] == -2)
            {
                OpenCellFunction(flat_ind, data, x_block, y_block, z_block, radius, size_x, size_y, size_z);
            }
        }
    }
    if (neighborhood_check_type == "by_eucledian")
    {
        time_t start, end;
        std::vector<unsigned int> identificator;
        double radius_eucledian = 0.0;
        int radius_new = 0;
        //double radius_eucledian_pre_processing = EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[radius], y_centers[radius], z_centers[radius]);
        //radius_eucledian_pre_processing = (1.0 + 0.01 / radius) * radius_eucledian_pre_processing;
        int neighborhood_elements_number = (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1);
        for (const auto& element : data) { identificator.push_back(0); }
        if (open_cell_structure_computation_type == "per_voxel")
        {
            radius_new = radius + 1;
            radius_eucledian = std::max(EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[radius_new], y_centers[0], z_centers[0]), std::max(EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[0], y_centers[radius_new], z_centers[0]), EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[0], y_centers[0], z_centers[radius_new])));
            radius_eucledian = (1.0 + 0.01 / radius_new) * radius_eucledian;
            time(&start);
            #pragma omp parallel for
            for (int flat_ind = 0; flat_ind < data.size(); flat_ind++)
            {
                if (data[flat_ind] == -1)
                {
                    data[flat_ind] = -2;
                }
                if (data[flat_ind] == -2)
                {
                    OpenCellFunctionStep1(flat_ind, identificator, data, x_block, y_block, z_block, x_centers, y_centers, z_centers, radius_eucledian, length_x, length_y, length_z, radius, size_x, size_y, size_z, neighborhood_elements_number);
                }
            }
            time(&end);
        }
        else
        {
            radius_new = (int)(radius / 2) + 2;
            radius_eucledian = std::max(EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[radius_new], y_centers[0], z_centers[0]), std::max(EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[0], y_centers[radius_new], z_centers[0]), EucledianDistance(x_centers[0], y_centers[0], z_centers[0], x_centers[0], y_centers[0], z_centers[radius_new])));
            radius_eucledian = (1.0 + 0.01 / radius_new) * radius_eucledian;
            time(&start);
            #pragma omp parallel for
            for (int flat_ind = 0; flat_ind < data.size(); flat_ind++)
            {
                if (data[flat_ind] == -1)
                {
                    data[flat_ind] = -2;
                }
                if (data[flat_ind] == -2)
                {
                    OpenCellFunctionStep1(flat_ind, identificator, data, x_block, y_block, z_block, x_centers, y_centers, z_centers, radius_eucledian, length_x, length_y, length_z, radius, size_x, size_y, size_z, neighborhood_elements_number);
                }
            }
            time(&end);
            //PostProcessing::SkeletonWriteToFile(identificator, x_centers, y_centers, z_centers);
            std::cout << "The evaluation time of the 1st stage is " << (double)(difftime(end, start)) << " seconds." << std::endl;
            time(&start);
            for (int flat_ind = 0; flat_ind < data.size(); flat_ind++)
            {
                if (identificator[flat_ind] == 1)
                {
                    OpenCellFunctionStep2(flat_ind, data, x_block, y_block, z_block, x_centers, y_centers, z_centers, radius_eucledian, length_x, length_y, length_z, radius, size_x, size_y, size_z);
                }
            }
            time(&end);
            std::cout << "The evaluation time of the 2nd stage is " << (double)(difftime(end, start)) << " seconds." << std::endl;
        }
    }
}

int CloseStructCheckInternal(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id)
{
    int result = 0;
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    int k = (int)(par_ind % (2 * radius + 1)) - radius;

    int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
    if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
    {
        if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
        {
            result = 1;
        }
    }

    return(result);
}

int CloseStructCheckBoundary(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id)
{
    int result = 0, index_x = 0, index_y = 0, index_z = 0, change_block_x = 0, change_block_y = 0, change_block_z = 0;
    std::array<bool, 3> change_block = { false, false, false };
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
    {
        index_x = abs(size_x - 1 - abs(current_index[0] + i));
        change_block_x = (int)copysign(1.0, current_index[0] + i);
        change_block[0] = true;
    }
    else
    {
        index_x = current_index[0] + i;
    }
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
    {
        index_y = abs(size_y - 1 - abs(current_index[1] + j));
        change_block_y = (int)copysign(1.0, current_index[1] + j);
        change_block[1] = true;
    }
    else
    {
        index_y = current_index[1] + j;
    }
    int k = (int)(par_ind % (2 * radius + 1)) - radius;
    if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
    {
        index_z = abs(size_z - 1 - abs(current_index[2] + k));
        change_block_z = (int)copysign(1.0, current_index[2] + k);
        change_block[2] = true;
    }
    else
    {
        index_z = current_index[2] + k;
    }

    int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
    if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
    {
        if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
        {
            if ((data[index] != current_cell_id) || (x_block[index] != (current_block[0] + change_block_x)) || (y_block[index] != (current_block[1] + change_block_y)) || (z_block[index] != (current_block[2] + change_block_z)))
            {
                result = 1;
            }
        }
    }
    else
    {
        if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
        {
            if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
            {
                result = 1;
            }
        }
    }

    return(result);
}

int CloseStructCheckInternalGeneral(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, long int current_cell_id, double minimum_neighborhood_radius, const std::string& neighborhood_check_type)
{
    int result = 0;
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    int k = (int)(par_ind % (2 * radius + 1)) - radius;

    if (neighborhood_check_type == "by_index")
    {
        int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
        if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
        {
            if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
            {
                result = 1;
            }
        }
    }
    else
    {
        if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[current_index[0] + i], y_centers[current_index[1] + j], z_centers[current_index[2] + k]) <= minimum_neighborhood_radius)
        {
            int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                {
                    result = 1;
                }
            }
        }
    }

    return(result);
}

int CloseStructCheckBoundaryGeneral(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, long int current_cell_id, double minimum_neighborhood_radius, double length_x, double length_y, double length_z, const std::string& neighborhood_check_type)
{
    int result = 0, index_x = 0, index_y = 0, index_z = 0, change_block_x = 0, change_block_y = 0, change_block_z = 0;
    std::array<bool, 3> change_block = { false, false, false };
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
    {
        index_x = abs(size_x - 1 - abs(current_index[0] + i));
        change_block_x = (int)copysign(1.0, current_index[0] + i);
        change_block[0] = true;
    }
    else
    {
        index_x = current_index[0] + i;
    }
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
    {
        index_y = abs(size_y - 1 - abs(current_index[1] + j));
        change_block_y = (int)copysign(1.0, current_index[1] + j);
        change_block[1] = true;
    }
    else
    {
        index_y = current_index[1] + j;
    }
    int k = (int)(par_ind % (2 * radius + 1)) - radius;
    if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
    {
        index_z = abs(size_z - 1 - abs(current_index[2] + k));
        change_block_z = (int)copysign(1.0, current_index[2] + k);
        change_block[2] = true;
    }
    else
    {
        index_z = current_index[2] + k;
    }

    if (neighborhood_check_type == "by_index")
    {
        int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
        if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != (current_block[0] + change_block_x)) || (y_block[index] != (current_block[1] + change_block_y)) || (z_block[index] != (current_block[2] + change_block_z)))
                {
                    result = 1;
                }
            }
        }
        else
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                {
                    result = 1;
                }
            }
        }
    }
    else
    {
        if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[index_x] + length_x * change_block_x, y_centers[index_y] + length_y * change_block_y, z_centers[index_z] + length_z * change_block_z) <= minimum_neighborhood_radius)
        {
            int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
            if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
                {
                    if ((data[index] != current_cell_id) || (x_block[index] != (current_block[0] + change_block_x)) || (y_block[index] != (current_block[1] + change_block_y)) || (z_block[index] != (current_block[2] + change_block_z)))
                    {
                        result = 1;
                    }
                }
            }
            else
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
                {
                    if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                    {
                        result = 1;
                    }
                }
            }
        }
    }
    

    return(result);
}

std::array<int, 5> OpenStructCheckInternalGeneral(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, long int current_cell_id, double minimum_neighborhood_radius, const std::string& neighborhood_check_type)
{
    std::array<int, 5> result = {0, 0, 0, 0, 0};
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    int k = (int)(par_ind % (2 * radius + 1)) - radius;

    if (neighborhood_check_type == "by_index")
    {
        int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
        if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
        {
            if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
            {
                result = {1, data[index], x_block[index], y_block[index], z_block[index]};
            }
        }
    }
    else
    {
        if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[current_index[0] + i], y_centers[current_index[1] + j], z_centers[current_index[2] + k]) <= minimum_neighborhood_radius)
        {
            int index = (size_z - 1) * (size_y - 1) * (current_index[0] + i) + (size_z - 1) * (current_index[1] + j) + (current_index[2] + k);
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                {
                    result = {1, data[index], x_block[index], y_block[index], z_block[index]};
                }
            }
        }
    }
    return(result);
}

std::array<int, 5> OpenStructCheckBoundaryGeneral(int par_ind, std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, long int current_cell_id, double minimum_neighborhood_radius, double length_x, double length_y, double length_z, const std::string& neighborhood_check_type)
{
    std::array<int, 5> result = {0, 0, 0, 0, 0};
    int index_x = 0, index_y = 0, index_z = 0, change_block_x = 0, change_block_y = 0, change_block_z = 0;
    std::array<bool, 3> change_block = { false, false, false };
    int i = (int)(floor(floor(par_ind / (2 * radius + 1)) / (2 * radius + 1))) - radius;
    if ((current_index[0] + i < 0) || (current_index[0] + i > size_x - 2))
    {
        index_x = abs(size_x - 1 - abs(current_index[0] + i));
        change_block_x = (int)copysign(1.0, current_index[0] + i);
        change_block[0] = true;
    }
    else
    {
        index_x = current_index[0] + i;
    }
    int j = (int)((int)(floor(par_ind / (2 * radius + 1))) % (2 * radius + 1)) - radius;
    if ((current_index[1] + j < 0) || (current_index[1] + j > size_y - 2))
    {
        index_y = abs(size_y - 1 - abs(current_index[1] + j));
        change_block_y = (int)copysign(1.0, current_index[1] + j);
        change_block[1] = true;
    }
    else
    {
        index_y = current_index[1] + j;
    }
    int k = (int)(par_ind % (2 * radius + 1)) - radius;
    if ((current_index[2] + k < 0) || (current_index[2] + k > size_z - 2))
    {
        index_z = abs(size_z - 1 - abs(current_index[2] + k));
        change_block_z = (int)copysign(1.0, current_index[2] + k);
        change_block[2] = true;
    }
    else
    {
        index_z = current_index[2] + k;
    }

    if (neighborhood_check_type == "by_index")
    {
        int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
        if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != (current_block[0] + change_block_x)) || (y_block[index] != (current_block[1] + change_block_y)) || (z_block[index] != (current_block[2] + change_block_z)))
                {
                    result = {1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z};
                }
            }
        }
        else
        {
            if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
            {
                if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                {
                    result = {1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z};
                }
            }
        }
    }
    else
    {
        if (EucledianDistance(x_centers[current_index[0]], y_centers[current_index[1]], z_centers[current_index[2]], x_centers[index_x] + length_x * change_block_x, y_centers[index_y] + length_y * change_block_y, z_centers[index_z] + length_z * change_block_z) <= minimum_neighborhood_radius)
        {
            int index = (size_z - 1) * (size_y - 1) * index_x + (size_z - 1) * index_y + index_z;
            if ((change_block[0] == true) || (change_block[1] == true) || (change_block[2] == true))
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
                {
                    if ((data[index] != current_cell_id) || (x_block[index] != (current_block[0] + change_block_x)) || (y_block[index] != (current_block[1] + change_block_y)) || (z_block[index] != (current_block[2] + change_block_z)))
                    {
                        result = {1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z};
                    }
                }
            }
            else
            {
                if ((data[index] != -1) && (data[index] != -2) && (data[index] != -3) && (data[index] != -4) && (data[index] != -5))
                {
                    if ((data[index] != current_cell_id) || (x_block[index] != current_block[0]) || (y_block[index] != current_block[1]) || (z_block[index] != current_block[2]))
                    {
                        result = {1, data[index], x_block[index] - change_block_x, y_block[index] - change_block_y, z_block[index] - change_block_z};
                    }
                }
            }
        }
    }
    return(result);
}

int DetectingNeighbourhood(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, int radius, int size_x, int size_y, int size_z, long int current_cell_id)
{
    int result = 0;
    bool error = false;
    //std::array<int, 728> index, i, j, k;
    int index = 0, i, j, k;
    int flat_current_index = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
    int index_x, index_y, index_z, cur_x_block, cur_y_block, cur_z_block;
    long int temp_cell_id = -1, temp_cell_id_nochange = -1;
    std::array<bool, 3> change_block = { false, false, false };
    int change_block_x = 0, change_block_y = 0, change_block_z = 0, temp_change_block_x = 0, temp_change_block_y = 0, temp_change_block_z = 0, temp_cell_block_x = 0, temp_cell_block_y = 0, temp_cell_block_z = 0;
    int temp_cell_block_nochange_x = 0, temp_cell_block_nochange_y = 0, temp_cell_block_nochange_z = 0;

    cur_x_block = x_block[flat_current_index];
    cur_y_block = y_block[flat_current_index];
    cur_z_block = z_block[flat_current_index];

    std::vector<std::array<int, 4>> surrounding_neighborhood;
    std::array<int, 4> temp;

    if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
    {
#pragma omp parallel for
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            if (result != 1)
            {
                if (CloseStructCheckInternal(par_ind, data, x_block, y_block, z_block, current_block, current_index, radius, size_x, size_y, size_z, current_cell_id) == 1)
                {
                    result = 1;
                }
            }
        }
    }
    else
    {
#pragma omp parallel for
        for (int par_ind = 0; par_ind < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); par_ind++)
        {
            if (result != 1)
            {
                if (CloseStructCheckBoundary(par_ind, data, x_block, y_block, z_block, current_block, current_index, radius, size_x, size_y, size_z, current_cell_id) == 1)
                {
                    result = 1;
                }
            }
        }
    }

    return(result);
}

int DetectingNeighborhoodGeneral(std::vector<long int>& data, std::vector<int>& x_block, std::vector<int>& y_block, std::vector<int>& z_block, std::array<int, 3>& current_block, std::array<int, 3>& current_index, std::vector<double>& x_centers, std::vector<double>& y_centers, std::vector<double>& z_centers, int radius, int size_x, int size_y, int size_z, long int current_cell_id, double minimum_neighborhood_radius, double length_x, double length_y, double length_z, const std::string& structure_type, const std::string& neighborhood_check_type)
{
    int result = 0;
    unsigned int neighborhood_elements_number = (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1);
    bool error = false;
    //std::array<int, 728> index, i, j, k;
    int index = 0, i, j, k, counter = 0;
    int flat_current_index = (size_z - 1) * (size_y - 1) * current_index[0] + (size_z - 1) * current_index[1] + current_index[2];
    int index_x, index_y, index_z, cur_x_block, cur_y_block, cur_z_block;
    long int temp_cell_id = -1, temp_cell_id_nochange = -1;
    std::array<bool, 3> change_block = { false, false, false };
    int change_block_x = 0, change_block_y = 0, change_block_z = 0, temp_change_block_x = 0, temp_change_block_y = 0, temp_change_block_z = 0, temp_cell_block_x = 0, temp_cell_block_y = 0, temp_cell_block_z = 0;
    int temp_cell_block_nochange_x = 0, temp_cell_block_nochange_y = 0, temp_cell_block_nochange_z = 0;

    cur_x_block = x_block[flat_current_index];
    cur_y_block = y_block[flat_current_index];
    cur_z_block = z_block[flat_current_index];

    //for (unsigned int i = 0; i < (2 * radius + 1) * (2 * radius + 1) * (2 * radius + 1); i++) { neighbor_cells_set.push_back({ 0, 0, 0, 0, 0 }); }

    //std::vector<std::array<int, 4>> surrounding_neighborhood;
    std::array<int, 4> temp;
    if (structure_type == "closed")
    {
        if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
        {
#pragma omp parallel for
            for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
            {
                if (result != 1)
                {
                    if (CloseStructCheckInternalGeneral(par_ind, data, x_block, y_block, z_block, current_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, current_cell_id, minimum_neighborhood_radius, neighborhood_check_type) == 1)
                    {
                        result = 1;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        else
        {
#pragma omp parallel for
            for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
            {
                if (result != 1)
                {
                    if (CloseStructCheckBoundaryGeneral(par_ind, data, x_block, y_block, z_block, current_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, current_cell_id, minimum_neighborhood_radius, length_x, length_y, length_z, neighborhood_check_type) == 1)
                    {
                        result = 1;
                    }
                }
                else
                {
                    break;
                }
            }
        }
    }
    else
    {
        std::vector<std::array<int, 5>> neighbor_cells_set;
        std::array<int, 5> neighbor_temp_1, neighbor_temp_2;
        for (unsigned int i = 0; i < neighborhood_elements_number; i++) { neighbor_cells_set.push_back({ 0,0,0,0,0 }); }
        if ((current_index[0] >= radius) && (current_index[0] <= (size_x - radius - 2)) && (current_index[1] >= radius) && (current_index[1] <= (size_y - radius - 2)) && (current_index[2] >= radius) && (current_index[2] <= (size_z - radius - 2)))
        {
#pragma omp parallel for
            for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
            {
                neighbor_cells_set[par_ind] = OpenStructCheckInternalGeneral(par_ind, data, x_block, y_block, z_block, current_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, current_cell_id, minimum_neighborhood_radius, neighborhood_check_type);
            }
        }
        else
        {
#pragma omp parallel for
            for (int par_ind = 0; par_ind < neighborhood_elements_number; par_ind++)
            {
                neighbor_cells_set[par_ind] = OpenStructCheckBoundaryGeneral(par_ind, data, x_block, y_block, z_block, current_block, current_index, x_centers, y_centers, z_centers, radius, size_x, size_y, size_z, current_cell_id, minimum_neighborhood_radius, length_x, length_y, length_z, neighborhood_check_type);
            }
        }
        for (unsigned int i = 0; i < neighborhood_elements_number; i++)
        { 
            neighbor_temp_1 = neighbor_cells_set[i];
            if (neighbor_temp_1[0] == 1)
            {
                for (unsigned int i = 0; i < neighborhood_elements_number; i++)
                {
                    neighbor_temp_2 = neighbor_cells_set[i];
                    if (neighbor_temp_2[0] != 0)
                    {
                        if ((neighbor_temp_1[1] != neighbor_temp_2[1]) || (neighbor_temp_1[2] != neighbor_temp_2[2]) || (neighbor_temp_1[3] != neighbor_temp_2[3]) || (neighbor_temp_1[4] != neighbor_temp_2[4]))
                        {
                            result = 1;
                        }
                    }
                    if (result == 1) { break; }
                }
                break;
            }
        }
    }

    return(result);
}

std::priority_queue<std::array<double, 9>> ListInitialization(std::vector<std::array<double, 3>>& seeds_centers, std::vector<double>& centers_x, std::vector<double>& centers_y, std::vector<double>& centers_z, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double phi, double theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    int seeds_number = seeds_centers.size();
    int x_number = x.size();
    int y_number = y.size();
    int z_number = z.size();
    std::priority_queue<std::array<double, 9>> Queue;
    std::array<double, 9> temp_for_push = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::array<int, 3> corresponding_center = { 0, 0, 0 };
    int seed_id = -1;
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
            dist = Distance(seed[0], seed[1], seed[2], centers_x[corresponding_center[0]], centers_y[corresponding_center[1]], centers_z[corresponding_center[2]], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
            temp_for_push = { -dist, (double)voxel_counter, (double)corresponding_center[0], (double)corresponding_center[1], (double)corresponding_center[2], (double)seed_id, 0.0, 0.0, 0.0 };
            Queue.push(temp_for_push);
            //std::cout << Queue.size() << std::endl;
            voxel_counter = voxel_counter + 1;
        }
    }

    return Queue;
}

long int GetCellIdentifier(std::array<int, 3>& block_identifier, int seed_identifier_1, int seeds_number)
{
    long int cell_identifier = 0;
    int max_block = 15;

    if ((block_identifier[0] == 0) && (block_identifier[1] == 0) && (block_identifier[2] == 0))
    {
        cell_identifier = seed_identifier_1;
    }
    else
    {
        double seed_identifier = seed_identifier_1 + 5;
        if ((block_identifier[0] >= 0) && (block_identifier[1] >= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = (max_block * (max_block - 1) * block_identifier[0] + (max_block - 1) * block_identifier[1] + block_identifier[2] + seeds_number) * seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] <= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block * (max_block - 1) * block_identifier[0] + (max_block - 1) * block_identifier[1] + block_identifier[2] - seeds_number) * seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] >= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] <= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = -(max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seed_identifier;
        }

        if ((block_identifier[0] >= 0) && (block_identifier[1] <= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = (max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seeds_number * seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] >= 0) && (block_identifier[2] <= 0))
        {
            cell_identifier = -(max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seeds_number * seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] >= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = (max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seeds_number * seeds_number * seed_identifier;
        }

        if ((block_identifier[0] <= 0) && (block_identifier[1] <= 0) && (block_identifier[2] >= 0))
        {
            cell_identifier = -(max_block * (max_block - 1) * abs(block_identifier[0]) + (max_block - 1) * abs(block_identifier[1]) + abs(block_identifier[2]) + seeds_number) * seeds_number * seeds_number * seeds_number * seed_identifier;
        }
    }

    return cell_identifier;
}

std::list<std::vector<long int>> MethodGrowth(std::vector<std::array<double, 3>>& seeds, double phi, double theta, double length_x, double length_y, double length_z, int size_x, int size_y, int size_z, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    std::vector<double> grid_x = GridComputation(length_x, size_x);
    std::vector<double> grid_y = GridComputation(length_y, size_y);
    std::vector<double> grid_z = GridComputation(length_z, size_z);
    std::vector<double> grid_centers_x = GridCentersComputation(grid_x);
    std::vector<double> grid_centers_y = GridCentersComputation(grid_y);
    std::vector<double> grid_centers_z = GridCentersComputation(grid_z);
    std::vector<long int> result_full;
    std::vector<int> x_block, y_block, z_block;
    std::list<std::vector<long int>> final, temporary;

    for (int i = 0; i < grid_centers_x.size(); i++)
    {
        for (int j = 0; j < grid_centers_y.size(); j++)
        {
            for (int k = 0; k < grid_centers_z.size(); k++)
            {
                result_full.push_back(-1);
                x_block.push_back(0);
                y_block.push_back(0);
                z_block.push_back(0);
            }
        }
    }

    std::priority_queue<std::array<double, 9>> Queue = ListInitialization(seeds, grid_centers_x, grid_centers_y, grid_centers_z, grid_x, grid_y, grid_z, phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);

    std::array<int, 3> left_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> right_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> front_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> back_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> up_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> down_neighbour_block = { 0, 0, 0 };

    std::array<int, 3> left_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> right_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> front_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> back_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> up_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> down_neighbour_center = { 0, 0, 0 };

    std::array<double, 9> current_info, temp_to_push;
    std::array<double, 3> temp_point, seed = { 0.0, 0.0, 0.0 };
    std::array<int, 3> current_cell;
    std::array<int, 3> current_block;
    double dist = 0.0;
    int current_seed_id, current_index, index, neighbourhood, current_neighbourhood = 0, iterator = 0;
    long int voxel_id = Queue.size();
    int initial_id = voxel_id;

    double x_last = grid_x.back();
    double y_last = grid_y.back();
    double z_last = grid_z.back();

    while (Queue.size() != 0)
    {
        current_info = Queue.top();
        Queue.pop();
        current_cell[0] = (int)current_info[2];
        current_cell[1] = (int)current_info[3];
        current_cell[2] = (int)current_info[4];
        current_index = (size_z - 1) * (size_y - 1) * current_cell[0] + (size_z - 1) * current_cell[1] + current_cell[2]; //look here again
        current_seed_id = (int)current_info[5];
        current_block[0] = (int)current_info[6];
        current_block[1] = (int)current_info[7];
        current_block[2] = (int)current_info[8];
        seed = seeds[current_seed_id];
        //std::cout << "result = " << current_index << std::endl;

        //current_neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, current_cell, edge_radius, size_x, size_y, size_z, result_full[current_index]);
        if (result_full[current_index] == -1)
        {
            current_neighbourhood = DetectingNeighbourhood(result_full, x_block, y_block, z_block, current_block, current_cell, edge_radius, size_x, size_y, size_z, current_seed_id);

            if (current_neighbourhood == 0)
            {
                result_full[current_index] = current_seed_id;
                x_block[current_index] = current_block[0];
                y_block[current_index] = current_block[1];
                z_block[current_index] = current_block[2];

                //CheckContiniouty(result_full, x_block, y_block, z_block, current_index, size_x, size_y, size_z, voxel_id);

                //index = (size_z - 1)*(size_y - 1)*0 + (size_z - 1)*74 + 75;
                //std::cout << "well  (" << result_full[index] << ", " << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << std::endl;

                //index = (size_z - 1)*(size_y - 1)*1 + (size_z - 1)*74 + 75;
                //std::cout << "well  (" << result_full[index] << ", " << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << std::endl;

                left_neighbour_block[0] = current_block[0];
                left_neighbour_block[1] = current_block[1];
                left_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * left_neighbour_center[0] + (size_z - 1) * left_neighbour_center[1] + left_neighbour_center[2];
                //x_block[index] = left_neighbour_block[0];
                //y_block[index] = left_neighbour_block[1];
                //z_block[index] = left_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[left_neighbour_center[0]] + left_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[left_neighbour_center[1]] + left_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[left_neighbour_center[2]] + left_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)left_neighbour_center[0], (double)left_neighbour_center[1], (double)left_neighbour_center[2], (double)current_seed_id, (double)left_neighbour_block[0], (double)left_neighbour_block[1], (double)left_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "left yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                right_neighbour_block[0] = current_block[0];
                right_neighbour_block[1] = current_block[1];
                right_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * right_neighbour_center[0] + (size_z - 1) * right_neighbour_center[1] + right_neighbour_center[2];
                //x_block[index] = right_neighbour_block[0];
                //y_block[index] = right_neighbour_block[1];
                //z_block[index] = right_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[right_neighbour_center[0]] + right_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[right_neighbour_center[1]] + right_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[right_neighbour_center[2]] + right_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)right_neighbour_center[0], (double)right_neighbour_center[1], (double)right_neighbour_center[2], (double)current_seed_id, (double)right_neighbour_block[0], (double)right_neighbour_block[1], (double)right_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "right yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                front_neighbour_block[0] = current_block[0];
                front_neighbour_block[1] = current_block[1];
                front_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * front_neighbour_center[0] + (size_z - 1) * front_neighbour_center[1] + front_neighbour_center[2];
                //x_block[index] = front_neighbour_block[0];
                //y_block[index] = front_neighbour_block[1];
                //z_block[index] = front_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[front_neighbour_center[0]] + front_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[front_neighbour_center[1]] + front_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[front_neighbour_center[2]] + front_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)front_neighbour_center[0], (double)front_neighbour_center[1], (double)front_neighbour_center[2], (double)current_seed_id, (double)front_neighbour_block[0], (double)front_neighbour_block[1], (double)front_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "front yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                back_neighbour_block[0] = current_block[0];
                back_neighbour_block[1] = current_block[1];
                back_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * back_neighbour_center[0] + (size_z - 1) * back_neighbour_center[1] + back_neighbour_center[2];
                //x_block[index] = back_neighbour_block[0];
                //y_block[index] = back_neighbour_block[1];
                //z_block[index] = back_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[back_neighbour_center[0]] + back_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[back_neighbour_center[1]] + back_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[back_neighbour_center[2]] + back_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)back_neighbour_center[0], (double)back_neighbour_center[1], (double)back_neighbour_center[2], (double)current_seed_id, (double)back_neighbour_block[0], (double)back_neighbour_block[1], (double)back_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "back yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                up_neighbour_block[0] = current_block[0];
                up_neighbour_block[1] = current_block[1];
                up_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * up_neighbour_center[0] + (size_z - 1) * up_neighbour_center[1] + up_neighbour_center[2];
                //x_block[index] = up_neighbour_block[0];
                //y_block[index] = up_neighbour_block[1];
                //z_block[index] = up_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[up_neighbour_center[0]] + up_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[up_neighbour_center[1]] + up_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[up_neighbour_center[2]] + up_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)up_neighbour_center[0], (double)up_neighbour_center[1], (double)up_neighbour_center[2], (double)current_seed_id, (double)up_neighbour_block[0], (double)up_neighbour_block[1], (double)up_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "up yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")"<< index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                down_neighbour_block[0] = current_block[0];
                down_neighbour_block[1] = current_block[1];
                down_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * down_neighbour_center[0] + (size_z - 1) * down_neighbour_center[1] + down_neighbour_center[2];
                //x_block[index] = down_neighbour_block[0];
                //y_block[index] = down_neighbour_block[1];
                //z_block[index] = down_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[down_neighbour_center[0]] + down_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[down_neighbour_center[1]] + down_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[down_neighbour_center[2]] + down_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)down_neighbour_center[0], (double)down_neighbour_center[1], (double)down_neighbour_center[2], (double)current_seed_id, (double)down_neighbour_block[0], (double)down_neighbour_block[1], (double)down_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }

            }
            else
            {
                if (current_neighbourhood == 1)
                {
                    result_full[current_index] = -2;
                }
                else
                {
                    result_full[current_index] = -(current_neighbourhood + 1);
                }
            }

            //std::cout << "down yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;

            //std::cout << "current block number yo  (" << x_block[current_index] << ", " << y_block[current_index] << ", " << z_block[current_index] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << left_neighbour_block[0] << ", " << left_neighbour_block[1] << ", " << left_neighbour_block[2] << ")   " << "(" << right_neighbour_block[0] << ", " << right_neighbour_block[1] << ", " << right_neighbour_block[2] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << front_neighbour_block[0] << ", " << front_neighbour_block[1] << ", " << front_neighbour_block[2] << ")   " << "(" << back_neighbour_block[0] << ", " << back_neighbour_block[1] << ", " << back_neighbour_block[2] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << down_neighbour_block[0] << ", " << down_neighbour_block[1] << ", " << down_neighbour_block[2] << ")   " << "(" << up_neighbour_block[0] << ", " << up_neighbour_block[1] << ", " << up_neighbour_block[2] << ")" << std::endl;

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

    OpenStructCopmutation(result_full, x_block, y_block, z_block, grid_centers_x, grid_centers_y, grid_centers_z, length_x, length_y, length_z, (int)(edge_radius / 2) + 1, size_x, size_y, size_z, neighborhood_check_type_const);

    //TestIndexation(result_full, x_block, y_block, z_block);

    final.push_back(result_full);

    return final;
}

std::list<std::vector<long int>> NewMethodGrowth(std::vector<std::array<double, 3>>& seeds, double phi, double theta, double length_x, double length_y, double length_z, int size_x, int size_y, int size_z, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry, const std::string& structure_type, const std::string& neighborhood_check_type)
{
    std::vector<double> grid_x = GridComputation(length_x, size_x);
    std::vector<double> grid_y = GridComputation(length_y, size_y);
    std::vector<double> grid_z = GridComputation(length_z, size_z);
    std::vector<double> grid_centers_x = GridCentersComputation(grid_x);
    std::vector<double> grid_centers_y = GridCentersComputation(grid_y);
    std::vector<double> grid_centers_z = GridCentersComputation(grid_z);
    std::vector<long int> result_full;
    std::vector<int> x_block, y_block, z_block;
    std::list<std::vector<long int>> final, temporary;
    double minimum_neighborhood_radius = 0;
    int counter = 0;

    for (int i = 0; i < grid_centers_x.size(); i++)
    {
        for (int j = 0; j < grid_centers_y.size(); j++)
        {
            for (int k = 0; k < grid_centers_z.size(); k++)
            {
                result_full.push_back(-1);
                x_block.push_back(0);
                y_block.push_back(0);
                z_block.push_back(0);
            }
        }
    }

    minimum_neighborhood_radius = std::max(EucledianDistance(grid_centers_x[0], grid_centers_y[0], grid_centers_z[0], grid_centers_x[edge_radius], grid_centers_y[0], grid_centers_z[0]), std::max(EucledianDistance(grid_centers_x[0], grid_centers_y[0], grid_centers_z[0], grid_centers_x[0], grid_centers_y[edge_radius], grid_centers_z[0]), EucledianDistance(grid_centers_x[0], grid_centers_y[0], grid_centers_z[0], grid_centers_x[0], grid_centers_y[0], grid_centers_z[edge_radius])));
    minimum_neighborhood_radius = (1.0 + 0.1 / edge_radius) * minimum_neighborhood_radius;

    std::priority_queue<std::array<double, 9>> Queue = ListInitialization(seeds, grid_centers_x, grid_centers_y, grid_centers_z, grid_x, grid_y, grid_z, phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);

    std::array<int, 3> left_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> right_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> front_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> back_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> up_neighbour_block = { 0, 0, 0 };
    std::array<int, 3> down_neighbour_block = { 0, 0, 0 };

    std::array<int, 3> left_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> right_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> front_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> back_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> up_neighbour_center = { 0, 0, 0 };
    std::array<int, 3> down_neighbour_center = { 0, 0, 0 };

    std::array<double, 9> current_info, temp_to_push;
    std::array<double, 3> temp_point, seed = { 0.0, 0.0, 0.0 };
    std::array<int, 3> current_cell;
    std::array<int, 3> current_block;
    double dist = 0.0;
    int current_seed_id, current_index, index, neighbourhood, current_neighbourhood = 0, iterator = 0;
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
        current_index = (size_z - 1) * (size_y - 1) * current_cell[0] + (size_z - 1) * current_cell[1] + current_cell[2]; //look here again
        current_seed_id = (int)current_info[5];
        current_block[0] = (int)current_info[6];
        current_block[1] = (int)current_info[7];
        current_block[2] = (int)current_info[8];
        seed = seeds[current_seed_id];
        //std::cout << "result = " << current_index << std::endl;

        //current_neighbourhood = DetectingNeighbourhood(result_full, grid_x, grid_y, grid_z, current_cell, edge_radius, size_x, size_y, size_z, result_full[current_index]);
        if (result_full[current_index] == -1)
        {
            current_neighbourhood = DetectingNeighborhoodGeneral(result_full, x_block, y_block, z_block, current_block, current_cell, grid_centers_x, grid_centers_y, grid_centers_z, edge_radius, size_x, size_y, size_z, current_seed_id, minimum_neighborhood_radius, length_x, length_y, length_z, structure_type, neighborhood_check_type);

            if (current_neighbourhood == 0)
            {
                result_full[current_index] = current_seed_id;
                x_block[current_index] = current_block[0];
                y_block[current_index] = current_block[1];
                z_block[current_index] = current_block[2];

                left_neighbour_block[0] = current_block[0];
                left_neighbour_block[1] = current_block[1];
                left_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * left_neighbour_center[0] + (size_z - 1) * left_neighbour_center[1] + left_neighbour_center[2];
                //x_block[index] = left_neighbour_block[0];
                //y_block[index] = left_neighbour_block[1];
                //z_block[index] = left_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[left_neighbour_center[0]] + left_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[left_neighbour_center[1]] + left_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[left_neighbour_center[2]] + left_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)left_neighbour_center[0], (double)left_neighbour_center[1], (double)left_neighbour_center[2], (double)current_seed_id, (double)left_neighbour_block[0], (double)left_neighbour_block[1], (double)left_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "left yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                right_neighbour_block[0] = current_block[0];
                right_neighbour_block[1] = current_block[1];
                right_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * right_neighbour_center[0] + (size_z - 1) * right_neighbour_center[1] + right_neighbour_center[2];
                //x_block[index] = right_neighbour_block[0];
                //y_block[index] = right_neighbour_block[1];
                //z_block[index] = right_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[right_neighbour_center[0]] + right_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[right_neighbour_center[1]] + right_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[right_neighbour_center[2]] + right_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)right_neighbour_center[0], (double)right_neighbour_center[1], (double)right_neighbour_center[2], (double)current_seed_id, (double)right_neighbour_block[0], (double)right_neighbour_block[1], (double)right_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "right yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                front_neighbour_block[0] = current_block[0];
                front_neighbour_block[1] = current_block[1];
                front_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * front_neighbour_center[0] + (size_z - 1) * front_neighbour_center[1] + front_neighbour_center[2];
                //x_block[index] = front_neighbour_block[0];
                //y_block[index] = front_neighbour_block[1];
                //z_block[index] = front_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[front_neighbour_center[0]] + front_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[front_neighbour_center[1]] + front_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[front_neighbour_center[2]] + front_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)front_neighbour_center[0], (double)front_neighbour_center[1], (double)front_neighbour_center[2], (double)current_seed_id, (double)front_neighbour_block[0], (double)front_neighbour_block[1], (double)front_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "front yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                back_neighbour_block[0] = current_block[0];
                back_neighbour_block[1] = current_block[1];
                back_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * back_neighbour_center[0] + (size_z - 1) * back_neighbour_center[1] + back_neighbour_center[2];
                //x_block[index] = back_neighbour_block[0];
                //y_block[index] = back_neighbour_block[1];
                //z_block[index] = back_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[back_neighbour_center[0]] + back_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[back_neighbour_center[1]] + back_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[back_neighbour_center[2]] + back_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)back_neighbour_center[0], (double)back_neighbour_center[1], (double)back_neighbour_center[2], (double)current_seed_id, (double)back_neighbour_block[0], (double)back_neighbour_block[1], (double)back_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "back yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                up_neighbour_block[0] = current_block[0];
                up_neighbour_block[1] = current_block[1];
                up_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * up_neighbour_center[0] + (size_z - 1) * up_neighbour_center[1] + up_neighbour_center[2];
                //x_block[index] = up_neighbour_block[0];
                //y_block[index] = up_neighbour_block[1];
                //z_block[index] = up_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[up_neighbour_center[0]] + up_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[up_neighbour_center[1]] + up_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[up_neighbour_center[2]] + up_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)up_neighbour_center[0], (double)up_neighbour_center[1], (double)up_neighbour_center[2], (double)current_seed_id, (double)up_neighbour_block[0], (double)up_neighbour_block[1], (double)up_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }
                //std::cout << "up yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")"<< index << std::endl;
            //-------------------------------------------------------------------------------------------------------------------------------------
                down_neighbour_block[0] = current_block[0];
                down_neighbour_block[1] = current_block[1];
                down_neighbour_block[2] = current_block[2];
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
                index = (size_z - 1) * (size_y - 1) * down_neighbour_center[0] + (size_z - 1) * down_neighbour_center[1] + down_neighbour_center[2];
                //x_block[index] = down_neighbour_block[0];
                //y_block[index] = down_neighbour_block[1];
                //z_block[index] = down_neighbour_block[2];
                if (result_full[index] == -1)
                {
                    temp_point[0] = grid_centers_x[down_neighbour_center[0]] + down_neighbour_block[0] * x_last;
                    temp_point[1] = grid_centers_y[down_neighbour_center[1]] + down_neighbour_block[1] * y_last;
                    temp_point[2] = grid_centers_z[down_neighbour_center[2]] + down_neighbour_block[2] * z_last;
                    dist = Distance(seed[0], seed[1], seed[2], temp_point[0], temp_point[1], temp_point[2], phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
                    temp_to_push = { -dist, (double)voxel_id, (double)down_neighbour_center[0], (double)down_neighbour_center[1], (double)down_neighbour_center[2], (double)current_seed_id, (double)down_neighbour_block[0], (double)down_neighbour_block[1], (double)down_neighbour_block[2] };
                    Queue.push(temp_to_push);
                    voxel_id = voxel_id + 1;
                }

            }
            else
            {
                if (current_neighbourhood == 1)
                {
                    result_full[current_index] = -2;
                }
                else
                {
                    result_full[current_index] = -(current_neighbourhood + 1);
                }
            }

            //std::cout << "down yo  (" << x_block[index] << ", " << y_block[index] << ", " << z_block[index] << ")" << index << std::endl;

            //std::cout << "current block number yo  (" << x_block[current_index] << ", " << y_block[current_index] << ", " << z_block[current_index] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << left_neighbour_block[0] << ", " << left_neighbour_block[1] << ", " << left_neighbour_block[2] << ")   " << "(" << right_neighbour_block[0] << ", " << right_neighbour_block[1] << ", " << right_neighbour_block[2] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << front_neighbour_block[0] << ", " << front_neighbour_block[1] << ", " << front_neighbour_block[2] << ")   " << "(" << back_neighbour_block[0] << ", " << back_neighbour_block[1] << ", " << back_neighbour_block[2] << ")" << std::endl;
            //std::cout << "current neighborhood  (" << down_neighbour_block[0] << ", " << down_neighbour_block[1] << ", " << down_neighbour_block[2] << ")   " << "(" << up_neighbour_block[0] << ", " << up_neighbour_block[1] << ", " << up_neighbour_block[2] << ")" << std::endl;

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
        counter = counter + 1;
    }

    std::cout << "open-cell structure's thickness radius = " << edge_radius + 1 << std::endl << std::endl;

    //std::vector<long int> result_temp = result_full;
    //TestIndexation(result_temp, x_block, y_block, z_block);

    //temporary.push_back(result_temp);

    OpenStructCopmutation(result_full, x_block, y_block, z_block, grid_centers_x, grid_centers_y, grid_centers_z, length_x, length_y, length_z, edge_radius, size_x, size_y, size_z, neighborhood_check_type);
    
    final.push_back(result_full);

    //temporary.push_back(result_full);
    //PrintToFileLongInt(temporary, std::string("lolkekcheburek"));

    return final;
}

void InitializeGrowth(double phi, double theta, double max_length, int max_dimension, int edge_radius, double regular_distance_factor, double semiregular_distance_factor, const std::string& seed_distribution, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    std::string name_structure_sample = NameInitialisation(phi, theta, edge_radius, regular_distance_factor, semiregular_distance_factor, seed_distribution, distance_type, periodic_distance_type, periodic_distance_symmetry);

    DistancePrintToFile(phi, theta, max_length, max_dimension, 400, 200, regular_distance_factor, semiregular_distance_factor, seed_distribution, distance_type, periodic_distance_type, periodic_distance_symmetry, folder, name_structure_sample);

    std::vector<std::array<double, 3>> seeds;
    std::array<double, 3> temp;
    std::array<int, 3> sample_size = { 0, 0, 0 };
    std::array<double, 3> sample_length = { 0.0, 0.0, 0.0 };

    std::cout << "          ________________________________________________________________________________" << std::endl << std::endl;

    std::cout << "the distance position in the space: (rotational angle around z-axis = " << phi << " rad, rotational angle around y-axis = " << theta << " rad)," << std::endl
       << std::endl <<"max sample length = " << max_length << ", max sample dimension = " << max_dimension << ", edge radius (structure thickness) = " << edge_radius << ", seed distribution type: " << seed_distribution << std::endl << std::endl;

    std::cout << "regular distance sample multiplication factor = " << regular_distance_factor << ", semiregular distance multiplication factor = " << semiregular_distance_factor << ", distance type: " << distance_type << ", periodic_distance_type: " << periodic_distance_type << ", periodic_distance_symmetry: " << periodic_distance_symmetry << std::endl << std::endl;

    CellGridInitialization(seeds, sample_length, sample_size, max_length, max_dimension, seed_distribution, distance_type);


    std::cout << "the sample size is (" << sample_size[0] << " x " << sample_size[1] << " x " << sample_size[2] << ") elements" << std::endl << std::endl;
    std::cout << "the sample length is (" << sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << ") units" << std::endl << std::endl;
    std::cout << "the seeds coordinates:" << std::endl;
    for (auto element : seeds)
    {
        temp = element;
        std::cout << "{ " << temp[0] << ", " << temp[1] << ", " << temp[2] << " }" << std::endl;
    }

    //CellGridExpansion(seeds, sample_length);

    std::cout << std::endl << "Computing the structure..." << std::endl << std::endl;

    //std::list<std::vector<long int>> solution = MethodGrowth(seeds, phi, theta, sample_length[0], sample_length[1], sample_length[2], sample_size[0], sample_size[1], sample_size[2], edge_radius, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
    std::list<std::vector<long int>> solution = NewMethodGrowth(seeds, phi, theta, sample_length[0], sample_length[1], sample_length[2], sample_size[0], sample_size[1], sample_size[2], edge_radius, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry, structure_type_const, neighborhood_check_type_const);
    std::list<std::vector<unsigned int>> processed_solution = SolutionPostProcessing(solution);

    
    PrintToFileUnsignedInt(processed_solution, folder, name_structure_sample + std::string("_hom"));
    std::cout << std::endl << "Writing original data into the file..." << std::endl << std::endl;

    //PrintToFileLongInt(solution, std::string("test"));
    //PrintToFileUnsignedInt(processed_solution, file_id);
    /*std::list<std::vector<unsigned int>> solution_list_test;
    std::vector<unsigned int> solution_test;
    std::vector<long int> solution_for_test = solution.front();
    for (const auto& element : solution_for_test)
    {
        solution_test.push_back((unsigned int)element);
    }
    solution_list_test.push_back(solution_test);
    solution_list_test.push_back(solution_test);
    
    WriteMetaData(solution_list_test, folder, std::string("test"), sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);*/
    WriteMetaData(processed_solution, folder, name_structure_sample, sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);

    if (distance_type == std::string("multiple"))
    {
        PrintMultiDistInterpolationToFile(400, 200, 25, folder, name_structure_sample);
        std::cout << std::endl << "Writing extended data into the file..." << std::endl << std::endl;
        std::vector<unsigned int> closed_structure = processed_solution.front(), open_structure = processed_solution.back();
        processed_solution.pop_back();
        processed_solution.pop_front();
        closed_structure = StructureExpansion(closed_structure, structure_expansion_rate_x, structure_expansion_rate_y, structure_expansion_rate_z, sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);
        open_structure = StructureExpansion(open_structure, structure_expansion_rate_x, structure_expansion_rate_y, structure_expansion_rate_z, sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);
        //closed_structure = StructureEdging(closed_structure, 3 * (sample_size[0] - 1), 3 * (sample_size[1] - 1), 3 * (sample_size[2] - 1));
        //open_structure = StructureEdging(open_structure, 3 * (sample_size[0] - 1), 3 * (sample_size[1] - 1), 3 * (sample_size[2] - 1));
        //std::cout << std::endl << open_structure.size() << std::endl;
        closed_structure = StructureEdgingManual(closed_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, structure_expansion_rate_x * (sample_size[0] - 1), structure_expansion_rate_y * (sample_size[1] - 1), structure_expansion_rate_z * (sample_size[2] - 1));
        open_structure = StructureEdgingManual(open_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, structure_expansion_rate_x * (sample_size[0] - 1), structure_expansion_rate_y * (sample_size[1] - 1), structure_expansion_rate_z * (sample_size[2] - 1));
        closed_structure = StructureEdging(closed_structure, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        open_structure = StructureEdging(open_structure, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        //closed_structure = StructureEdgingManual(closed_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, (sample_size[0] - 1), (sample_size[1] - 1), (sample_size[2] - 1));
        //open_structure = StructureEdgingManual(open_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, (sample_size[0] - 1), (sample_size[1] - 1), (sample_size[2] - 1));
        //closed_structure = StructureEdging(closed_structure, (sample_size[0] - 1) + 2 * structure_edge_expansion_x, (sample_size[1] - 1) + 2 * structure_edge_expansion_y, (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        //open_structure = StructureEdging(open_structure, (sample_size[0] - 1) + 2 * structure_edge_expansion_x, (sample_size[1] - 1) + 2 * structure_edge_expansion_y, (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        processed_solution.push_front(closed_structure);
        processed_solution.push_back(open_structure);

        name_structure_sample = "expanded_" + name_structure_sample;
        WriteMetaData(processed_solution, folder + "/multi_distance", name_structure_sample, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x + 2, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y + 2, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z + 2);
        //WriteMetaData(processed_solution, folder + "/multi_distance", name_structure_sample, (sample_size[0] - 1) + 2 * structure_edge_expansion_x + 2, (sample_size[1] - 1) + 2 * structure_edge_expansion_y + 2, (sample_size[2] - 1) + 2 * structure_edge_expansion_z + 2);
    }
    else
    {
        std::cout << std::endl << "Writing extended data into the file..." << std::endl << std::endl;

        std::vector<unsigned int> closed_structure = processed_solution.front(), open_structure = processed_solution.back();
        processed_solution.pop_back();
        processed_solution.pop_front();
        closed_structure = StructureExpansion(closed_structure, structure_expansion_rate_x, structure_expansion_rate_y, structure_expansion_rate_z, sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);
        open_structure = StructureExpansion(open_structure, structure_expansion_rate_x, structure_expansion_rate_y, structure_expansion_rate_z, sample_size[0] - 1, sample_size[1] - 1, sample_size[2] - 1);
        //closed_structure = StructureEdging(closed_structure, 3 * (sample_size[0] - 1), 3 * (sample_size[1] - 1), 3 * (sample_size[2] - 1));
        //open_structure = StructureEdging(open_structure, 3 * (sample_size[0] - 1), 3 * (sample_size[1] - 1), 3 * (sample_size[2] - 1));
        //std::cout << std::endl << open_structure.size() << std::endl;
        closed_structure = StructureEdgingManual(closed_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, structure_expansion_rate_x * (sample_size[0] - 1), structure_expansion_rate_y * (sample_size[1] - 1), structure_expansion_rate_z * (sample_size[2] - 1));
        open_structure = StructureEdgingManual(open_structure, structure_edge_expansion_x, structure_edge_expansion_y, structure_edge_expansion_z, structure_expansion_rate_x * (sample_size[0] - 1), structure_expansion_rate_y * (sample_size[1] - 1), structure_expansion_rate_z * (sample_size[2] - 1));
        closed_structure = StructureEdging(closed_structure, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        open_structure = StructureEdging(open_structure, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z);
        //std::cout << std::endl << open_structure.size() << std::endl;
        processed_solution.push_front(closed_structure);
        processed_solution.push_back(open_structure);

        name_structure_sample = "expanded_" + name_structure_sample;
        WriteMetaData(processed_solution, folder + "/expanded_structure", name_structure_sample, structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x + 2, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y + 2, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z + 2);

        std::cout << std::endl << "          ________________________________________________________________________________" << std::endl << std::endl;

        /*std::array<int, 3> resulting_structure_size = { structure_expansion_rate_x * (sample_size[0] - 1) + 2 * structure_edge_expansion_x + 2, structure_expansion_rate_y * (sample_size[1] - 1) + 2 * structure_edge_expansion_y + 2, structure_expansion_rate_z * (sample_size[2] - 1) + 2 * structure_edge_expansion_z + 2 };
        std::vector<unsigned int> shaped_structure = PostProcessing::StructureIntoVolume(processed_solution, sample_length, resulting_structure_size, edge_radius, std::string("sphere"));
        std::list<std::vector<unsigned int>> list_shaped_structure;
        list_shaped_structure.push_back(shaped_structure);
        std::cout << processed_solution.front().size() << " " << list_shaped_structure.front().size() << std::endl;
        WriteMetaData(list_shaped_structure, folder + "/expanded_structure", "shaped" + name_structure_sample, resulting_structure_size[0], resulting_structure_size[1], resulting_structure_size[2]);*/
    }
}

void InitializeRandomExploration(double max_length, int max_dimension, int edge_radius, unsigned int ensemble_size, const std::string& seed_distribution, const std::string& distance_type)
{
    std::vector<double> random_phi, random_theta, random_distance_factor;
    std::vector<unsigned int> random_lattice_type, random_periodic_distance_type, random_periodic_distance_symmetry;
    double min_factor_value = 0.5;

    srand(time(NULL));

    for (unsigned int index = 0; index < ensemble_size; index++)
    {
        random_phi.push_back(2.0 * pi * ((double)rand() / (double)RAND_MAX));
        random_theta.push_back(pi * ((double)rand() / (double)RAND_MAX));
        random_distance_factor.push_back(std::max(min_factor_value, (double)(rand() % 3)));
        random_periodic_distance_symmetry.push_back(rand() % 8);
        random_lattice_type.push_back((int)(rand() % 4));
        if (random_periodic_distance_symmetry[index] <= 4) { random_periodic_distance_type.push_back(0); }
        else { random_periodic_distance_type.push_back(1); }
    }

    for (unsigned int index = 0; index < ensemble_size; index++)
    {
        std::cout << std::endl << std::endl << "*****   Structure number " << (index + 1) << " out of " << ensemble_size << "   *****" << std::endl << std::endl;
        InitializeGrowth(random_phi[index], random_theta[index], max_length, max_dimension, edge_radius, random_distance_factor[index], random_distance_factor[index], seed_distribution_type_list[(int)(3 + random_lattice_type[index])], distance_type, periodic_distance_type_list[random_periodic_distance_type[index]], periodic_distance_symmetry_list[random_periodic_distance_symmetry[index]]);
    }
}