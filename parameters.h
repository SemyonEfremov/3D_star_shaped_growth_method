#pragma once
#include "Classes.h"
#include <array>


const std::string folder = "metadata/SymmetryTests"; // + additional internal folders: .../expanded_structure, .../distance, .../multi_distance and .../homogenization

const float pi = 3.141592653589793238462643;
const double golden_ratio = (1.0 + sqrt(5.0)) / 2.0, tribonacci_constant = 1.83929;
const double regular_distance_factor_const = 5.0, semiregular_distance_factor_const = 5.0;
const int total_number_of_distances = 2;
const std::array<double, 10> distance_sample_face_1_period = { pi / 2.5, pi / 3.0, pi / 3.0, pi / 4.0, pi / 5.0, pi / 6.0, pi / 6.0, pi / 8.0, pi / 8.0, pi / 4.0 };		// the size of the periodic distance sample on the face of the 1st type
const std::array<double, 10> distance_sample_face_2_period = { 0.0, 0.0, 0.0, 0.0, 0.0, pi / 3.0, pi / 4.0, pi / 6.0, pi / 8.0, 0.0 };		// the size of the periodic distance sample on the face of the 1st type
const std::array<double, 10> distance_sample_face_1_frequency = { 1.5, 2.0, 1.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.0 };		// the rotational symmetry of the periodic distance sample of the 1nd type
const std::array<double, 10> distance_sample_face_2_frequency = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 2.0, 2.5, 2.0, 1.0 };		// the rotational symmetry of the periodic distance sample of the 2nd type
const std::array<std::string, 3> distance_type_list = { "analytical", "periodic", "multiple" };
const std::array<std::string, 2> periodic_distance_type_list = { "regular", "semiregular" };
const std::array<std::string, 10> periodic_distance_symmetry_list = { "tetrahedral", "cubic", "octahedral", "dodecahedral", "icosahedral", "truncated_tetrahedral", "cuboctahedral", "icosidodecahedral", "snub_cubic", "diploidal" };
const std::string distance_type_const = distance_type_list[0];	//	analytical, periodic, multiple
const std::string periodic_distance_type_const = periodic_distance_type_list[0];	// regular or semiregular
const std::string periodic_distance_symmetry_const = periodic_distance_symmetry_list[3];	// options: (regular: tetrahedral, cubic, octahedral, dodecahedral, icosahedral), (semiregular: truncated_tetrahedral, cuboctahedral, icosidodecahedral, snub_cubic)

const int edge_radius_value = 8;	//	the final thickness of structure's elements is edge_radius_value
const int rotation_discretization = 1;
const std::string structure_type_const = "closed";	// closed, open

const unsigned int expansion_rate = 4;
const unsigned int expansion_rate_x = 1;
const unsigned int expansion_rate_y = 1;
const unsigned int expansion_rate_z = 1;
const unsigned int structure_expansion_rate_x = 3;
const unsigned int structure_expansion_rate_y = 3;
const unsigned int structure_expansion_rate_z = 3;
const unsigned int structure_edge_expansion_x = edge_radius_value;
const unsigned int structure_edge_expansion_y = edge_radius_value;
const unsigned int structure_edge_expansion_z = edge_radius_value;

const int max_dimension_value = 81;	// minimum number of voxels along one of three axes
const double max_length_value = 1.0;	// minimum length of the computational domain along one of three axes
const std::array<std::string, 7> periodic_lattice_symmetry_group_list = { "triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "cubic", "hexagonal" };
const std::array<std::string, 15> seed_distribution_type_list = { "triclinic", "monoclinic_primitive", "monoclinic_base-centered", "orthorhombic_primitive", "orthorhombic_base-centered", "orthorhombic_body-centered", "orthorhombic_face-centered", "tetragonal_primitive", "tetragonal_body-centered", "trigonal", "hexagonal", "cubic_primitive", "cubic_body-centered", "cubic_face-centered", "diamond_cubic" };
const std::array<std::string, 8> neighborhood_check_type_list = { "by_index", "by_eucledian" };
const std::string seed_distribution_type_const = seed_distribution_type_list[12]; // options: cubic, triangular_prismatic, hexagonal_prismatic, alternated_cubic, bitruncated_cubic, body-centered_cubic, face-centered_cubic, diamond_cubic
const std::string neighborhood_check_type_const = neighborhood_check_type_list[1]; // options: by_index, by_eucledian
const std::string open_cell_structure_computation_type = "per_voxel"; // options: per_voxel, per_neighborhood

const double rt = 1.0 / sqrt(3.0);
const double ro = 1.0 / sqrt(2.0);
const double rd1 = 1.0 / sqrt(3.0), rd2 = (0.5 + sqrt(5.0 / 4.0)) / sqrt(3.0), rd3 = 1.0 / ((0.5 + sqrt(5.0 / 4.0)) * sqrt(3.0));
const double ri1 = 1.0 / sqrt(2.5 + 0.5 * sqrt(5.0)), ri2 = (0.5 + 0.5 * sqrt(5.0)) / sqrt(2.5 + 0.5 * sqrt(5.0));
const std::array<std::array<double, 3>, 4> regular_tetrahedron_vertices = { {{rt, rt, -rt}, {-rt, -rt, -rt}, {rt, -rt, rt}, {-rt, rt, rt}} };
const std::array<std::array<double, 3>, 4> regular_dual_tetrahedron_vertices = { {{rt, -rt, -rt}, {-rt, rt, -rt}, {-rt, -rt, rt}, {rt, rt, rt}} };
const std::array<std::array<double, 3>, 8> regular_hexahedron_vertices = { {{rt, rt, -rt}, {rt, -rt, -rt}, {-rt, rt, -rt}, {-rt, -rt, -rt}, {-rt, -rt, rt}, {rt, -rt, rt}, {-rt, rt, rt}, {rt, rt, rt}} };
const std::array<std::array<double, 3>, 6> regular_octahedron_vertices = { {{0.0, -1.0, 0.0}, {0.0,1.0,0.0}, {-1.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,0.0,-1.0}} };
const std::array<std::array<double, 3>, 20> regular_dodecahedron_vertices = { {{rd1, rd1, -rd1}, {rd1, -rd1, -rd1}, {-rd1, rd1, -rd1}, {-rd1, -rd1, -rd1}, {-rd1, -rd1, rd1}, {rd1, -rd1, rd1}, {-rd1, rd1, rd1}, {rd1, rd1, rd1}, {0.0, -rd3, -rd2}, {0.0, rd3, -rd2}, {0.0, -rd3, rd2}, {0.0, rd3, rd2}, {-rd3, -rd2, 0.0}, {rd3, -rd2, 0.0}, {-rd3, rd2, 0.0}, {rd3, rd2, 0.0}, {-rd2, 0.0, -rd3}, {-rd2, 0.0, rd3}, {rd2, 0.0, -rd3}, {rd2, 0.0, rd3}} };
const std::array<std::array<double, 3>, 12> regular_icosahedron_vertices = { {{ri1, 0.0, -ri2}, {-ri1, 0.0, ri2}, {-ri1, 0.0, -ri2}, {ri1, 0.0, ri2}, {ri2, -ri1, 0.0}, {-ri2, ri1, 0.0}, {-ri2, -ri1, 0.0}, {ri2, ri1, 0.0}, {0.0, ri2, -ri1}, {0.0, -ri2, ri1}, {0.0, -ri2, -ri1}, {0.0, ri2, ri1}} };

const double regular_tetrahedron_face_radius = sqrt(2.0 * (sqrt(3.0) - 1.0) / sqrt(3.0));
const double regular_hexahedron_face_radius = sqrt(2.0 - sqrt(2.0));
const double regular_octahedron_face_radius = sqrt(2.0 * (1.0 - sqrt(2.0 / 3.0)));
const double regular_dodecahedron_face_radius = sqrt(ri1 * ri1 + (1.0 - ri2) * (1.0 - ri2));
const double regular_icosahedron_face_radius = sqrt(rd3 * rd3 + (1.0 - rd2) * (1.0 - rd2));

const double tt1 = 1.0 / sqrt(11.0), tt2 = 3.0 / sqrt(11.0);
const double ttd = 1.0 / sqrt(3.0);
const double co = 1.0 / sqrt(2.0);
const double cd1 = 1.0 / sqrt(2.0), cd2 = 1.0 / sqrt(3.0);
const double id_norm = 0.5 * sqrt(1.0 / (golden_ratio * golden_ratio) + 1.0 + golden_ratio * golden_ratio), id1 = 0.5 / (id_norm * golden_ratio), id2 = 0.5 / id_norm, id3 = 0.5 * golden_ratio / id_norm;
const double rtc_norm1 = sqrt(1.0 + golden_ratio * golden_ratio), rtc_norm2 = sqrt(3.0), rtc_norm3 = sqrt(golden_ratio * golden_ratio + 1.0 / (golden_ratio * golden_ratio)), rtc1 = 1.0 / rtc_norm1, rtc2 = golden_ratio / rtc_norm1, rtc3 = 1.0 / rtc_norm2, rtc4 = golden_ratio / rtc_norm3, rtc5 = 1.0 / (golden_ratio * rtc_norm3);
const double sc_norm = sqrt(1.0 + 1.0 / (tribonacci_constant * tribonacci_constant) + tribonacci_constant * tribonacci_constant), sc1 = 1.0 / sc_norm, sc2 = 1.0 / (tribonacci_constant * sc_norm), sc3 = tribonacci_constant / (sc_norm);
const double pi_norm1 = sqrt(1.0 + (2.0 * tribonacci_constant + 1.0) * (2.0 * tribonacci_constant + 1.0) + pow(tribonacci_constant, 4.0)), pi_norm2 = sqrt(3.0) * tribonacci_constant * tribonacci_constant, pi1 = 1.0 / pi_norm1, pi2 = (2.0 * tribonacci_constant + 1.0) / pi_norm1, pi3 = tribonacci_constant * tribonacci_constant / pi_norm1, pi4 = tribonacci_constant * tribonacci_constant / pi_norm2;
const std::array<std::array<double, 3>, 12> truncated_tetrahedron_vertices = { {{tt1,-tt1,-tt2}, {-tt1,-tt2,tt1}, {-tt2,tt1,-tt1}, {-tt1,tt1,-tt2}, {tt1,-tt2,-tt1}, {-tt2,-tt1,tt1}, {-tt1,-tt1,tt2}, {-tt1,tt2,-tt1}, {tt2,-tt1,-tt1}, {tt1,tt1,tt2}, {tt1,tt2,tt1}, {tt2,tt1,tt1}} };
const std::array<std::array<double, 3>, 8> truncated_tetrahedron_dual_vertices = { {{-ttd,-ttd,-ttd}, {ttd,-ttd,-ttd}, {-ttd,ttd,-ttd}, {-ttd,-ttd,ttd}, {ttd,ttd,-ttd}, {ttd,-ttd,ttd}, {-ttd,ttd,ttd}, {ttd,ttd,ttd}} };
const std::array<std::array<double, 3>, 12> cuboctahedron_vertices = { {{-co,-co,0.0}, {-co,co,0.0}, {co,-co,0.0}, {co,co,0.0}, {-co,0.0,-co}, {-co,0.0,co}, {co,0.0,-co}, {co,0.0,co}, {0.0,-co,-co}, {0.0,-co,co}, {0.0,co,-co}, {0.0,co,co}} };
const std::array<std::array<double, 3>, 14> rhombic_dodecahedron_vertices = { {{0.0,0.0,-1.0}, {0.0,0.0,1.0}, {1.0,0.0,0.0}, {-1.0,0.0,0.0}, {0.0,-1.0,0.0}, {0.0,1.0,0.0}, {-cd2,-cd2,-cd2}, {-cd2,-cd2,cd2}, {-cd2,cd2,-cd2}, {cd2,-cd2,-cd2}, {-cd2,cd2,cd2}, {cd2,-cd2,cd2}, {cd2,cd2,-cd2}, {cd2,cd2,cd2}} };
const std::array<std::array<double, 3>, 30> icosidodecahedron_vertices = { {{1.0,0.0,0.0}, {-1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,-1.0,0.0}, {0.0,0.0,1.0}, {0.0,0.0,-1.0}, {id1,id2,id3}, {-id1,id2,id3}, {id1,-id2,id3}, {id1,id2,-id3}, {-id1,-id2,id3}, {-id1,id2,-id3}, {id1,-id2,-id3}, {-id1,-id2,-id3}, {id2,id3,id1}, {-id2,id3,id1}, {id2,-id3,id1}, {id2,id3,-id1}, {-id2,-id3,id1}, {-id2,id3,-id1}, {id2,-id3,-id1}, {-id2,-id3,-id1}, {id3,id1,id2}, {-id3,-id1,id2}, {-id3,id1,-id2}, {id3,-id1,-id2}, {-id3,-id1,-id2}, {-id3,id1,id2}, {id3,-id1,id2}, {id3,id1,-id2}} };
const std::array<std::array<double, 3>, 32> rhombic_triacontahedron_vertices = { {{rtc1,0.0,rtc2}, {-rtc1,0.0,rtc2}, {rtc1,0.0,-rtc2}, {-rtc1,0.0,-rtc2}, {rtc2,rtc1,0.0}, {-rtc2,rtc1,0.0}, {rtc2,-rtc1,0.0}, {-rtc2,-rtc1,0.0}, {0.0,rtc2,rtc1}, {0.0,-rtc2,rtc1}, {0.0,rtc2,-rtc1}, {0.0,-rtc2,-rtc1}, {-rtc3,-rtc3,-rtc3}, {rtc3,-rtc3,-rtc3}, {-rtc3,rtc3,-rtc3}, {-rtc3,-rtc3,rtc3}, {rtc3,rtc3,-rtc3}, {rtc3,-rtc3,rtc3}, {-rtc3,rtc3,rtc3}, {rtc3,rtc3,rtc3}, {rtc4,0.0,rtc5}, {-rtc4,0.0,rtc5}, {rtc4,0.0,-rtc5}, {-rtc4,0.0,-rtc5}, {rtc5,rtc4,0.0}, {-rtc5,rtc4,0.0}, {rtc5,-rtc4,0.0}, {-rtc5,-rtc4,0.0}, {0.0,-rtc5,-rtc4}, {0.0,rtc5,-rtc4}, {0.0,-rtc5,rtc4}, {0.0,rtc5,rtc4}} };
const std::array<std::array<double, 3>, 24> snub_cube_vertices = { {{-sc1,-sc2,-sc3}, {-sc1,sc2,sc3}, {sc1,-sc2,sc3}, {sc1,sc2,-sc3}, {-sc2,-sc3,-sc1}, {-sc2,sc3,sc1}, {sc2,-sc3,sc1}, {sc2,sc3,-sc1}, {-sc3,-sc1,-sc2}, {-sc3,sc1,sc2}, {sc3,-sc1,sc2}, {sc3,sc1,-sc2}, {-sc2,-sc1,sc3}, {sc2,-sc1,-sc3}, {-sc2,sc1,-sc3}, {sc2,sc1,sc3}, {-sc1,-sc3,sc2}, {sc1,-sc3,-sc2}, {-sc1,sc3,-sc2}, {sc1,sc3,sc2}, {-sc3,-sc2,sc1}, {sc3,-sc2,-sc1}, {-sc3,sc2,-sc1}, {sc3,sc2,sc1}} };
const std::array<std::array<double, 3>, 38> pentagonal_icositetrahedron_vertices { {{-pi1,-pi2,pi3}, {pi1,-pi2,-pi3}, {-pi1,pi2,-pi3}, {pi1,pi2,pi3}, {-pi2,-pi3,pi1}, {pi2,-pi3,-pi1}, {-pi2,pi3,-pi1}, { pi2,pi3,pi1 }, {-pi3,-pi1,pi2}, {pi3,-pi1,-pi2}, {-pi3,pi1,-pi2}, {pi3,pi1,pi2}, {-pi2,pi1,pi3}, {pi2,-pi1,pi3}, {pi2,pi1,-pi3}, {-pi2,-pi1,-pi3}, {-pi1,pi3,pi2}, {pi1,-pi3,pi2}, {pi1,pi3,-pi2}, {-pi1,-pi3,-pi2}, {-pi3,-pi2,-pi1}, {-pi3,pi2,pi1}, {pi3,-pi2,pi1}, {pi3,pi2,-pi1}, {-1.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,-1.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,-1.0}, {0.0,0.0,1.0}, {-pi4,-pi4,-pi4}, {pi4,-pi4,-pi4}, {-pi4,-pi4,pi4}, {-pi4,pi4,-pi4}, {-pi4,pi4,pi4}, {pi4,-pi4,pi4}, {pi4,pi4,-pi4}, {pi4,pi4,pi4}} };

const double truncated_tetrahedron_face_radius_1 = sqrt(2.0 - (10.0 / 9.0) * sqrt(3.0)), truncated_tetrahedron_face_radius_2 = sqrt(2.0 - (2.0 / 3.0) * sqrt(3.0));
const double truncated_tetrahedron_distance_factor = -(3.0 / 3.0) * (8.0 / 11.0);
const double cuboctahedron_face_radius_1 = sqrt(2.0 - (4.0 / 3.0) * sqrt(2.0)), cuboctahedron_face_radius_2 = sqrt(2.0 - (2.0 / 3.0) * sqrt(6.0));
const double icosidodecahedron_face_radius_1 = 0.1884306430661105, icosidodecahedron_face_radius_2 = 0.45950584109472237;

const std::array<double, 8> truncated_tetrahedron_face_types_radius = { truncated_tetrahedron_face_radius_2, truncated_tetrahedron_face_radius_1, truncated_tetrahedron_face_radius_1, truncated_tetrahedron_face_radius_1, truncated_tetrahedron_face_radius_2, truncated_tetrahedron_face_radius_2, truncated_tetrahedron_face_radius_2, truncated_tetrahedron_face_radius_1 };
const std::array<double, 14> cuboctahedron_face_types_radius = { cuboctahedron_face_radius_2, cuboctahedron_face_radius_2, cuboctahedron_face_radius_2, cuboctahedron_face_radius_2, cuboctahedron_face_radius_2, cuboctahedron_face_radius_2, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1, cuboctahedron_face_radius_1 };
const std::array<double, 32> icosidodecahedron_face_types_radius = { icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_2,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1,icosidodecahedron_face_radius_1 };

const double truncated_tetrahedron_compensation = 0.6963106238227914;
const double cuboctahedron_compensation = 0.21877959948235687;
const double icosidodecahedron_compensation = 0.1670431012213515;
const double snub_cube_compensation = 0.10529166910111312;

const std::array<double, 8> truncated_tetrahedron_face_types_distance_compensation = { 0.0, truncated_tetrahedron_compensation, truncated_tetrahedron_compensation, truncated_tetrahedron_compensation, 0.0, 0.0, 0.0, truncated_tetrahedron_compensation };
const std::array<double, 14> cuboctahedron_face_types_distance_compensation = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation, cuboctahedron_compensation };
const std::array<double, 32> icosidodecahedron_face_types_distance_compensation = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation, icosidodecahedron_compensation };
const std::array<double, 38> snub_cube_face_types_distance_compensation = { snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation, snub_cube_compensation };

//const std::array<double, 3> array_data = {1.0,2.0,3.0};