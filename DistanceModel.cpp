#include <math.h>
#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <fstream>
#include <omp.h>
#include "GridComputationMethods.h"
#include "parameters.h"
#include "MetaDataWrite.h"






std::array<double, 3> CartesianToSpherical(double x, double y, double z)            //  result[0] = radius, result[1] = phi_angle, result[2] = theta_angle;
{
    std::array<double, 3> spherical_point;

    spherical_point[0] = sqrt(x * x + y * y + z * z);
    if (y >= 0.0) { spherical_point[1] = atan2(y, x); }
    else { spherical_point[1] = 2.0 * pi + atan2(y, x); }
    if (z != 0.0)
    {
        if (z > 0)
        {
            spherical_point[2] = atan(sqrt(x * x + y * y) / z);
        }
        else
        {
            spherical_point[2] = pi + atan(sqrt(x * x + y * y) / z);
        }
    }
    else
    {
        spherical_point[2] = 0.5 * pi;
    }

    return(spherical_point);
}

std::vector<std::array<double, 3>> vertices_initialization(std::string polyhedron_family, std::string polyhedron_type)
{
    std::vector<std::array<double, 3>> result;
    std::array<double, 3> vertice;

    if (polyhedron_family == "regular")
    {
        if (polyhedron_type == "tetrahedral")
        {
            vertice = { 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
        }
        if (polyhedron_type == "cubic")
        {
            vertice = { 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
            result.push_back(vertice);
        }
        if (polyhedron_type == "octahedral")
        {
            vertice = { 1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0 };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0 };
            result.push_back(vertice);
            vertice = { 1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0 };
            result.push_back(vertice);
            vertice = { -1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0 };
            result.push_back(vertice);
            vertice = { 0.0, 0.0, 1.0 };
            result.push_back(vertice);
            vertice = { 0.0, 0.0, -1.0 };
            result.push_back(vertice);
        }
        if (polyhedron_type == "dodecahedral")
        {

        }
    }
    else
    {
        if (polyhedron_family == "semiregular")
        {
            if (polyhedron_type == "truncated_tetrahedral")
            {

            }
        }
        else { std::cout << "ERROR: incorrect polyhedron family or type name." << std::endl; }
    }

    return(result);
}

double SphericalHarmonicsReal(double phi, double theta, int degree, int order)          //  Y_degree^order(phi_angle, theta_angle), order <= degree;
{
    return(std::sph_legendre(degree, order, theta) * std::cos(order * phi));
}

std::vector<double> InterpolationParameter(double x, double y, double z, const int indicator)
{
    std::vector<double> result;
    if (indicator == 2)
    {
        /*if ((x >= 5) && (x <= 15.0))
        {
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            result.push_back(0.0);
            result.push_back(1.0);
        }*/

        /*if (x <= 10.0)
        {
            result.push_back(x / 10.0);
            result.push_back(1.0 - (x / 10.0));
        }
        else
        {
                    result.push_back(1.0 - ((x - 10.0) / 10.0));
                    result.push_back((x - 10.0) / 10.0);
        }*/

        /*std::array<double, 3> spherical_coord;
        std::array<double, 3> radius_vector;
        radius_vector[0] = x - 0.5;
        radius_vector[1] = y - 0.5;
        radius_vector[2] = z - 0.5;
        spherical_coord = CartesianToSpherical(radius_vector[0], radius_vector[1], radius_vector[2]);
        if (spherical_coord[0] <= (0.5 * (0.9 * fabs(cos(2.0 * spherical_coord[2])) + 0.1)))
        {
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            result.push_back(0.0);
            result.push_back(1.0);
        }*/

        //complex BCC with two seed types
        /*if ((x == 0.5) && (y = 0.5) && (z == 0.5))
        {
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            result.push_back(0.0);
            result.push_back(1.0);
        }*/

        //complex tetragonal lattice with two seed types
        /*if ((y == 1.0) || ((x == 0.5) && (y == 1.5)) || ((x == 1.5) && (y == 0.5)))
        {
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            result.push_back(0.0);
            result.push_back(1.0);
        }*/

        //complex hexagonal lattice with two seed types
        if (((y == 0.0) && ((x == 2.0) || (x == 8.0))) || ((y != 0.0) && ((x == 5.0) || (x == 11.0))))
        {
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            result.push_back(0.0);
            result.push_back(1.0);
        }
    }
    if (indicator == 3)
    {
        /*if (((x - 4.0) * (x - 4.0) + (y - 4.0) * (y - 4.0)) <= 6.0)
        {
            result.push_back(0.0);
            result.push_back(1.0);
            result.push_back(0.0);
        }
        else
        {
            if (((x < 4.0) && (y < 4.0)) || ((x >= 4.0) && (y >= 4.0)))
            {
                result.push_back(0.0);
                result.push_back(0.0);
                result.push_back(1.0);
            }
            else
            {
                result.push_back(1.0);
                result.push_back(0.0);
                result.push_back(0.0);
            }
        }*/

        /*if ((x <= 3.5) || (x >= 17.5))
        {
            result.push_back(1.0);
            result.push_back(0.0);
            result.push_back(0.0);
        }
        else
        {
            if (x <= 10.5)
            {
                result.push_back(0.0);
                result.push_back(1.0);
                result.push_back(0.0);
            }
            else
            {
                result.push_back(0.0);
                result.push_back(0.0);
                result.push_back(1.0);
            }
        }*/

        if ((z == 1.0) || (z == 3.0))
        {
            result.push_back(1.0);
            result.push_back(0.0);
            result.push_back(0.0);
        }
        else
        {
            if (((z == 0.0) && (((x == 0.0) && (y == 0.0)) || ((x == 1.25) && (y == 0.5)))) || ((z == 2.0) && (((x == 0.0) && (y == 0.5)) || ((x == 1.25) && (y == 0.0)))))
            {
                result.push_back(0.0);
                result.push_back(1.0);
                result.push_back(0.0);
            }
            else
            {
                result.push_back(0.0);
                result.push_back(0.0);
                result.push_back(1.0);
            }
        }
    }
    if (indicator == 4)
    {
        if ((x <= 4.0) || (x > 16.0))
        {
            result.push_back(1.0);
            result.push_back(0.0);
            result.push_back(0.0);
            result.push_back(0.0);
        }
        else
        {
            if (x <= 8.0)
            {
                result.push_back(0.0);
                result.push_back(1.0);
                result.push_back(0.0);
                result.push_back(0.0);
            }
            else
            {
                if (x <= 12.0)
                {
                    result.push_back(0.0);
                    result.push_back(0.0);
                    result.push_back(1.0);
                    result.push_back(0.0);
                }
                else
                {
                    result.push_back(0.0);
                    result.push_back(0.0);
                    result.push_back(0.0);
                    result.push_back(1.0);
                }
            }
        }
    }

    return(result);
}

double DistanceLinearInterpolation(std::vector<double> function, std::vector<double> a)
{
    double check = 0.0;
    for (int i = 0; i < a.size(); i++) { check = check + a[i]; }
    if ((check < 0.0) || (check > 1.0)) { std::cout << "ERROR: FAULSE DISTANCE INTERPOLATION VALUE.\n" << std::endl; }
    if (function.size() != a.size()) { std::cout << "ERROR: FAULSE DISTANCE INTERPOLATION PARAMETERS.\n" << std::endl; }
    double result = 0.0;
    for (int i = 0; i < a.size(); i++) { result = result + function[i] * a[i]; }
    return(result);
}

double SquaredEucledianDistance(double x0, double y0, double z0, double x1, double y1, double z1)
{
    return((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
}

double EucledianDistance(double x0, double y0, double z0, double x1, double y1, double z1)
{
    return(sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1)));
}

double RegularDistanceSampleFunction(double phi, double theta, double regular_distance_factor, unsigned int period_index)      //  phi_angle varies from 0 to 2.0 * pi, theta_angle varies from 0 to 0.25 * pi;
{
    double result = 1.0, period = distance_sample_face_1_period[period_index];

    if (theta <= period)
    { 
        result = result + regular_distance_factor * std::max(1.3 * std::exp(-(100.0 * theta * theta)), fabs(std::cos(0.5 * (pi / period) * theta) * std::cos(distance_sample_face_1_frequency[period_index] * phi)));
        //result = result + regular_distance_factor * fabs(std::sin((pi / period) * theta)) * std::sin(distance_sample_face_1_frequency[period_index] * phi);
    }

    return(result);
}

double SemiregularDistanceSampleFunction(double phi, double theta, unsigned int face_type, double semiregular_distance_factor, unsigned int period_index)
{
    double result = 1.0, period1 = distance_sample_face_1_period[period_index], period2 = distance_sample_face_2_period[period_index];

    if (face_type == 1)
    {
        if (theta <= period1)
        {
            result = result + semiregular_distance_factor * std::max(fabs(std::cos(0.5 * (pi / period1) * theta) * std::cos(distance_sample_face_1_frequency[period_index] * phi)), 1.3 * std::exp(-(100.0 * theta * theta)));
            //result = result + semiregular_distance_factor * std::cos(0.5 * (pi / period1) * theta);
        }
    }
    if (face_type == 2)
    {
        if (theta <= period2)
        {
            result = result + std::max(fabs(std::cos(0.5 * (pi / period2) * theta) * std::cos(distance_sample_face_2_frequency[period_index] * phi)), 1.3 * std::exp(-(100.0 * theta * theta)));
            //result = result + fabs(std::cos(0.5 * (pi / period2) * theta));
        }
    }

    return(result);
}

std::array<double, 3> PointRotationCartesian(double x, double y, double z, double angle, const std::string& axis)
{
    std::array<double, 3> rotated_point = {x, y, z};
    if (axis == "X")
    {
        rotated_point[0] = x;
        rotated_point[1] = cos(angle) * y - sin(angle) * z;
        rotated_point[2] = sin(angle) * y + cos(angle) * z;
    }
    if (axis == "Y")
    {
        rotated_point[0] = cos(angle) * x + sin(angle) * z;
        rotated_point[1] = y;
        rotated_point[2] = - sin(angle) * x + cos(angle) * z;
    }
    if (axis == "Z")
    {
        rotated_point[0] = cos(angle) * x - sin(angle) * y;
        rotated_point[1] = sin(angle) * x + cos(angle) * y;
        rotated_point[2] = z;
    }
    return(rotated_point);
}

//  As one need to rotate a whole shape, at first the later is rotated by a given phi and then by a prescribed theta
std::array<double, 3> PointRotationSpherical(double x, double y, double z, double phi, double theta)

{
    std::array<double, 3> rotated_point, temp_point;
    //double phi_current;
    //if (y >= 0.0) { phi_current = atan2(y, x); }
    //else { phi_current = 2.0 * pi + atan2(y, x); }

    //temp_point = PointRotationCartesian(x, y, z, (-phi_current), "Z");
    temp_point = PointRotationCartesian(x, y, z, phi, "Z");

    rotated_point = PointRotationCartesian(temp_point[0], temp_point[1], temp_point[2], theta, "Y");
    
    return(rotated_point);
}

double StarShapedAnalyticalDistance(double phi, double theta)
{
    double d = 1.0;

    //return (0.25*sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist) / (4.0 + sin(6.0*theta)*cos(2.0*phi)));
    //return (sqrt(x_dist*x_dist / 10.0 + y_dist*y_dist / 10.0 + z_dist*z_dist));
    //return (sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist) / (1.0 + 2.5*fabs(sin(2.0*theta)*sin(2.0*phi))));
    //return (sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist));

    /*if (sin(2.0 * phi) * sin(theta) * sin(theta) * cos(theta) < 0.0)
    {
        d = 0.2;
    }
    else
    {
        d = 0.2 + fabs(1.0 / 4.0 * sqrt(105.0 / (2.0 * pi)) * sin(2.0 * phi) * sin(theta) * sin(theta) * cos(theta));
    }*/

    //d = 1.0 + SphericalHarmonicsReal(phi, theta, 7, 5);

    //d = 1.0 + 0.8 * (SphericalHarmonicsReal(phi + 0.5 * pi, theta, 2, 1) + SphericalHarmonicsReal(phi, theta, 3, 1) + SphericalHarmonicsReal(phi + 0.35 * pi, theta, 2, 2));  // triclinic distance (1 point group)

    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 4, 1) + SphericalHarmonicsReal(phi, theta, 6, 4));  // monoclinic distance (2/m point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 2, 1) + SphericalHarmonicsReal(phi - 0.25 * pi, theta, 3, 2) + SphericalHarmonicsReal(phi, theta, 4, 3)); // monoclinic distance (2 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 5, 2) + SphericalHarmonicsReal(phi, theta, 4, 1) + SphericalHarmonicsReal(phi, theta, 3, 1)); // monoclinic distance (m point group)

    //d = 1.0 + (SphericalHarmonicsReal(phi - 0.25 * pi, theta, 3, 2) + SphericalHarmonicsReal(phi - 0.375 * pi, theta, 5, 4));    // orthorhombic distance (222 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 3, 2) + SphericalHarmonicsReal(phi, theta, 5, 4));    // orthorhombic distance (mm2 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 4, 2) + SphericalHarmonicsReal(phi, theta, 6, 4));    // orthorhombic distance (mmm point group)

    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 5, 4) + SphericalHarmonicsReal(phi, theta, 6, 4));    // tetragonal (4mm point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi - 0.12 * pi, theta, 4, 4) + SphericalHarmonicsReal(phi, theta, 5, 4) + SphericalHarmonicsReal(phi, theta, 6, 4));    // tetragonal (4 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 8, 4) + SphericalHarmonicsReal(phi, theta, 6, 4));    // tetragonal (4/mmm point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi - 0.12 * pi, theta, 4, 4) + SphericalHarmonicsReal(phi, theta, 8, 4) + SphericalHarmonicsReal(phi, theta, 6, 4));    // tetragonal (4/m point group)
    d = 1.0 + (SphericalHarmonicsReal(phi - 0.125 * pi, theta, 5, 4) + SphericalHarmonicsReal(phi - 0.250 * pi, theta, 6, 4));    // tetragonal (422 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 3, 2) + SphericalHarmonicsReal(phi, theta, 4, 4));    // tetragonal (-4m2 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 3, 2) + SphericalHarmonicsReal(phi, theta, 4, 4) + SphericalHarmonicsReal(phi - 0.08 * pi, theta, 6, 4));    // tetragonal (-4 point group)

    //d = 1.3 + (SphericalHarmonicsReal(phi, theta, 4, 3) + SphericalHarmonicsReal(phi, theta, 5, 3) + SphericalHarmonicsReal(phi - 0.08 * pi, theta, 6, 3)); // trigonal (3 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 4, 3) + SphericalHarmonicsReal(phi, theta, 6, 6) + SphericalHarmonicsReal(phi - 0.08 * pi, theta, 8, 6)); // trigonal (-3 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi - pi / 6.0, theta, 4, 3) + SphericalHarmonicsReal(phi, theta, 5, 3) + SphericalHarmonicsReal(phi - pi / 6.0, theta, 6, 3)); // trigonal (32 point group)
    //d = 1.3 + (SphericalHarmonicsReal(phi, theta, 5, 3) + SphericalHarmonicsReal(phi, theta, 7, 6) + SphericalHarmonicsReal(phi, theta, 6, 3)); // trigonal (3m point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 4, 3) + SphericalHarmonicsReal(phi, theta, 6, 6)); // trigonal (-3m point group)

    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 7, 6) + SphericalHarmonicsReal(phi, theta, 8, 6) + SphericalHarmonicsReal(phi - 0.1 * pi, theta, 9, 6));    // hexagonal (6 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 5, 3) + SphericalHarmonicsReal(phi, theta, 7, 3) + SphericalHarmonicsReal(phi - 0.1 * pi, theta, 8, 6));    // hexagonal (-6 point group)    
    //d = 1.1 + (SphericalHarmonicsReal(phi - 0.08 * pi, theta, 6, 6) + SphericalHarmonicsReal(phi, theta, 8, 6) + SphericalHarmonicsReal(phi - 0.15 * pi, theta, 10, 6));    // hexagonal (6/m point group)
    //d = 1.2 + (SphericalHarmonicsReal(phi - pi / 12.0, theta, 7, 6) + SphericalHarmonicsReal(phi, theta, 6, 6) + SphericalHarmonicsReal(phi, theta, 8, 6) + SphericalHarmonicsReal(phi - pi / 12.0, theta, 11, 6));    // hexagonal (622 point group)
    //d = 1.3 + (SphericalHarmonicsReal(phi, theta, 6, 6) + SphericalHarmonicsReal(phi, theta, 7, 6) + SphericalHarmonicsReal(phi, theta, 8, 6) + SphericalHarmonicsReal(phi, theta, 10, 6));    // hexagonal (6mm point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 5, 3) + SphericalHarmonicsReal(phi, theta, 10, 6));    // hexagonal (-6m2 point group)
    //d = 1.0 + (SphericalHarmonicsReal(phi, theta, 6, 6) + SphericalHarmonicsReal(phi, theta, 8, 6) + SphericalHarmonicsReal(phi, theta, 10, 6));    // hexagonal (6/mmm point group)

    //d = 0.7 * d1 + 0.3 * d2;

    return d;
}

double RegularDistance(double x, double y, double z, double regular_distance_factor, const std::string& periodic_distance_symmetry)
{
    unsigned int order_phi = 1, quadrant_phi = 0, period_index = 0;
    double result = 1.0, phi_step = 1.0, theta_step = 1.0, phi_central = 0.0, norm = sqrt(x * x + y * y + z * z), closest_distance = 100.0;
    std::array<double, 3> spherical_coordinates = CartesianToSpherical(x, y, z), face_center_closest, cartesian_coordinates, spherical_coordinates_face;

    if (periodic_distance_symmetry == "tetrahedral")
    {
        //#pragma omp parallel for
        for (const auto& point : regular_dual_tetrahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (spherical_coordinates_face[2] > 0.5 * pi)
        {
            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    if (periodic_distance_symmetry == "cubic")
    {
        period_index = 1;
        //#pragma omp parallel for
        for (const auto& point : regular_octahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (spherical_coordinates_face[2] > 0.5 * pi)
        {
            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.75 * pi, "Z");
        }
        else
        {
            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.25 * pi, "Z");
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    if (periodic_distance_symmetry == "octahedral")
    {
        period_index = 2;
        //#pragma omp parallel for
        for (const auto& point : regular_hexahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (spherical_coordinates_face[2] < 0.5 * pi)
        {
            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    if (periodic_distance_symmetry == "dodecahedral")
    {
        period_index = 3;
        //#pragma omp parallel for
        for (const auto& point : regular_icosahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (face_center_closest[2] == 0.0)
        {
            if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 1.1 * pi, "Z");
            }
            else
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.9 * pi, "Z");
            }
        }
        else
        {
            if (spherical_coordinates_face[2] < 0.5 * pi)
            {
                if (face_center_closest[2] == ri1)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
                }
            }
            else
            {
                if (face_center_closest[2] == -ri2)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
                }
            }
        }
        

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    if (periodic_distance_symmetry == "icosahedral")
    {
        period_index = 4;
        //#pragma omp parallel for
        for (const auto& point : regular_dodecahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (face_center_closest[2] == 0.0)
        {
            if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -pi / 6.0, "Z");
            }
            else
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi / 6.0, "Z");
            }
        }
        else
        {
            if ((face_center_closest[2] != rd2) && (face_center_closest[2] != -rd3))
            {
                if (spherical_coordinates_face[2] < 0.5 * pi)
                {
                    if (face_center_closest[2] == rd3)
                    {
                        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
                    }
                    else
                    {
                        if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
                        {
                            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi / 10.0, "Z");
                        }
                        else
                        {
                            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -pi / 10.0, "Z");
                        }
                    }
                }
                else
                {
                    if (face_center_closest[2] == -rd2)
                    {
                        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
                    }
                    else
                    {
                        if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
                        {
                            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 9.0 / 10.0 * pi, "Z");
                        }
                        else
                        {
                            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 11.0 / 10.0 * pi, "Z");
                        }
                    }
                }
            }
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    if (periodic_distance_symmetry == "diploidal")
    {
        period_index = 9;
        //#pragma omp parallel for
        for (const auto& point : regular_octahedron_vertices)
        {
            if (EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) < closest_distance)
            {
                face_center_closest = point;
                closest_distance = EucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm);
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if ((face_center_closest[1] == 0.0) && (face_center_closest[2] == 0.0))
        {
            cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.5 * pi, "Z");
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = RegularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], regular_distance_factor, period_index);
    }

    return(result);
}

double SemiRegularDistance(double x, double y, double z, double semiregular_distance_factor, const std::string& periodic_distance_symmetry)
{
    unsigned int order_phi = 1, quadrant_phi = 0, face_type = 0, period_index = 0;
    double result = 1.0, phi_step = 1.0, theta_step = 1.0, phi_central = 0.0, norm = sqrt(x * x + y * y + z * z), closest_distance = 100.0, check_radius = 0.0;
    std::array<double, 3> spherical_coordinates = CartesianToSpherical(x, y, z), face_center_closest, cartesian_coordinates, spherical_coordinates_face, point;

    if (periodic_distance_symmetry == "truncated_tetrahedral")
    {
        period_index = 5;
        //#pragma omp parallel for
        for (unsigned int index = 0; index < truncated_tetrahedron_dual_vertices.size(); index++)
        {
            point = truncated_tetrahedron_dual_vertices[index];
            
            if ((SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + truncated_tetrahedron_face_types_distance_compensation[index]) < closest_distance)
            {
                if (((point[0] * point[1] > 0.0) && (point[2] > 0.0)) || ((point[0] * point[1] < 0.0) && (point[2] < 0.0)))
                {
                    face_type = 1;
                }
                else
                {
                    face_type = 2;
                }
                face_center_closest = point;
                closest_distance = SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + truncated_tetrahedron_face_types_distance_compensation[index];
            }
        }
        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (face_type == 1)
        {
            if (spherical_coordinates_face[2] < 0.5 * pi)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
        }
        else
        {
            if (spherical_coordinates_face[2] < 0.5 * pi)
            {
                if (face_center_closest[0] > 0.0)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -(7.0 / 6.0) * pi, "Z");
                }
                else
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], (5.0 / 6.0) * pi, "Z");
                }
            }
            else
            {
                if (face_center_closest[0] > 0.0)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -(1.0 / 6.0) * pi, "Z");
                }
                else
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -(1.0 / 6.0) * pi, "Z");
                }
            }
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = SemiregularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], face_type, semiregular_distance_factor, period_index);
    }
    if (periodic_distance_symmetry == "cuboctahedral")
    {
        period_index = 6;
        //#pragma omp parallel for
        for (unsigned int index = 0; index < rhombic_dodecahedron_vertices.size(); index++)
        {
            point = rhombic_dodecahedron_vertices[index];
            if ((SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + cuboctahedron_face_types_distance_compensation[index]) < closest_distance)
            {
                if (cuboctahedron_face_types_radius[index] == cuboctahedron_face_radius_1) { face_type = 1; }
                else { face_type = 2; }
                face_center_closest = point;
                closest_distance = SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + cuboctahedron_face_types_distance_compensation[index];
            }
        }

        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (face_type == 1)
        {
            if (spherical_coordinates_face[2] > 0.5 * pi)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
        }
        else
        {
            if (face_center_closest[2] == -1.0)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = SemiregularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], face_type, semiregular_distance_factor, period_index);
    }
    if (periodic_distance_symmetry == "icosidodecahedral")
    {
        period_index = 7;
        //#pragma omp parallel for
        for (unsigned int index = 0; index < rhombic_triacontahedron_vertices.size(); index++)
        {
            point = rhombic_triacontahedron_vertices[index];
            if ((SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + icosidodecahedron_face_types_distance_compensation[index]) < closest_distance)
            {
                if (icosidodecahedron_face_types_radius[index] == icosidodecahedron_face_radius_1) { face_type = 1; }
                else { face_type = 2; }
                face_center_closest = point;
                closest_distance = SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + icosidodecahedron_face_types_distance_compensation[index];
            }
        }

        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        if (face_type == 1)
        {
            if ((face_center_closest[2] == rtc4) || (face_center_closest[2] == -rtc5))
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
            if (fabs(face_center_closest[2]) == rtc3)
            {
                if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.5 * pi, "Z");
                }
                else
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -0.5 * pi, "Z");
                }
            }
            if (face_center_closest[2] == 0.0)
            {
                if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi / 6.0, "Z");
                }
                else
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -pi / 6.0, "Z");
                }
            }
        }
        else
        {
            if (face_center_closest[2] == rtc2)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
            if (face_center_closest[2] == -rtc1)
            {
                cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], pi, "Z");
            }
            if (face_center_closest[2] == 0.0)
            {
                if ((face_center_closest[0] * face_center_closest[1]) > 0.0)
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], 0.09999991553 * pi, "Z");
                }
                else
                {
                    cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -0.09999991553 * pi, "Z");
                }
            }
        }

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = SemiregularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], face_type, semiregular_distance_factor, period_index);
    }
    if (periodic_distance_symmetry == "snub_cubic")
    {
        period_index = 8;
        //#pragma omp parallel for
        for (unsigned int index = 0; index < pentagonal_icositetrahedron_vertices.size(); index++)
        {
            point = pentagonal_icositetrahedron_vertices[index];
            if ((SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + snub_cube_face_types_distance_compensation[index]) < closest_distance)
            {
                if (snub_cube_face_types_distance_compensation[index] == snub_cube_compensation) { face_type = 1; }
                else { face_type = 2; }
                face_center_closest = point;
                closest_distance = SquaredEucledianDistance(point[0], point[1], point[2], x / norm, y / norm, z / norm) + snub_cube_face_types_distance_compensation[index];
            }
        }

        spherical_coordinates_face = CartesianToSpherical(face_center_closest[0], face_center_closest[1], face_center_closest[2]);

        cartesian_coordinates = PointRotationCartesian(x, y, z, -spherical_coordinates_face[1], "Z");
        cartesian_coordinates = PointRotationCartesian(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2], -spherical_coordinates_face[2], "Y");

        spherical_coordinates = CartesianToSpherical(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2]);         //  Now theta_angle can vary within [0, 0.25 * pi] and phi angle lies within [0, 2.0 * pi];
        result = SemiregularDistanceSampleFunction(spherical_coordinates[1], spherical_coordinates[2], face_type, semiregular_distance_factor, period_index);
    }
    return(result);
}

double StarShapedPeriodicDistance(double x, double y, double z, double regular_distance_factor, double semiregular_distance_factor, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    double result = 1.0;
    if (periodic_distance_type == "regular")
    {
        result = RegularDistance(x, y, z, regular_distance_factor, periodic_distance_symmetry);
    }
    if (periodic_distance_type == "semiregular")
    {
        result = SemiRegularDistance(x, y, z, semiregular_distance_factor, periodic_distance_symmetry);
    }
    return(result);
}

std::vector<double> MultipleDistanceFunction(double x, double y, double z, const int indicator)
{
    std::array<double, 3> point_spherical = CartesianToSpherical(x, y, z);
    std::vector<double> result;

    if (indicator == 2)
    {
        //result = 1.0 + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 0, 0);
        //tetragonal
        //result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 2) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 4, 4) + SphericalHarmonicsReal(point_spherical[1] - 0.08 * pi, point_spherical[2], 6, 4)));
        //orthorhombic
        //result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 4, 2) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 6, 4)));
        //hexagonal
        result.push_back(1.1 + (SphericalHarmonicsReal(point_spherical[1] - 0.08 * pi, point_spherical[2], 6, 6) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 8, 6) + SphericalHarmonicsReal(point_spherical[1] - 0.15 * pi, point_spherical[2], 10, 6)));
        //trigonal
        result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1] - pi / 6.0, point_spherical[2], 4, 3) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 5, 3) + SphericalHarmonicsReal(point_spherical[1] - pi / 6.0, point_spherical[2], 6, 3)));
    }
    else
    {
        if (indicator == 3)
        {
            /*result.push_back(StarShapedPeriodicDistance(x, y, z, 1.0, 1.0, periodic_distance_type_list[0], periodic_distance_symmetry_list[0]));
            result.push_back(StarShapedPeriodicDistance(x, y, z, 1.0, 1.0, periodic_distance_type_list[0], periodic_distance_symmetry_list[3]));
            std::array<double, 3> point = { x, y, z };
            point = PointRotationSpherical(x, y, z, pi / 6.0, pi / 6.0);
            result.push_back(StarShapedPeriodicDistance(point[0], point[1], point[2], 1.0, 1.0, periodic_distance_type_list[1], periodic_distance_symmetry_list[7]));*/

            //result.push_back(1.0 + 0.8 * (SphericalHarmonicsReal(point_spherical[1] + 0.5 * pi, point_spherical[2], 2, 1) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 1) + SphericalHarmonicsReal(point_spherical[1] + 0.35 * pi, point_spherical[2], 2, 2)));
            //result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 2) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 4, 4) + SphericalHarmonicsReal(point_spherical[1] - 0.08 * pi, point_spherical[2], 6, 4)));
            //result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 5, 4) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 6, 4) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 1)));

            //orthorhombic
            result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1] - 0.25 * pi, point_spherical[2], 3, 2) + SphericalHarmonicsReal(point_spherical[1] - 0.375 * pi, point_spherical[2], 5, 4)));
            //tetragonal
            result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1] - 0.125 * pi, point_spherical[2], 5, 4) + SphericalHarmonicsReal(point_spherical[1] - 0.250 * pi, point_spherical[2], 6, 4)));
            //euclidean
            result.push_back(1.0);
        }
        else
        {
            if (indicator == 4)
            {
                //triclinic
                result.push_back(1.0 + 0.8 * (SphericalHarmonicsReal(point_spherical[1] + 0.5 * pi, point_spherical[2], 2, 1) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 1) + SphericalHarmonicsReal(point_spherical[1] + 0.35 * pi, point_spherical[2], 2, 2)));
                //monoclinic
                result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 4, 1) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 6, 4)));
                //tetragonal
                result.push_back(1.0 + (SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 3, 2) + SphericalHarmonicsReal(point_spherical[1], point_spherical[2], 4, 4) + SphericalHarmonicsReal(point_spherical[1] - 0.08 * pi, point_spherical[2], 6, 4)));
                //euclidean
                result.push_back(1.0);
                
            }
        }
    }

    return(result);
}

double StarShapedMultipleDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    std::vector<double> interpolation_coefficient_values;
    std::vector<double> distance_values;
    double result = 1.0;
    double x = x2 - x1;
    double y = y2 - y1;
    double z = z2 - z1;

    interpolation_coefficient_values = InterpolationParameter(x1, y1, z1, total_number_of_distances);
    distance_values = MultipleDistanceFunction(x, y, z, total_number_of_distances);
    result = DistanceLinearInterpolation(distance_values, interpolation_coefficient_values);

    return(result);
}

double PeriodicDistanceFunction(double theta, double phi, double d_theta)
{
    double result = 0.0;

    if (sin(phi) * sin(phi) * sin(theta + d_theta) * sin(theta + d_theta) + cos(theta + d_theta) * cos(theta + d_theta) < 0.04)
    {
        result = 1.0;
    }
    else
    {
        result = 0.8;
    }

    return result;
}

double PeriodicDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    int k = 4, n = 4;

    double x_dist = x2 - x1;
    double y_dist = y2 - y1;
    double z_dist = z2 - z1;

    double d = 1.0;

    double phi = pi + atan2(y_dist, x_dist);
    double theta = 0.0, d_theta = 0.0;

    if (z_dist != 0.0)
    {
        if (z_dist > 0)
        {
            theta = atan(sqrt(x_dist * x_dist + y_dist * y_dist) / z_dist);
        }
        else
        {
            theta = pi + atan(sqrt(x_dist * x_dist + y_dist * y_dist) / z_dist);
        }
    }
    else
    {
        theta = 0.5 * pi;
    }

    if (theta <= 0.5 * pi)
    {
        d_theta = acos(1.0 / sqrt(3.0)) - 0.25 * pi;
    }
    else
    {
        d_theta = -(acos(1.0 / sqrt(3.0)) - 0.25 * pi);
    }

    int i = (int)(0.5 * k * phi / pi);
    int j = (int)(0.5 * n * theta / pi);

    double phi_rotation = (2 * i + 1) * pi / k, theta_rotation = 0.0;

    if (n > 2)
    {
        theta_rotation = 0.5 * pi - (2 * j + 1) * pi / n;
    }
    else
    {
        theta_rotation = 0;
    }

    double x_dist_temp = cos(-phi_rotation) * x_dist - sin(-phi_rotation) * y_dist;
    double y_dist_temp = sin(-phi_rotation) * x_dist + cos(-phi_rotation) * y_dist;
    double z_dist_temp = z_dist;

    x_dist = cos(-theta_rotation) * x_dist_temp + sin(-theta_rotation) * z_dist_temp;
    y_dist = y_dist_temp;
    z_dist = -sin(-theta_rotation) * x_dist_temp + cos(-theta_rotation) * z_dist_temp;

    if (z_dist != 0.0)
    {
        if (z_dist > 0)
        {
            theta = atan(sqrt(x_dist * x_dist + y_dist * y_dist) / z_dist);
        }
        else
        {
            theta = pi + atan(sqrt(x_dist * x_dist + y_dist * y_dist) / z_dist);
        }
    }
    else
    {
        theta = 0.5 * pi;
    }

    if ((fabs(phi) > 0.25 * pi) || (theta > 0.75 * pi) || (theta < 0.25 * pi))
    {
        std::cout << "ALARM!" << std::endl;
    }

    d = PeriodicDistanceFunction(theta, phi, d_theta);

    return (sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist) / d);
}

double Distance(double x1, double y1, double z1, double x2, double y2, double z2, double phi, double theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    std::array<double, 3> radius_vector_cartesian;
    double d = 1.0;

    radius_vector_cartesian[0] = x2 - x1;
    radius_vector_cartesian[1] = y2 - y1;
    radius_vector_cartesian[2] = z2 - z1;
    radius_vector_cartesian = PointRotationSpherical(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2], phi, theta);

    if (distance_type == "multiple")
    {
        d = StarShapedMultipleDistance(x1, y1, z1, x2, y2, z2); // Which point to take???
    }

    if (distance_type == "analytical")
    {
        std::array<double, 3> radius_vector_spherical = CartesianToSpherical(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2]);
        d = StarShapedAnalyticalDistance(radius_vector_spherical[1], radius_vector_spherical[2]);
    }

    if (distance_type == "periodic")
    {
        d = StarShapedPeriodicDistance(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2], regular_distance_factor, semiregular_distance_factor, periodic_distance_type, periodic_distance_symmetry);
    }

    return (sqrt(radius_vector_cartesian[0] * radius_vector_cartesian[0] + radius_vector_cartesian[1] * radius_vector_cartesian[1] + radius_vector_cartesian[2] * radius_vector_cartesian[2]) / d);
}

double DistanceForVisualization(double x1, double y1, double z1, double x2, double y2, double z2, double phi, double theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    std::array<double, 3> radius_vector_cartesian;
    double d = 1.0;

    radius_vector_cartesian[0] = x2 - x1;
    radius_vector_cartesian[1] = y2 - y1;
    radius_vector_cartesian[2] = z2 - z1;
    radius_vector_cartesian = PointRotationSpherical(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2], phi, theta);

    if (distance_type == "multiple")
    {
        d = StarShapedMultipleDistance(x1, y1, z1, x2, y2, z2); // Which point to take???
    }

    if (distance_type == "analytical")
    {
        std::array<double, 3> radius_vector_spherical = CartesianToSpherical(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2]);
        d = StarShapedAnalyticalDistance(radius_vector_spherical[1], radius_vector_spherical[2]);
    }

    if (distance_type == "periodic")
    {
        d = StarShapedPeriodicDistance(radius_vector_cartesian[0], radius_vector_cartesian[1], radius_vector_cartesian[2], regular_distance_factor, semiregular_distance_factor, periodic_distance_type, periodic_distance_symmetry);
    }

    return (d);
}

void DistancePrintToFile(double phi, double theta, double max_length, int max_dimension, unsigned int num_phi, unsigned int num_theta, double regular_distance_factor, double semiregular_distance_factor, const std::string& grid_type, const std::string& distance_type, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry, const std::string& save_folder, const std::string& name)
{
    std::vector<std::array<double, 3>> grid_sample, vertices;
    std::vector<std::array<unsigned int, 3>> triangles;
    double x, y, z, radius, rad_max_val = 0.0, dist = 0.0, dist_min_val = 1000000000000000000000.0;
    std::vector<double> x_set, y_set, z_set;
    std::ofstream myfile_distance;
    myfile_distance.open(save_folder + "/distance/" + name + ".csv");
    myfile_distance << "X, Y, Z, R\n";

    for (unsigned int i = 0; i < (2 * (num_phi - 1) * (num_theta - 1)); i++)
    {
        triangles.push_back({ 0, 0, 0 });
    }

    for (unsigned int i = 0; i < num_phi; i++)
    {
        for (unsigned int j = 0; j < num_theta; j++)
        {
            x = sin(pi * j / (num_theta - 1)) * cos(2.0 * pi * i / (num_phi - 1));
            y = sin(pi * j / (num_theta - 1)) * sin(2.0 * pi * i / (num_phi - 1));
            z = cos(pi * j / (num_theta - 1));
            radius = DistanceForVisualization(0.0, 0.0, 0.0, x, y, z, phi, theta, regular_distance_factor, semiregular_distance_factor, distance_type, periodic_distance_type, periodic_distance_symmetry);
            if (radius > rad_max_val) { rad_max_val = radius; }
            x = radius * x;
            y = radius * y;
            z = radius * z;
            x_set.push_back(x);
            y_set.push_back(y);
            z_set.push_back(z);
            vertices.push_back({ x, y, z });
            myfile_distance << x << ", " << y << ", " << z << ", " << radius << "\n";
        }
    }
    myfile_distance.close();

    for (int i = 0; i < num_phi - 1; i++)
    {
        for (int j = 0; j < num_theta - 1; j++)
        {
            triangles[i * ((int)num_theta - 1) + j] = { (unsigned int)((i + 1) * num_theta + j + 1), (unsigned int)((i + 1) * num_theta + j), (unsigned int)(i * num_theta + j) };
            triangles[i * ((int)num_theta - 1) + j + ((int)num_theta - 1) * ((int)num_phi - 1)] = { (unsigned int)(i * num_theta + j + 1), (unsigned int)((i + 1) * num_theta + j + 1), (unsigned int)(i * num_theta + j) };
        }
    }

    for (int j = 0; j < (num_theta - 1); j++)
    {
        triangles.push_back({ (unsigned int)(j + 1), (unsigned int)j, (unsigned int)((num_phi - 1) * num_theta + j) });
        triangles.push_back({ (unsigned int)(j + 1), (unsigned int)((num_phi - 1) * num_theta + j + 1), (unsigned int)((num_phi - 1) * num_theta + j) });
    }

    SaveToSTL(vertices, triangles, save_folder + std::string("/distance"), name);

    grid_sample = PrintGridToFile(max_length, max_dimension, grid_type);

    for (const auto& point1 : grid_sample)
    {
        for (const auto& point2 : grid_sample)
        {
            dist = EucledianDistance(point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]);
            if ((dist < dist_min_val) && (dist > 0.0)) { dist_min_val = dist; }
        }
    }

    myfile_distance.open(save_folder + "/distance/with_seeds_grid_" + name + ".csv");
    myfile_distance << "X, Y, Z\n";
    for (const auto& grid_point : grid_sample)
    {
        for (unsigned int i = 0; i < x_set.size(); i++)
        {
            myfile_distance << grid_point[0] + 0.25 * dist_min_val * x_set[i] / rad_max_val << ", " << grid_point[1] + 0.25 * dist_min_val * y_set[i] / rad_max_val << ", " << grid_point[2] + 0.25 * dist_min_val * z_set[i] / rad_max_val << "\n";
        }
    }
    myfile_distance.close();

    std::vector<std::array<double, 3>> computation_grid_sample;
    std::array<double, 3> sample_length;
    std::array<int, 3> sample_size;
    CellGridInitialization(computation_grid_sample, sample_length, sample_size, max_length, max_dimension, grid_type, distance_type);
    if (distance_type != "multiple")
    {
        CellGridExpansion(computation_grid_sample, sample_length, sample_size);
    }
}

void PrintMultiDistInterpolationToFile(unsigned int num_phi, unsigned int num_theta, unsigned int discretization_number, const std::string& save_folder, const std::string& name)
{
    std::ofstream myfile_distance;
    std::ostringstream id;
    std::string id_string;
    double x, y, z, radius, coefficient = 0.0;

    for (unsigned int i = 0; i < discretization_number; i++)
    {
        id << i;
        id_string = id.str();

        coefficient = (double)i / (double)(discretization_number - 1);

        myfile_distance.open(save_folder + "/multi_distance/" + name + ".csv" + id_string);
        myfile_distance << "X, Y, Z\n";
        for (unsigned int i = 0; i < num_phi; i++)
        {
            for (unsigned int j = 0; j < num_theta; j++)
            {
                x = sin(pi * j / (num_theta - 1)) * cos(2.0 * pi * i / (num_phi - 1));
                y = sin(pi * j / (num_theta - 1)) * sin(2.0 * pi * i / (num_phi - 1));
                z = cos(pi * j / (num_theta - 1));
                radius = StarShapedMultipleDistance(0.0, 0.0, 0.0, x, y, z);
                x = radius * x;
                y = radius * y;
                z = radius * z;
                myfile_distance << x << ", " << y << ", " << z << ", " << radius << "\n";
            }
        }
        myfile_distance.close();
        id.str("");
        id.clear();
    }
}

void PrintPolyhedron(unsigned int discretization_number, const std::string& periodic_distance_type, const std::string& periodic_distance_symmetry)
{
    std::vector<std::array<double, 3>> points;
    std::vector<std::array<double, 3>> vertices, face_centers;
    std::array<double, 3> temporary_point, face_closest;
    std::vector<double> coverage_radius_distribution;
    double neighbour_radius = 0, coverage_radius = 0.0, epsilon = 0.01, x, y, z, minimum_radius = 100.0, identificator;

    if (periodic_distance_type == "regular")
    {
        if (periodic_distance_symmetry == "tetrahedral")
        {
            for (const auto& point : regular_tetrahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_dual_tetrahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = 2.0 + epsilon;
            coverage_radius = regular_tetrahedron_face_radius;
        }
        if (periodic_distance_symmetry == "cubic")
        {
            for (const auto& point : regular_hexahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_octahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = sqrt(2.0) + epsilon;
            coverage_radius = regular_hexahedron_face_radius;
        }
        if (periodic_distance_symmetry == "octahedral")
        {
            for (const auto& point : regular_octahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_hexahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = sqrt(2.0) + epsilon;
            coverage_radius = regular_octahedron_face_radius;
        }
        if (periodic_distance_symmetry == "dodecahedral")
        {
            for (const auto& point : regular_dodecahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_icosahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = (sqrt(5.0) - 1.0) / sqrt(3.0) + epsilon;
            coverage_radius = regular_dodecahedron_face_radius;
        }
        if (periodic_distance_symmetry == "icosahedral")
        {
            for (const auto& point : regular_icosahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_dodecahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = 2.0 / sqrt(2.5 + 0.5 * sqrt(5.0)) + epsilon;
            coverage_radius = regular_icosahedron_face_radius;
        }
        if (periodic_distance_symmetry == "diploidal")
        {
            for (const auto& point : regular_hexahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : regular_octahedron_vertices) { face_centers.push_back(point); }
            neighbour_radius = sqrt(2.0) + epsilon;
            coverage_radius = regular_hexahedron_face_radius;
        }
    }
    if (periodic_distance_type == "semiregular")
    {
        if (periodic_distance_symmetry == "truncated_tetrahedral")
        {
            for (const auto& point : truncated_tetrahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : truncated_tetrahedron_dual_vertices) { face_centers.push_back(point); }
            for (const auto& element : truncated_tetrahedron_face_types_distance_compensation) { coverage_radius_distribution.push_back(element); }
            neighbour_radius = sqrt(8.0 / 11.0) + epsilon;
        }
        if (periodic_distance_symmetry == "cuboctahedral")
        {
            for (const auto& point : cuboctahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : rhombic_dodecahedron_vertices) { face_centers.push_back(point); }
            for (const auto& element : cuboctahedron_face_types_distance_compensation) { coverage_radius_distribution.push_back(element); }
            neighbour_radius = 1.0 + epsilon;
        }
        if (periodic_distance_symmetry == "icosidodecahedral")
        {
            for (const auto& point : icosidodecahedron_vertices) { vertices.push_back(point); }
            for (const auto& point : rhombic_triacontahedron_vertices) { face_centers.push_back(point); }
            for (const auto& element : icosidodecahedron_face_types_distance_compensation) { coverage_radius_distribution.push_back(element); }
            neighbour_radius = sqrt(id_norm * id_norm + 0.25 * (1.0 - 2.0 * golden_ratio)) / id_norm + epsilon;
        }
        if (periodic_distance_symmetry == "snub_cubic")
        {
            for (const auto& point : snub_cube_vertices) { vertices.push_back(point); }
            for (const auto& point : pentagonal_icositetrahedron_vertices) { face_centers.push_back(point); }
            for (const auto& element : snub_cube_face_types_distance_compensation) { coverage_radius_distribution.push_back(element); }
            neighbour_radius = sqrt(2.0 + 4.0 * tribonacci_constant - 2.0 * tribonacci_constant * tribonacci_constant) / sqrt(1.0 + 1.0 / (tribonacci_constant * tribonacci_constant) + tribonacci_constant * tribonacci_constant) + epsilon;
        }
    }

    for (const auto& point : vertices) { points.push_back(point); }

    for (const auto& point1 : vertices)
    {
        for (const auto& point2 : vertices)
        {
            if (point1 != point2)
            {
                if (EucledianDistance(point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]) <= neighbour_radius)
                {
                    for (unsigned int i = 1; i < (discretization_number - 1); i++)
                    {
                        temporary_point[0] = point1[0] + (point2[0] - point1[0]) * (1.0 / (double)(discretization_number - 1)) * (double)i;
                        temporary_point[1] = point1[1] + (point2[1] - point1[1]) * (1.0 / (double)(discretization_number - 1)) * (double)i;
                        temporary_point[2] = point1[2] + (point2[2] - point1[2]) * (1.0 / (double)(discretization_number - 1)) * (double)i;
                        points.push_back(temporary_point);
                    }
                }
            }
        }
    }

    std::ofstream myfile_distance;
    myfile_distance.open("polyhedron_ordinary.csv");
    myfile_distance << "X, Y, Z\n";

    for (const auto& point : points)
    {
        myfile_distance << point[0] << ", " << point[1] << ", " << point[2] << "\n";
    }

    myfile_distance.close();

    myfile_distance.open("polyhedron_spherical.csv");
    myfile_distance << "X, Y, Z\n";

    for (const auto& point : points)
    {
        myfile_distance << point[0] / sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]) << ", " << point[1] / sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]) << ", " << point[2] / sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]) << "\n";
    }

    myfile_distance.close();

    myfile_distance.open("polyhedron_coverage.csv.");
    myfile_distance << "X, Y, Z, ID\n";
    unsigned int coverage_discretization = 2 * discretization_number;
    if (periodic_distance_type == "regular")
    {
        for (unsigned int i = 0; i < coverage_discretization; i++)
        {
            for (unsigned int j = 0; j < coverage_discretization; j++)
            {
                minimum_radius = 100.0;
                x = sin(pi * j / (coverage_discretization - 1)) * cos(2.0 * pi * i / (coverage_discretization - 1));
                y = sin(pi * j / (coverage_discretization - 1)) * sin(2.0 * pi * i / (coverage_discretization - 1));
                z = cos(pi * j / (coverage_discretization - 1));
                for (const auto& point : face_centers)
                {
                    if (EucledianDistance(point[0], point[1], point[2], x, y, z) < minimum_radius)
                    {
                        face_closest = point;
                        minimum_radius = EucledianDistance(point[0], point[1], point[2], x, y, z);
                    }
                }
                for (unsigned int l = 0; l < face_centers.size(); l++)
                {
                    if (face_centers[l] == face_closest) { identificator = (double)l; }
                }
                myfile_distance << x << ", " << y << ", " << z << ", " << identificator << "\n";
            }
        }
    }
    if (periodic_distance_type == "semiregular")
    {
        for (unsigned int i = 0; i < coverage_discretization; i++)
        {
            for (unsigned int j = 0; j < coverage_discretization; j++)
            {
                minimum_radius = 100.0;
                x = sin(pi * j / (coverage_discretization - 1)) * cos(2.0 * pi * i / (coverage_discretization - 1));
                y = sin(pi * j / (coverage_discretization - 1)) * sin(2.0 * pi * i / (coverage_discretization - 1));
                z = cos(pi * j / (coverage_discretization - 1));
                for (unsigned int l = 0; l < face_centers.size(); l++)
                {
                    temporary_point = face_centers[l];
                    if (SquaredEucledianDistance(temporary_point[0], temporary_point[1], temporary_point[2], x, y, z) + coverage_radius_distribution[l] < minimum_radius)
                    {
                        face_closest = temporary_point;
                        minimum_radius = SquaredEucledianDistance(temporary_point[0], temporary_point[1], temporary_point[2], x, y, z) + coverage_radius_distribution[l];
                        identificator = (double)l;
                    }
                }
                myfile_distance << x << ", " << y << ", " << z << ", " << identificator << "\n";
            }
        }
    }
    
    myfile_distance.close();
}