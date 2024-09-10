#ifndef MATHUTILS_HPP_
#define MATHUTILS_HPP_

#include<string.h>
#include<sstream>
#include<vector>
//#include"mathutils.cpp"

using std::string;

//just a header file for miscellaneous functions to avoid cluttering class-specific files
double cubic_interp(double r, double f0, double f1, double f2, double f3, int j0, double dr);
int bound(int n, int lower, int upper);
double fivePointDeriv(double step, int order, double f1, double f2, double f3, double f4, double f5);
double sevenPointDeriv(double step, int order, double f1, double f2, double f3, double f4, double f5, double f6, double f7);

double lagrange_interp(double x, const std::vector<double>& X, const std::vector<double>& Y);

//Below are templated, so we put these in the header directly

//stencil for a derivative evaluated about the 4th position on uniform grid
template <typename T>
T p4_stencil(double step, T f1, T f2, T f3, T f4, T f5)
{
    return (-(1.0 / 12.0) * f1 + (1.0 / 2.0) * f2 - (3.0 / 2.0) * f3 + (5.0 / 6.0) * f4 + (1.0 / 4.0) * f5) / step;
}
template <typename T>
T p5_stencil(double step, T f1, T f2, T f3, T f4, T f5)
{
    return ((1.0 / 4.0) * f1 - (4.0 / 3.0) * f2 + (3.0) * f3 - (4.0) * f4 + (25.0 / 12.0) * f5) / step;
}

template <typename T>
T sevenPointDeriv(double step, int order, T f1, T f2, T f3, T f4, T f5, T f6, T f7)
{
    switch (order)
    {
        case 0:
            return f4;
        case 1:
            return ((-1.0/60.0)*f1 + (3.0/20.0)*f2 - (3.0/4.0)*f3 + (3.0/4.0)*f5 -(3.0/20.0)*f6 + (1.0/60.0)*f7 )/step;
        case 2:
            return ((1.0/90.0)*f1 -(3.0/20.0)*f2 + (3.0/2.0)*f3 - (49.0/18.0)*f4 + (3.0/2.0)*f5 -(3.0/20.0)*f6 + (1.0/90.0)*f7)/pow(step,2);
        case 3:
            return ((1.0/8.0)*f1 -(1)*f2 + (13.0/8.0)*f3 - (13.0/8.0)*f5 +(1)*f6 - (1.0/8.0)*f7)/pow(step,3);
        case 4:
            return (-(1.0/6.0)*f1 +(2.0)*f2 - (13.0/2.0)*f3 + (28.0/3.0)*f4 - (13.0/2.0)*f5 +(2.0)*f6 - (1.0/6.0)*f7)/pow(step,4);
        case 5:
            return ((-1.0/2.0)*f1 +(2.0)*f2 - (5.0/2.0)*f3  + (5.0/2.0)*f5 -(2.0)*f6 + (1.0/2.0)*f7)/pow(step,5);
        case 6:
            return ((1.0)*f1 -(6.0)*f2 + (15.0)*f3 - (20.0)*f4 + (15.0)*f5 -(6.0)*f6 + (1.0)*f7)/pow(step,6);
        default:
            printf("ERROR: invalid derivative order requested in sevenPointDeriv");
            abort();

    }
}


template <typename T>
void fill_parameter (string& current_line, string line_start, T& parameter, bool quiet)
{
    if (current_line.find(line_start) != string::npos)
    {
        // Create a substring starting from the position right after line_start
        size_t pos = current_line.find(line_start);
        string rest_of_line = current_line.substr(pos + line_start.length());

        // Create a stringstream from the rest_of_line
        std::stringstream ss(rest_of_line);

        // Extract the double value from the stringstream
        if (ss >> parameter) {
            if (!quiet) std::cout << "Read in " << line_start << parameter << std::endl;
        } else if (!(ss >> parameter)) {
            std::cout << "WARNING: Failed to extract value for parameter " << line_start << std::endl;
        }
    }
}
#endif /* MATHUTILS_HPP_ */
