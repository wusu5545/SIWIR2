#pragma once
#include <map>
#include <vector>

using std::map;
using std::vector;

namespace FiniteElement
{
    struct stclMass
    {
        double stencil;
        double mass;
    };

    struct vertex
    {
        double x_c;
        double y_c;
        double kSqr;
        map<int, stclMass> neighbour;
    };

    struct face
    {
        int v1;
        int v2;
        int v3;
    };

    double kSqrValue(double, double);
    vector<double> multiplyA(const vector<vertex>&, const vector<double>&);
    vector<double> multiplyM(const vector<vertex>&, const vector<double>&);
    vector<double> operator+(const vector<double>&, const vector<double>&);
    vector<double> operator-(const vector<double>&, const vector<double>&);
    vector<double> operator*(const double&, const vector<double>&);
    double operator*(const vector<double>&, const vector<double>&);
    double norm(const vector<double>&);
    vector<double> CG(const vector<vertex>&, const vector<double>&);
    vector<double> IPI(const vector<vertex>&, double&);
}
