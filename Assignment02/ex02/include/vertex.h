#ifndef VERTEX_H
#define VERTEX_H

#include <map>
#include <vector>
#include <math.h>

using namespace std;
struct matrix
{
  double stiffness;
  double mass;
};

struct face
{
  size_t v0,v1,v2;
};

struct matrixvector
{
  vector<double> stiff;
  vector<double> mass;
};

class Vertex
{
  public:
    explicit Vertex();
    ~ Vertex();
    inline double x() const;
    double & x();
    inline double y() const;
    double & y();
    inline double k() const;
    double & k();
    inline matrix operator[](const size_t j) const;
    matrix &operator[](const size_t j);
    bool operator()(const size_t j);
    map<size_t,matrix>::const_iterator b() const;
    map<size_t,matrix>::const_iterator e() const;
  private:
    double x_,y_,k_;
    map<size_t, matrix> neighbour_;
};

#endif //VERTEX_H
