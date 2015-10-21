#include "vertex.h"

Vertex::Vertex()
{
  
}

Vertex::~Vertex()
{
}

inline double Vertex::x() const
{
  return x_;
}

double& Vertex::x()
{
  return x_;
}

inline double Vertex::y() const
{
  return y_;
}

double& Vertex::y()
{
  return y_;
}

inline double Vertex::k() const
{
  return k_;
}

double & Vertex::k()
{
  return k_;
}

inline matrix Vertex::operator[](const size_t j) const
{
  return neighbour_.at(j);
}

matrix & Vertex::operator[](const size_t j)
{
  return neighbour_[j];
}

bool Vertex::operator()(const size_t j)
{
  if (neighbour_.find(j)==neighbour_.end())
    return true;
  else 
    return false;
}

map<size_t,matrix>::const_iterator Vertex::b() const
{
  return neighbour_.begin();
}

map<size_t,matrix>::const_iterator Vertex::e() const
{
  return neighbour_.end();
}