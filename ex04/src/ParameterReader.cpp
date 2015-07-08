#include "ParameterReader.h"

#include <iostream>
#include <fstream>

ParameterReader::ParameterReader(){}
ParameterReader::~ParameterReader(){}

void ParameterReader::readParameters(const string & filename)
{
  string key;
  string value;
  
  ifstream input(filename);
  
  while (!input.eof())
  {
    input >> key;
    input >> value;
    parameters[key] = value;
  }
}

inline bool ParameterReader::IsDefined(const string & key) const
{
  if (parameters.find(key)!=parameters.end())
    return true;
  else
    return false;
}