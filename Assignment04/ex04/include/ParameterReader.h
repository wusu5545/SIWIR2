#pragma once

#include <string>
#include <map>
#include <sstream>

using namespace std;

class ParameterReader
{
  public:
    //constructor
    ParameterReader();
    //deconstructor
    ~ParameterReader();
    
    void readParameters(const string & filename);
    
    inline bool IsDefined(const string & key) const;
    
    template<typename Type>
    inline void GetParameter(const string & key,Type &value) const
    {
      string temp = parameters.at(key);
      stringstream ParOut (temp);
      ParOut >> value;
    }
  private:
    map<string,string> parameters;
};