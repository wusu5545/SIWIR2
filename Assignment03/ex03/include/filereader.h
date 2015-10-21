#ifndef FILEREADER_H
#define FILEREADER_H

#include<string>
#include<map>
#include <typeinfo>

using std::string;
using std::map;
using std::stod;
using std::stoul;

class FileReader
{
  public:
    void readParameters(string filename);
    template<typename Type>
    Type getParameter(const string str);
    template<typename Type>
    Type getFileName(const string str);
  private:
    map<string,string> parameter_;
};

template<typename Type>
Type FileReader::getParameter(const string str)
{
  Type res;
  if (typeid(Type) == typeid(size_t))
    res = stoul(parameter_.at(str));
  else if (typeid(Type) == typeid(double))
    res = stod(parameter_.at(str));
  return res;
}
template<typename Type>
Type FileReader::getFileName(const string str)
{
  if (parameter_.find(str) == parameter_.end())
    return "";
  else return parameter_.at(str);
}

#endif//FILEREADER_H
