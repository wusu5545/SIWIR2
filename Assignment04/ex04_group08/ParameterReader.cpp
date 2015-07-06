#include "ParameterReader.hpp"

#include <sstream>
#include <ios>
#include <iterator>
#include <cmath>
#include <string>
#include <map>
#include <typeinfo>


Parameter_Reader::~Parameter_Reader()
{
  
}

Parameter_Reader::Parameter_Reader()
{
}

void Parameter_Reader::readParameters(const std::string& filename)
{
	std::string key;
	std::string value;

	std::ifstream par_file(filename.c_str());
	
	while(!par_file.eof())
	{
		par_file >> key;

		par_file >> value;
		
		parameters.insert(std::pair<std::string, std::string>(key,value));
	}
}

bool Parameter_Reader::IsDefined(const std::string& key) const
{
	std::map<std::string, std::string>::iterator it;
	if(parameters.find(key)!=parameters.end())
	{return true;}
	else {
	  return false;
	}
}

