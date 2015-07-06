#ifndef ParameterReader
#define ParameterReader

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>



class Parameter_Reader 

{ 
public:
	std::map<std::string, std::string> parameters; 
	
	//CONSTRUCTOR
	Parameter_Reader();

	//DESTRUCTOR
	~Parameter_Reader();

	void readParameters(const std::string& filename);
	
	inline bool IsDefined(const std::string& key) const;
	
	template<typename Type> inline void GetParameter(const std::string& key, Type &value) const {
	  std::string temp1 = parameters.at(key);
	  std::stringstream temp (temp1);
	  temp >> value;
	}
};


#endif