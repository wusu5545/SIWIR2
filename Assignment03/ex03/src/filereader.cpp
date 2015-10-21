#include "filereader.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using std::fstream;

void FileReader::readParameters(string filename)
{
  fstream input;
  input.open(filename,fstream::in);
  string str1,str2;
  while (!input.eof()){
    input>>str1>>str2;
    (*this).parameter_[str1] = str2;
  }
  input.close();
}
