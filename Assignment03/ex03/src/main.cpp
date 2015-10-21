#include "functions.h"

int main(int argc, char* argv[])
{
  if (argc <=1)
  {
    cout<<"Please input Parameters File Name in input folder"<<endl;
    exit(0);
  }
  // Parsing the parameter file
  FileReader reader;
  reader.readParameters(argv[1]);

  Solver(reader);
  return 0;
}
