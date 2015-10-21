#include"functions.h"

double delta;
double sigma;
size_t lvlref;

int main(int argc, char* argv[])
{
  if (argc == 1)
  {
    cout<< "Args: delta and eps" <<endl;
    exit(0);
  }
  else{
    delta = strtod(argv[1],NULL);
    sigma = strtod(argv[2],NULL);
    if (argv[3]==NULL)
      lvlref = 0;
    else 
      lvlref = stoi(argv[3]);
    Waveguide("input/unit_circle.txt");
  }
  return 0;
}
