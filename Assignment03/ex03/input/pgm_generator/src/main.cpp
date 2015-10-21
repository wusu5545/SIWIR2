#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc <=1){
    cout<<"Please input sizex and sizey"<<endl;
    exit(0);
  }
  ofstream output;
  output.open("../random.pgm");
  
  size_t sizex,sizey;
  
  sizex = atoi(argv[1]);
  sizey = atoi(argv[2]);
  
  output<<"P2"<<endl;
  output<<"# random.pgm"<<endl;
  output<<sizex<<" "<<sizey<<endl;
  output<<255<<endl;
  
  for (size_t i=0;i<sizey;++i){
    for (size_t j=0;j<sizex;++j){
      if (double(i)/sizey>=0.25&&double(i)/sizey<=0.75 && double(j)/sizex>=0.25 && double(j)/sizex<=0.75)
	output<<setw(4)<<(rand()%2 == 0?0:255);
      else
	output<<setw(4)<<255;
    }
    output<<endl;
  }
  
  output.close();
  
  return 0;
}
