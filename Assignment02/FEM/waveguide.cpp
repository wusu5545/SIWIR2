#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "FiniteElement.h"
#include "Colsamm.h"

using namespace ::_COLSAMM_;

using std::string;
using std::strtod;
using std::fstream;
using std::ios;
using std::vector;
using std::map;
using std::pair;
using FiniteElement::stclMass;
using FiniteElement::vertex;
using FiniteElement::face;
using FiniteElement::kSqrValue;
using FiniteElement::IPI;
using ELEMENTS::Triangle;
using std::endl;

const double Euler=2.718281828459045;
double delta;
double sigma;

int main(int argc, char *argv[], char *envp[])
{
    const string inputFile="inputs/unit_circle.txt";
    const string outputFile1="ksq.txt";
    const string outputFile2="A.txt";
    const string outputFile3="M.txt";
    const string outputFile4="eigenmode.txt";
    const string outputFile5="lambda.txt";
    
    delta=strtod(argv[1],NULL);
    sigma=strtod(argv[2],NULL);
    
    fstream in;
    vector<vertex> vertexVector;
    vector<face>   faceVector;
    in.open(inputFile.c_str(),fstream::in);
    
    if (in.is_open())
    {
        in.ignore(76);
        //Read in vertex data.
        for (int i=1; i<=1039; ++i)
        {
            vertex v;
            in>>v.x_c;
            in>>v.x_c;
            in>>v.y_c;
            v.kSqr=kSqrValue(v.x_c, v.y_c);
            vertexVector.push_back(v);
        }

        in.ignore(86);
        //Read in face data and build neighbour data for each vertex.
        for (int i=1; i<=1976; ++i)
        {
            //Read in face data.
            face f;
            in>>f.v1;
            in>>f.v2;
            in>>f.v3;

            //Compute stencil matrix and mass matrix.
            Triangle faceElement;
            vector<vector<double> > localStcl;
            vector<vector<double> > localKsqr;
            vector<vector<double> > localMass;
            vector<double> position(6, 0.0);

            position[0]=vertexVector[f.v1].x_c;
            position[1]=vertexVector[f.v1].y_c;
            position[2]=vertexVector[f.v2].x_c;
            position[3]=vertexVector[f.v2].y_c;
            position[4]=vertexVector[f.v3].x_c;
            position[5]=vertexVector[f.v3].y_c;

            faceElement(position);
            localStcl=faceElement.integrate(grad(v_())*grad(w_()));
            localKsqr=faceElement.integrate(func<double>(kSqrValue)*v_()*w_());
            localMass=faceElement.integrate(v_()*w_());
            
            //Build neighbour data.
            //00eigenVector
            if (vertexVector.at(f.v1).neighbour.find(f.v1)==
                    vertexVector.at(f.v1).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[0][0]-localKsqr[0][0];
                n.mass=localMass[0][0];
                vertexVector.at(f.v1).neighbour.insert(pair<int, stclMass>(f.v1, n));
            }
            else
            {
                vertexVector.at(f.v1).neighbour.at(f.v1).stencil+=(localStcl[0][0]-localKsqr[0][0]);
                vertexVector.at(f.v1).neighbour.at(f.v1).mass+=localMass[0][0];
            }
            //01
            if (vertexVector.at(f.v1).neighbour.find(f.v2)==
                    vertexVector.at(f.v1).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[0][1]-localKsqr[0][1];
                n.mass=localMass[0][1];
                vertexVector.at(f.v1).neighbour.insert(pair<int, stclMass>(f.v2, n));
            }
            else
            {
                vertexVector.at(f.v1).neighbour.at(f.v2).stencil+=(localStcl[0][1]-localKsqr[0][1]);
                vertexVector.at(f.v1).neighbour.at(f.v2).mass+=localMass[0][1];
            }
            //02
            if (vertexVector.at(f.v1).neighbour.find(f.v3)==
                    vertexVector.at(f.v1).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[0][2]-localKsqr[0][2];
                n.mass=localMass[0][2];
                vertexVector.at(f.v1).neighbour.insert(pair<int, stclMass>(f.v3, n));
            }
            else
            {
                vertexVector.at(f.v1).neighbour.at(f.v3).stencil+=(localStcl[0][2]-localKsqr[0][2]);
                vertexVector.at(f.v1).neighbour.at(f.v3).mass+=localMass[0][2];
            }
            //10
            if (vertexVector.at(f.v2).neighbour.find(f.v1)==
                    vertexVector.at(f.v2).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[1][0]-localKsqr[1][0];
                n.mass=localMass[1][0];
                vertexVector.at(f.v2).neighbour.insert(pair<int, stclMass>(f.v1, n));
            }
            else
            {
                vertexVector.at(f.v2).neighbour.at(f.v1).stencil+=(localStcl[1][0]-localKsqr[1][0]);
                vertexVector.at(f.v2).neighbour.at(f.v1).mass+=localMass[1][0];
            }
            //11
            if (vertexVector.at(f.v2).neighbour.find(f.v2)==
                    vertexVector.at(f.v2).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[1][1]-localKsqr[1][1];
                n.mass=localMass[1][1];
                vertexVector.at(f.v2).neighbour.insert(pair<int, stclMass>(f.v2, n));
            }
            else
            {
                vertexVector.at(f.v2).neighbour.at(f.v2).stencil+=(localStcl[1][1]-localKsqr[1][1]);
                vertexVector.at(f.v2).neighbour.at(f.v2).mass+=localMass[1][1];
            }
            //12
            if (vertexVector.at(f.v2).neighbour.find(f.v3)==
                    vertexVector.at(f.v2).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[1][2]-localKsqr[1][2];
                n.mass=localMass[1][2];
                vertexVector.at(f.v2).neighbour.insert(pair<int, stclMass>(f.v3, n));
            }
            else
            {
                vertexVector.at(f.v2).neighbour.at(f.v3).stencil+=(localStcl[1][2]-localKsqr[1][2]);
                vertexVector.at(f.v2).neighbour.at(f.v3).mass+=localMass[1][2];
            }
            //20
            if (vertexVector.at(f.v3).neighbour.find(f.v1)==
                    vertexVector.at(f.v3).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[2][0]-localKsqr[2][0];
                n.mass=localMass[2][0];
                vertexVector.at(f.v3).neighbour.insert(pair<int, stclMass>(f.v1, n));
            }
            else
            {
                vertexVector.at(f.v3).neighbour.at(f.v1).stencil+=(localStcl[2][0]-localKsqr[2][0]);
                vertexVector.at(f.v3).neighbour.at(f.v1).mass+=localMass[2][0];
            }
            //21
            if (vertexVector.at(f.v3).neighbour.find(f.v2)==
                    vertexVector.at(f.v3).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[2][1]-localKsqr[2][1];
                n.mass=localMass[2][1];
                vertexVector.at(f.v3).neighbour.insert(pair<int, stclMass>(f.v2, n));
            }
            else
            {
                vertexVector.at(f.v3).neighbour.at(f.v2).stencil+=(localStcl[2][1]-localKsqr[2][1]);
                vertexVector.at(f.v3).neighbour.at(f.v2).mass+=localMass[2][1];
            }
            //22
            if (vertexVector.at(f.v3).neighbour.find(f.v3)==
                    vertexVector.at(f.v3).neighbour.end())
            {
                stclMass n;
                n.stencil=localStcl[2][2]-localKsqr[2][2];
                n.mass=localMass[2][2];
                vertexVector.at(f.v3).neighbour.insert(pair<int, stclMass>(f.v3, n));
            }
            else
            {
                vertexVector.at(f.v3).neighbour.at(f.v3).stencil+=(localStcl[2][2]-localKsqr[2][2]);
                vertexVector.at(f.v3).neighbour.at(f.v3).mass+=localMass[2][2];
            }
        }
    }
    
    fstream out;
    //Write out k square value.
    out.open(outputFile1.c_str(),fstream::out);
    if (out.is_open())
    {
	for (int i=0; i<=vertexVector.size()-1; ++i)
	    out<<vertexVector[i].x_c<<' '<<vertexVector[i].y_c<<' '<<vertexVector[i].kSqr<<' '<<endl;
    }
    out.close();
    
    //Write out stencil matrix.
    out.open(outputFile2.c_str(),fstream::out);
    if (out.is_open())
    {
	for (int i=0; i<=vertexVector.size()-1; ++i)
	{
	    for (map<int,stclMass>::const_iterator itr=vertexVector[i].neighbour.begin(); itr!=vertexVector[i].neighbour.end(); ++itr)
	        out<<i<<' '<<(itr->first)<<' '<<(itr->second).stencil<<endl;
	}
    }
    out.close();
    
    //Write out mass matrix.
    out.open(outputFile3.c_str(),fstream::out);
    if (out.is_open())
    {
	for (int i=0; i<=vertexVector.size()-1; ++i)
	{
	    for (map<int,stclMass>::const_iterator itr=vertexVector[i].neighbour.begin(); itr!=vertexVector[i].neighbour.end(); ++itr)
	        out<<i<<' '<<(itr->first)<<' '<<(itr->second).mass<<endl;
	}
    }
    out.close();
    
    //Calculate the smallest eigen value.
    vector<double> eigenVector;
    double eigenValue=0.0;
    eigenVector=IPI(vertexVector, eigenValue);
    
    //Write out eigen vector.
    out.open(outputFile4.c_str(),fstream::out);
    if (out.is_open())
    {
	for (int i=0; i<=vertexVector.size()-1; ++i)
	    out<<vertexVector[i].x_c<<' '<<vertexVector[i].y_c<<' '<<eigenVector[i]<<endl;
    }
    out.close();
    
    //Write out eigen value.
    out.open(outputFile5.c_str(),fstream::out);
    out<<eigenValue<<endl;
    out.close();
    
    return 0;
}
