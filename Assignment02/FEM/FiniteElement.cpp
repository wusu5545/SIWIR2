#include "FiniteElement.h"
#include <vector>
#include <cmath>

#include <iostream>

using std::cout;

using std::vector;

extern const double Euler=2.718281828459045;
extern double delta;
extern double sigma;


double FiniteElement::kSqrValue(double x, double y)
{
    double kSqr=0.0;
    
    kSqr=(100+delta)*pow(Euler, -50.0*(x*x+y*y))-100.0;

    return kSqr;
}

vector<double> FiniteElement::multiplyA(const vector<vertex> &A, const vector<double> &u)
{
    vector<double> result(u.size(), 0.0);
    for (int i=0; i<=A.size()-1; ++i)
    {
	for (map<int,stclMass>::const_iterator itr=A[i].neighbour.begin(); itr!=A[i].neighbour.end(); ++itr)
	    result[i]+=u[itr->first]*((itr->second).stencil);
    }
    return result;
}

vector<double> FiniteElement::multiplyM(const vector<vertex> &M, const vector<double> &u)
{
    vector<double> result(u.size(), 0.0);
    for (int i=0; i<=M.size()-1; ++i)
    {
	for (map<int,stclMass>::const_iterator itr=M[i].neighbour.begin(); itr!=M[i].neighbour.end(); ++itr)
	    result[i]+=u[itr->first]*((itr->second).mass);
    }
    return result;
}

vector<double> FiniteElement::operator+(const vector<double> &l, const vector<double> &r)
{
    vector<double> result=l;
    for (int i=0; i<=r.size()-1; ++i)
        result[i]+=r[i];
    return result;
}

vector<double> FiniteElement::operator-(const vector<double> &l, const vector<double> &r)
{
    vector<double> result=l;
    for (int i=0; i<=r.size()-1; ++i)
        result[i]-=r[i];
    return result;
}

vector<double> FiniteElement::operator*(const double &factor, const vector<double> &r)
{
    vector<double> result=r;
    for (int i=0; i<=r.size()-1; ++i)
	    result[i]=factor*r[i];
    return result;
}
double FiniteElement::operator*(const vector<double> &l, const vector<double> &r)
{
    double result=0.0;
    for (int i=0; i<=l.size()-1; ++i)
	result+=(l[i]*r[i]);
	return result;
}

double FiniteElement::norm(const vector<double> &v)
{
	double result=0.0;
	for (int i=0; i<=v.size()-1; ++i)
	result+=(v[i]*v[i]);
	result=sqrt(result);
	return result;
}


vector<double> FiniteElement::CG(const vector<vertex> &A, const vector<double> &b)
{
    vector<double> x(b.size(), 0.0);
    vector<double> r_old;
    vector<double> r_new=b;
    vector<double> p=r_new;
    vector<double> tempMultiply;
    double alpha=0.0;
    double beta=0.0;
    double normX=0.0;

    do
    {
	tempMultiply=multiplyA(A, p);
        r_old=r_new;
        alpha=(r_old*r_old)/(p*tempMultiply);
        x=x+alpha*p;
        r_new=r_old-alpha*tempMultiply;
        beta=(r_new*r_new)/(r_old*r_old);
        p=r_new+beta*p;
        normX=norm(b-multiplyA(A, x));
    }
    while (normX>=1.0e-14);

    return x;
}

vector<double> FiniteElement::IPI(const vector<vertex> &A, double &egv)
{
    vector<double> eigenVector(A.size(), 0.0);
    vector<double> tempMultiplyA;
    vector<double> tempMultiplyM;
    double eigenValue=0.0;
    eigenVector[0]=1.0;

    do
    {
	tempMultiplyM=multiplyM(A, eigenVector);
        eigenVector=CG(A, tempMultiplyM);
        eigenVector=(1/norm(eigenVector))*eigenVector;
	tempMultiplyA=multiplyA(A, eigenVector);
	tempMultiplyM=multiplyM(A, eigenVector);
        eigenValue=(eigenVector*tempMultiplyA)/(eigenVector*tempMultiplyM);
    }
    while (norm(tempMultiplyA-eigenValue*tempMultiplyM)>=sigma);
    
    egv=eigenValue;

    return eigenVector;
}