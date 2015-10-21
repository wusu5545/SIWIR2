#include"functions.h"

#include"Source/Colsamm.h"

using namespace ::_COLSAMM_;

extern size_t lvlref;

void Waveguide( const string filename)
{

  vector<Vertex> vertex;
  fstream input;
  input.open(filename,fstream::in);
  
  size_t N;
  input>>N;
  
  input.ignore(70);
  //vertex = boost::make_shared<Vertex[]>;
  
  size_t index;
  Vertex v;
  for (size_t i=0;i<N;++i){
    input>>index>>v.x()>>v.y();
    v.k() = ksq(v.x(),v.y());
    vertex.push_back(v);
    //cout << index<<' '<<vertex[index].x()<<' '<<vertex[index].y()<<endl;
  }


  //deal with faces
  size_t face_N;
  size_t v0,v1,v2;
  input >> face_N;
  vector<face> faces(face_N);
  input.ignore(80);
  for (size_t i=0;i<face_N;++i){
    input>>faces[i].v0>>faces[i].v1>>faces[i].v2;
  }

  //do refinement
  index = N;
  size_t v01,v02,v12;
  for (size_t lvl = 1; lvl<=lvlref;++lvl){
    vector<map<size_t,size_t>> midvertex(index);
    for (size_t i=0;i<face_N;++i)
    {
      v0 = faces[i].v0;
      v1 = faces[i].v1;
      v2 = faces[i].v2;
      //01
      if (midvertex[v0].find(v1)==midvertex[v0].end()){
	v.x() = (vertex[v0].x()+vertex[v1].x())/2.0;
	v.y() = (vertex[v0].y()+vertex[v1].y())/2.0;
	v.k() = ksq(v.x(),v.y());
	vertex.push_back(v);
	midvertex[v0][v1] = index;
	midvertex[v1][v0] = index;
	v01 = index;
	index++;
      }
      else{
	v01 = midvertex[v0].at(v1);
      }
      //02
      if (midvertex[v0].find(v2)==midvertex[v0].end()){
	v.x() = (vertex[v0].x()+vertex[v2].x())/2.0;
	v.y() = (vertex[v0].y()+vertex[v2].y())/2.0;
	v.k() = ksq(v.x(),v.y());
	vertex.push_back(v);
	midvertex[v0][v2] = index;
	midvertex[v2][v0] = index;
	v02 = index;
	index++;
      }
      else{
	v02 = midvertex[v0].at(v2);
      }
      //12
      if (midvertex[v1].find(v2)==midvertex[v1].end()){
	v.x() = (vertex[v1].x()+vertex[v2].x())/2.0;
	v.y() = (vertex[v1].y()+vertex[v2].y())/2.0;
	v.k() = ksq(v.x(),v.y());
	vertex.push_back(v);
	midvertex[v1][v2] = index;
	midvertex[v2][v1] = index;
	v12 = index;
	index++;
      }
      else{
	v12 = midvertex[v1].at(v2);
      }
      faces[i].v0 = v0;
      faces[i].v1 = v01;
      faces[i].v2 = v02;
      faces.push_back({v1,v01,v12});
      faces.push_back({v2,v12,v02});
      faces.push_back({v01,v02,v12});
    }
    face_N = faces.size();
    //cout<< faces.size()<< ' ' << vertex.size()<<endl;
  }
  //cout<< faces.size()<< ' ' << vertex.size()<<endl;
  //exit(0);
  
  vector<double> corners(6, 0.0);
  ELEMENTS::Triangle face_element;
  vector< vector< double > > local_stiff;
  vector< vector< double > > local_ksq;
  vector< vector< double > > local_mass;
  
  for (size_t i=0;i<face_N;++i){
    v0 = faces[i].v0;
    v1 = faces[i].v1;
    v2 = faces[i].v2;
    // array corners contains the x- and y-coordinates of the
    // triangle corners in the order x0, y0, x1, y1, x2, y2
    corners[0] = vertex[v0].x(); 
    corners[1] = vertex[v0].y();
    corners[2] = vertex[v1].x(); 
    corners[3] = vertex[v1].y();
    corners[4] = vertex[v2].x(); 
    corners[5] = vertex[v2].y();

    // pass the corners to the finite element
    face_element(corners);

    local_stiff = face_element.integrate(grad(v_()) * grad(w_()));
    local_ksq = face_element.integrate(func<double>(ksq) * v_() * w_());
    local_mass = face_element.integrate(v_() * w_());
    
    //build neighbour
    //refLvl = 0
    //00
    if (vertex[v0](v0)){
      matrix m;
      m.stiffness = local_stiff[0][0]-local_ksq[0][0];
      m.mass = local_mass[0][0];
      vertex[v0][v0] = m;
    }
    else{
      vertex[v0][v0].stiffness += local_stiff[0][0]-local_ksq[0][0];
      vertex[v0][v0].mass += local_mass[0][0];
    }
    //01
    if (vertex[v0](v1)){
      matrix m;
      m.stiffness = local_stiff[0][1]-local_ksq[0][1];
      m.mass = local_mass[0][1];
      vertex[v0][v1] = m;
    }
    else{
      vertex[v0][v1].stiffness += local_stiff[0][1]-local_ksq[0][1];
      vertex[v0][v1].mass += local_mass[0][1];
    }
    //02
    if (vertex[v0](v2)){
      matrix m;
      m.stiffness = local_stiff[0][2]-local_ksq[0][2];
      m.mass = local_mass[0][2];
      vertex[v0][v2] = m;
    }
    else{
      vertex[v0][v2].stiffness += local_stiff[0][2]-local_ksq[0][2];
      vertex[v0][v2].mass += local_mass[0][2];
    }
    
    //10
    if (vertex[v1](v0)){
      matrix m;
      m.stiffness = local_stiff[1][0]-local_ksq[1][0];
      m.mass = local_mass[1][0];
      vertex[v1][v0] = m;
    }
    else{
      vertex[v1][v0].stiffness += local_stiff[1][0]-local_ksq[1][0];
      vertex[v1][v0].mass += local_mass[1][0];
    }
    //11
    if (vertex[v1](v1)){
      matrix m;
      m.stiffness = local_stiff[1][1]-local_ksq[1][1];
      m.mass = local_mass[1][1];
      vertex[v1][v1] = m;
    }
    else{
      vertex[v1][v1].stiffness += local_stiff[1][1]-local_ksq[1][1];
      vertex[v1][v1].mass += local_mass[1][1];
    }
    //12
    if (vertex[v1](v2)){
      matrix m;
      m.stiffness = local_stiff[1][2]-local_ksq[1][2];
      m.mass = local_mass[1][2];
      vertex[v1][v2] = m;
    }
    else{
      vertex[v1][v2].stiffness += local_stiff[1][2]-local_ksq[1][2];
      vertex[v1][v2].mass += local_mass[1][2];
    }
    
    //20
    if (vertex[v2](v0)){
      matrix m;
      m.stiffness = local_stiff[2][0]-local_ksq[2][0];
      m.mass = local_mass[2][0];
      vertex[v2][v0] = m;
    }
    else{
      vertex[v2][v0].stiffness += local_stiff[2][0]-local_ksq[2][0];
      vertex[v2][v0].mass += local_mass[2][0];
    }
    //21
    if (vertex[v2](v1)){
      matrix m;
      m.stiffness = local_stiff[2][1]-local_ksq[2][1];
      m.mass = local_mass[2][1];
      vertex[v2][v1] = m;
    }
    else{
      vertex[v2][v1].stiffness += local_stiff[2][1]-local_ksq[2][1];
      vertex[v2][v1].mass += local_mass[2][1];
    }
    //22
    if (vertex[v2](v2)){
      matrix m;
      m.stiffness = local_stiff[2][2]-local_ksq[2][2];
      m.mass = local_mass[2][2];
      vertex[v2][v2] = m;
    }
    else{
      vertex[v2][v2].stiffness += local_stiff[2][2]-local_ksq[2][2];
      vertex[v2][v2].mass += local_mass[2][2];
    }

  }
  input.close();

  const string outputksq= "output/ksq.txt";
  const string outputA = "output/A.txt";
  const string outputM = "output/M.txt";
  const string outputEig = "output/eigenmode.txt";
  const string outputlambda = "output/lambda.txt";
  
  fstream output;
  //Write out k square value.
  output.open(outputksq,fstream::out);
  for (size_t i=0;i<vertex.size();++i)
    output<<setw(13)<<vertex[i].x()<<setw(13)<<vertex[i].y()<<setw(13)<<vertex[i].k()<<endl;
  output.close();
  
  //Write out Stiffness matrix
  output.open(outputA,fstream::out);
  for (size_t i=0;i<vertex.size();++i)
    for (auto itr = vertex[i].b();itr!= vertex[i].e();++itr)
      output<<setw(13)<<i<<setw(13)<<itr->first<<setw(13)<<itr->second.stiffness<<endl;
  output.close();
  
  //Write out Mass matrix
  output.open(outputM,fstream::out);
  for (size_t i=0;i<vertex.size();++i)
    for (auto itr = vertex[i].b();itr!= vertex[i].e();++itr)
      output<<setw(13)<<i<<setw(13)<<itr->first<<setw(13)<<itr->second.mass<<endl;
  output.close();
  
  double eigenvalue;
  vector<double> eigenvector;
  eigenvector = InversePowerItration(vertex,eigenvalue);
  output.open(outputEig,fstream::out);
  for (size_t i=0;i<vertex.size();++i){
    output<<setw(13)<<vertex[i].x()<<setw(13)<<vertex[i].y()<<setw(13)<<eigenvector[i]<<endl;
  }
  output.close();
  
  output.open(outputlambda,fstream::out);
  output << eigenvalue<<endl;
  output.close();
  
}

