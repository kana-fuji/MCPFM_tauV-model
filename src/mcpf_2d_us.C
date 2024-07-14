//////////////////////////////////////////////////////////////////////////////////////////
//Multicellular Systems using Phase-Field Models with anti-adhesive molecules in 2D
//Date   : since 2019/08/27  update 2022/03/01
//Auther : Kana Fuji
//////////////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>  // pf4hf
#include <string.h>  // pf4hf

#include <iomanip>   // setprecision()
#include <iostream>
#include <fstream>
#include <sstream>

#include <cmath>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#include "SFMT.h"
#include "vector3.h"

using namespace std;

//Mersenne twister---------------------------//
sfmt_t sfmt;
void my_srand( uint32_t seed ) {
    sfmt_init_gen_rand( &sfmt, seed );
}
float MT_rand(){
  return sfmt_genrand_res53(&sfmt);
}
//Mersenne twister---------------------------//

//-----------------------------------------------------------------------//
// class
//-----------------------------------------------------------------------//

class sparameter{
public:

  int xmax, ymax;            // system size
  float dx, dy;              //

  int tmax;                  // time step
  float dt;                  //
  float sampling;            // sampling rate 

  // parameters of Cells
  int xmax_cell, ymax_cell;  // system size of cells
  int init_num_cells;        // initial number of cells
  int MAX_num_cells;         // maximum number of cells
  int num_ctypes;
  float *rx;                 // position of initial cells
  float *ry;                 // position of initial cells
  float *init_r_u;           // radius of initial cells

  float D;                   // diffusion coefficient
  float V;                   // maximum volume of the cells
  float tau;                 // relaxation time of the cells
  float v_d;                 // volume cell division
  float ep_d;                //
  float alpha;               // volume consevation
  float beta;                // excluded volume
  float gamma;               // cell-cell adhesion
  float eta;                 // cell-cell adhesion
  float tol_u;               // reshaping

  float gamma_curv_u;        // curevature effects
  float gamma_curv_s;        // curevature effects
  float gamma_curv_c;        // curevature effects

  //CFG parameters
  float rho0;
  float rhoe;
  float dts;     
  float mu;
  float ls;
  float sigma;

  // lumen
  float xi0;           // lumen pressure
  float c_s;           // feedback effect of lumen pressure
  float beta_s;        // excluded volume between lumens & cells
  float tau_s;         // relaxation time of lumen
  float Ds;            //
  float eta_s, gamma_s;// cell-lumen adhesion
  float r_s;           //
  float tol_s;         // reshaping
  float alpha_s;       // local volume consevasion
  float V_s;           // local volume consevasion
  float v_t;           //
  float p_st;          //

  float init_theta;    // angle of the initial division plane

  float tau_V;         // relaxation time of the maximum cell volume
  float noise_tauV;    // noise of tauV
  float alpha_V;
  float td;            // division condtion
  float noise_td;      // noise of td
  float rseed;         // seed for noise of tauV

  float tau_c;
  float Dc;
  float alpha_c;
  float beta_cu,beta_cs;
  float eta_cu,eta_cs,gamma_c;
  float Vsys;
  float tol_c;         // reshaping
  float xi_c;          // pressure of medium

  float beta_w;        // excluded volume of wall
  float eta_w;         // adhesion on the wall

  float kappa_n;       // interaction with other cell polarities
  float lambda_n;      // interaction with lumens

  //anti-adhesive molecules
  float tau_p;
  float D_p;
  float alpha_p;
  float V_p;
  float C_p;
  float dwell;
  float rp;
  float kp;
  float p_th;
  float l_anti;

  float D_pol,tau_pol;
  float gamma_pol;

  string Dir;           // output directry 

  sparameter();
  ~sparameter();

  void paraminput(char *argv[]);
} param ;

sparameter::sparameter(){
  xmax=0; ymax=0;
  dx=0.0f; dy=0.0f;
  tmax=0;
  dt=0.0f;
  sampling=1.0f;

  D=0.0f;
  V=0.0f;
  tau=0.0f;
  xmax_cell=0; ymax_cell=0;
  init_num_cells=0;
  num_ctypes=0;
  v_d=0.0f;
  ep_d=0.0f;
  beta=0.0f;
  tol_u=0.0f;

  gamma_curv_u=0.0f;
  gamma_curv_s=0.0f;
  gamma_curv_c=0.0f;

  tau_V=0.0f;
  noise_tauV=0.0f;
  alpha_V=0.0f;
  rseed=0.0f;
  xi0=0.0f;
  c_s=0.0f;
  tol_s=0.0f;
  alpha_s=0.0f;
  V_s=0.0f;
  v_t=0.0f;
  p_st=0.0f;

  tau_c=0.0f;
  Dc=0.0f;
  alpha_c=0.0f;
  beta_cu=0.0f,beta_cs=0.0f;
  eta_cu=0.0f,eta_cs=0.0f,gamma_c=0.0f;
  Vsys=0.0f;
  xi_c=0.0f;
  tol_c=0.0f;

  beta_w=0.0f;
  eta_w=0.0f;

  kappa_n=0.0f;
  lambda_n=0.0f;

  tau_p=0.0f;
  D_p=0.0f;
  alpha_p=0.0f;
  V_p=0.0f;
  C_p=0.0f;
  rp=0.0f;
  dwell=0.0f;
  kp=0.0f;
  p_th=0.0f;
  l_anti=0.0f;
}

inline sparameter::~sparameter(){
  delete[] rx;
  delete[] ry;
  delete[] init_r_u;
}


void sparameter::paraminput(char *argv[]){

  //system
  //cin>>xmax>>dx>>ymax>>dy;
  xmax=2048; ymax=2048;
  dx=0.02; dy=0.02;
  cout<<"xmax="<<xmax<<" dx="<<dx<<" ymax="<<ymax<<" dy="<<dy<<endl;
  int tend=0;
  //cin>>tend>>dt>>sampling;
  tend=100000; dt=0.02; sampling=10;
  tmax=int(tend/dt);
  cout<<"tend="<<tend<<" tmax="<<tmax<<" dt="<<dt<<" sampling="<<sampling<<endl;

  //cin>>xmax_cell>>ymax_cell;
  xmax_cell=200; ymax_cell=200;
  cout<<"xmax_cell="<<xmax_cell<<" ymax_cell="<<ymax_cell<<endl;

  //cell
  //cin>>num_ctypes;
  //cin>>MAX_num_cells;
  //cin>>init_num_cells;
  num_ctypes=1;
  MAX_num_cells=300;
  init_num_cells=atoi(argv[2]);
  cout<<"num_ctypes="<<num_ctypes<<" Max_num_cells="<<MAX_num_cells<<" init_num_cells="<<init_num_cells<<endl;
  //cout<<argv[2]<<endl;

  rx= new float[init_num_cells];  
  ry= new float[init_num_cells];  
  init_r_u= new float[init_num_cells];
  // for(int n=0;n<init_num_cells;n++){
  //   cin>>rx[n]>>ry[n]>>init_r_u[n];
  // }
  
  string fname;
  fname="init_cells_"+to_string(init_num_cells);
  cout<<fname<<endl;
  ifstream fin(fname.c_str()); 
  if(!fin.is_open()){cerr<<"ERROR:Could not open input file "<<endl;exit(8);}
  for(int n=0;n<init_num_cells;n++){ fin>>rx[n]>>ry[n]>>init_r_u[n]; }
  fin.close();

  //cin>>V>>v_d;
  V=3.0; v_d=0.1;//atof(argv[8]);
  cout<<"V="<<V<<" v_d="<<v_d<<endl;

  //cin>>tau>>D;
  tau=1.0; D=0.001;
  cout<<"tau="<<tau<<" D="<<D<<endl;
  //20200807 reshaping & curvature effects
  //cin>>tol_u>>gamma_curv_u;
  tol_u=0.000; gamma_curv_u=0.00;
  cout<<"tol_u="<<tol_u<<" gamma_curv="<<gamma_curv_u<<endl;
  //cin>>alpha>>beta>>gamma>>eta;
  alpha=1.0f; beta=1.0; gamma=0.01; eta=0.008;
  cout<<"alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<" eta="<<eta<<endl;


  //division
  //cin>>ep_d;
  //cin>>init_theta;
  ep_d=0.1; init_theta=0.0;
  cout<<"ep_d="<<ep_d<<" init_theta="<<init_theta<<endl;

  //CFG
  //cin>>rho0>>rhoe;
  rho0=0.01; rhoe=5.0;
  cout<<"rho0="<<rho0<<" rhoe="<<rhoe<<endl;
  //cin>>dts;
  dts=1.0;
  cout<<"dts="<<dts<<endl;
  //cin>>sigma>>mu>>ls;
  sigma=0.001;mu=1.0;ls=0.0;
  cout<<"sigma="<<sigma<<" mu="<<mu<<" ls="<<ls<<endl;

  //Vmax
  //cin>>tau_V>>noise_tauV;
  //cin>>alpha_V;
  tau_V=atof(argv[3]); noise_tauV=0.25; alpha_V=1.0;
  cout<<"tau_V="<<tau_V<<" noise_tauV="<<noise_tauV<<" alpha_V="<<alpha_V<<endl;
  //cin>>td>>noise_td;
  td=1;//atof(argv[3]);//not used
  noise_td=0.10;
  cout<<"td="<<td<<" noise_td="<<noise_td<<endl;

  //lumen
  //cin>>tau_s>>Ds;
  tau_s=1.0; Ds=0.001;
  cout<<"tau_s="<<tau_s<<" Ds="<<Ds<<endl;
  //cin>>tol_s>>gamma_curv_s;
  tol_s=0.000; gamma_curv_s=0.00;
  cout<<"tol_s="<<tol_s<<" gamma_curv_s="<<gamma_curv_s<<endl;
  //cin>>beta_s;
  beta_s=1.0;
  //cin>>xi0>>c_s;
  xi0=atof(argv[4]);
  c_s=0.0;
  cout<<"beta_s="<<beta_s<<" xi0="<<xi0<<" c_s="<<c_s<<endl;
  //cin>>gamma_s>>eta_s;
  gamma_s=0.0; eta_s=0.000;
  cout<<"gamma_s="<<gamma_s<<" eta_s="<<eta_s<<endl;
  //cin>>r_s;
  r_s=0.8;//r_s=0.5;
  cout<<"rs="<<r_s<<endl;
  //cin>>alpha_s>>V_s;
  alpha_s=0.0; V_s=1.0;
  cout<<"alpha_s="<<alpha_s<<" V_s="<<V_s<<endl;
  //cin>>v_t>>p_st;
  v_t=0.3; p_st=0.4;
  cout<<"v_t="<<v_t<<" p_st="<<p_st<<endl;

  //20200526 ecm //not used
  //cin>>tau_c>>Dc;
  tau_c=1.0; Dc=0.001;
  cout<<"tau_c="<<tau_c<<" Dc="<<Dc<<endl;
  //cin>>tol_c>>gamma_curv_c;
  tol_c=0.000; gamma_curv_c=0.00;
  cout<<"tol_c="<<tol_s<<" gamma_curv_c="<<gamma_curv_c<<endl;
  //cin>>alpha_c>>beta_cu>>beta_cs>>xi_c;
  alpha_c=0.001; beta_cu=1.0; beta_cs=1.0; xi_c=0.000;
  cout<<"alpha_c="<<alpha_c<<" beta_cu="<<beta_cu<<" beta_cs="<<beta_cs<<" xi_c="<<xi_c<<endl;
  //cin>>eta_cu>>eta_cs>>gamma_c;
  eta_cu=0.000; eta_cs=0.00; gamma_c=0.00;
  cout<<"eta_cu="<<eta_cu<<" eta_cs="<<eta_cs<<" gamma_c="<<gamma_c<<endl;
  Vsys=xmax*dx*ymax*dy;

  // //wall
  // cin>>beta_w>>eta_w;
  // cout<<"beta_w="<<beta_w<<" eta_w="<<eta_w<<endl;

  // //cell polarity
  // cin>>kappa_n>>lambda_n;
  // cout<<"kappa_n="<<kappa_n<<" lambda_n="<<lambda_n<<endl;

  //random number
  //cin>>rseed;
  //rseed=atof(argv[6]);
  rseed=1;
  cout<<"rseed="<<rseed<<endl;

  //Ds=atof(argv[30]);
  Dir="DATA/"+string(argv[1]);
  cout<<Dir<<endl;
  //cin>>Dir;

  fname=Dir+"/param.txt";
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file number of cell data"<<endl;exit(8);}
  fout<<sampling<<endl;
  fout.close();
}

class Commonfunction{
public:
  float h(float u);
  float kdelta(int& i,int& j);
};

inline float Commonfunction::h(float u){
  return u*u*(3.0f-2.0f*u);
}

inline float Commonfunction::kdelta(int& i,int& j){
  if(i==j)
    return 1.0f;
  else
    return 0.0f;
}



class Field:public Commonfunction{
protected:
  int m_xmax, m_ymax;
  float **m_f;
public:
  Field();
  Field(const int xmax, const int ymax);
  ~Field();
public:
  int Get_xmax(), Get_ymax();
  float Get_value(const int xmax, const int ymax);
  void Set_value(const int xmax, const int ymax, const float value);

  void init_field();

  float laplacian4(const int i, const int j);
  float laplacian8(const int i, const int j);
  float laplacian8_h(const int i, const int j);
  vector3 nabla_h(const int i, const int j);
  vector3 nabla(const int i, const int j);

  float volume();

  void outputtofile(const string fname);
  void outputtofile2(const string fname);
};
inline Field::Field(){
  m_xmax=0; m_ymax=0;
}
inline Field::Field(const int xmax, const int ymax){
  m_xmax=xmax; m_ymax=ymax;

  m_f = new float*[xmax];
  for (int i=0;i<xmax;i++) {
    m_f[i] = new float[ymax];
  }
}
inline Field::~Field(){
  for (int i=0; i<m_xmax;i++) {
    delete[] m_f[i];
  }
  delete[] m_f;
}

inline int Field::Get_xmax(){return m_xmax;}
inline int Field::Get_ymax(){return m_ymax;}
inline float Field::Get_value(const int i, const int j){return m_f[i][j];}
inline void Field::Set_value(const int i, const int j, const float v){m_f[i][j]=v;}

inline void Field::init_field(){
#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
	m_f[i][j]=0.0f;
    }
  }
}

inline float Field::laplacian4(const int i,const int j){
  return 
    (m_f[i][j-1]-2.0f*m_f[i][j]+m_f[i][j+1])/(param.dx*param.dx)
    +(m_f[i-1][j]-2.0f*m_f[i][j]+m_f[i+1][j])/(param.dy*param.dy)
    ;
}

inline float Field::laplacian8(const int i,const int j){
  return (
	  (m_f[i][j-1]-2.0f*m_f[i][j]+m_f[i][j+1])
	  /(param.dx*param.dx)
	  +(m_f[i-1][j]-2.0f*m_f[i][j]+m_f[i+1][j])
	  /(param.dy*param.dy)
	  +(0.5f*
	    (m_f[i-1][j-1]+m_f[i-1][j+1]+m_f[i+1][j-1]+m_f[i+1][j+1])-2.0f*m_f[i][j])
	  /(param.dx*param.dy)
	  )*0.5f;
}

inline float Field::laplacian8_h(const int i,const int j){
  return (
	  (h(m_f[i][j-1])-2.0f*h(m_f[i][j])+h(m_f[i][j+1]))
	  /(param.dx*param.dx)
	  +(h(m_f[i-1][j])-2.0f*h(m_f[i][j])+h(m_f[i+1][j]))
	  /(param.dy*param.dy)
	  +(0.5f*
	    (h(m_f[i-1][j-1])+h(m_f[i-1][j+1])+h(m_f[i+1][j-1])+h(m_f[i+1][j+1]))-2.0f*h(m_f[i][j]))
	  /(param.dx*param.dy)
	  )*0.5f;
}

inline vector3 Field::nabla_h(const int i, const int j){
  vector3 v;
  v.set(
	(h(m_f[i+1][j])-h(m_f[i-1][j]))*0.5f/param.dx,
	(h(m_f[i][j+1])-h(m_f[i][j-1]))*0.5f/param.dy,
	0.0f
	);
  return  v;
}

inline vector3 Field::nabla(const int i, const int j){
  vector3 v;
  v.set(
	(m_f[i+1][j]-m_f[i-1][j])*0.5f/param.dx,
	(m_f[i][j+1]-m_f[i][j-1])*0.5f/param.dy,
	0.0f
	);
  return  v;
}


inline float Field::volume(){
  float v=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      v+=m_f[i][j]*param.dx*param.dy;
  }
  return v;
}


void Field::outputtofile(const string fname){
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file "<<fname<<endl;exit(8);}
  for(int j=m_ymax-2;j>=0;j=j-2){
    for(int i=0;i<m_xmax;i=i+2){
      fout<<m_f[i][j]<<" ";
    }
    fout<<endl;
  }
  fout.close();
}


void Field::outputtofile2(const string fname){
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file "<<fname<<endl;exit(8);}
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      fout<<m_f[i][j]<<" ";
    }
    fout<<endl;
  }
  fout.close();
}


class VField:public Commonfunction{
protected:
  int m_xmax, m_ymax;
  vector3 **m_f;
public:
  VField();
  VField(const int xmax, const int ymax);
  ~VField();
public:
  void   Set_V(const int xmax, const int ymax, const vector3 f);
  vector3 Get_V(const int xmax, const int ymax);

  void outputtofile(const string fname);
};
inline VField::VField(){
  m_xmax=0; m_ymax=0;
}
inline VField::VField(const int xmax, const int ymax){
  m_xmax=xmax; m_ymax=ymax;

  m_f = new vector3*[xmax];
  for (int i=0;i<xmax;i++) {
    m_f[i] = new vector3[ymax];
  }
}
inline VField::~VField(){
  for (int i=0; i<m_xmax;i++) {
    delete[] m_f[i];
  }
  delete[] m_f;
}
inline void VField::Set_V(const int i, const int j, const vector3 f){m_f[i][j]=f;}
inline vector3 VField::Get_V(const int i, const int j){return  m_f[i][j];}

void VField::outputtofile(const string fname){
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file VField data"<<endl;exit(8);}
  for(int i=0;i<m_xmax;i=i+20){
    for(int j=0;j<m_ymax;j=j+20)
      fout<<i*param.dx<<" "<<j*param.dy<<" "<<m_f[i][j].x<<" "<<m_f[i][j].y<<endl;
    fout<<endl;
  }
  fout.close();
}


class PhaseField:public Commonfunction{
protected:
  int m_xmax, m_ymax;
  float m_dx, m_dy;
  float **m_pf, **m_pf_temp;

  float m_dt; 
  float m_tau;                //relaxation time
  float m_D;                  //diffusion coefficient

  float m_volume;             //volume

  float m_tol;                //reshaping step
  float m_gamma_curv;         //curvature effects

public:
  PhaseField();
  PhaseField(const int xmax, const int ymax);
  ~PhaseField();

public:
  int Get_xmax(), Get_ymax();
  float Get_value(const int xmax, const int ymax);
  float Get_value_temp(const int xmax, const int ymax);
  float Get_tau(), Get_D();
  float Get_tol();
  float Get_gamma_curv();
  float Get_volume();
  void Set_volume(const float volume);
  void Set_value(const int xmax, const int ymax, const float v);
  void Set_value_temp(const int xmax, const int ymax, const float v);

  void boundary_periodic();
  void boundary_neumann();
  void boundary_dirichlet(const float v);

  void init_field();
  void update();

  float laplacian4(const int i,const int j);
  float laplacian8(const int i,const int j);
  float laplacian8_h(const int i, const int j);
  vector3 nabla(const int i, const int j);
  vector3 nabla_h(const int i, const int j);

  float volume_h();
  float volume();
  float surface_h();
  float surface();

  void outputtofile(const string fname);

  vector3 normal(const int i,const int j);
  void intermediate_step();
}pf;

inline PhaseField::PhaseField(){
  m_xmax=0; m_ymax=0;
}

inline PhaseField::PhaseField(const int xmax, const int ymax){
  m_xmax=xmax; m_ymax=ymax;

  m_pf = new float*[xmax];
  m_pf_temp = new float*[xmax];
  for (int i=0;i<xmax;i++) {
    m_pf[i] = new float[ymax];
    m_pf_temp[i] = new float[ymax];
  }
}

inline PhaseField::~PhaseField(){
  for (int i=0; i<m_xmax;i++) {
    delete[] m_pf[i];
    delete[] m_pf_temp[i];
  }
  delete[] m_pf;
  delete[] m_pf_temp;
}

inline int PhaseField::Get_xmax(){return m_xmax;}
inline int PhaseField::Get_ymax(){return m_ymax;}
inline float PhaseField::Get_value(const int i, const int j){return m_pf[i][j];}
inline float PhaseField::Get_value_temp(const int i, const int j){return m_pf_temp[i][j];}
inline float PhaseField::Get_volume(){return m_volume;}
inline float PhaseField::Get_tau(){return m_tau;}
inline float PhaseField::Get_D(){return m_D;}
inline float PhaseField::Get_tol(){return m_tol;}
inline float PhaseField::Get_gamma_curv(){return m_gamma_curv;}
inline void PhaseField::Set_value(const int i, const int j, const float v){m_pf[i][j]=v;}
inline void PhaseField::Set_value_temp(const int i, const int j, const float v){m_pf_temp[i][j]=v;}
inline void PhaseField::Set_volume(const float volume){m_volume=volume;}

inline void PhaseField::init_field(){
  //#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      m_pf[i][j]=0.0f;
  }
}


inline void PhaseField::update(){
  //#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      m_pf[i][j]=m_pf_temp[i][j];
      //cout<<m_pf[i][j]<<" "<<m_pf_temp[i][j]<<endl;
    }
  }
}

inline void PhaseField::boundary_periodic(){ 
  for(int i=0;i<m_xmax;i++){
    m_pf[0][i]=m_pf[m_xmax][i]; m_pf[m_xmax+1][i]=m_pf[1][i];
    m_pf[i][0]=m_pf[i][m_ymax]; m_pf[i][m_ymax+1]=m_pf[i][1];
  }
}

inline void PhaseField::boundary_neumann(){ 
  for(int i=0;i<m_xmax;i++){
    m_pf[0][i]=m_pf[1][i]; m_pf[m_xmax+1][i]=m_pf[m_xmax][i];
    m_pf[i][0]=m_pf[i][1]; m_pf[i][m_ymax+1]=m_pf[i][m_ymax];
  }
}

inline void PhaseField::boundary_dirichlet(const float v){ 
  for(int i=0;i<m_xmax;i++){
    m_pf[0][i]=v; m_pf[m_xmax-1][i]=v;
    m_pf[i][0]=v; m_pf[i][m_ymax-1]=v;
    m_pf_temp[0][i]=v; m_pf_temp[m_xmax-1][i]=v;
    m_pf_temp[i][0]=v; m_pf_temp[i][m_ymax-1]=v;
  }
}

inline float PhaseField::laplacian4(const int i,const int j){
  return 
    (m_pf[i][j-1]-2.0f*m_pf[i][j]+m_pf[i][j+1])/(param.dx*param.dx)
    +(m_pf[i-1][j]-2.0f*m_pf[i][j]+m_pf[i+1][j])/(param.dy*param.dy)
    ;
}

inline float PhaseField::laplacian8(const int i,const int j){
  return (
	  (m_pf[i][j-1]-2.0f*m_pf[i][j]+m_pf[i][j+1])
	  /(param.dx*param.dx)
	  +(m_pf[i-1][j]-2.0f*m_pf[i][j]+m_pf[i+1][j])
	  /(param.dy*param.dy)
	  +(0.5f*
	    (m_pf[i-1][j-1]+m_pf[i-1][j+1]+m_pf[i+1][j-1]+m_pf[i+1][j+1])-2.0f*m_pf[i][j])
	  /(param.dx*param.dy)
	  )*0.5f;
}

inline float PhaseField::laplacian8_h(const int i,const int j){
  return (
	  (h(m_pf[i][j-1])-2.0f*h(m_pf[i][j])+h(m_pf[i][j+1]))
	  /(param.dx*param.dx)
	  +(h(m_pf[i-1][j])-2.0f*h(m_pf[i][j])+h(m_pf[i+1][j]))
	  /(param.dy*param.dy)
	  +(0.5f*
	    (h(m_pf[i-1][j-1])+h(m_pf[i-1][j+1])+h(m_pf[i+1][j-1])+h(m_pf[i+1][j+1]))-2.0f*h(m_pf[i][j]))
	  /(param.dx*param.dy)
	  )*0.5f;
}

inline vector3 PhaseField::nabla(const int i, const int j){
  vector3 v;
  v.set(
	(m_pf[i+1][j]-m_pf[i-1][j])*0.5f/param.dx,
	(m_pf[i][j+1]-m_pf[i][j-1])*0.5f/param.dy,
	0.0f
	);
  return  v;
}

inline vector3 PhaseField::nabla_h(const int i, const int j){
  vector3 v;
  v.set(
	(h(m_pf[i+1][j])-h(m_pf[i-1][j]))*0.5f/param.dx,
	(h(m_pf[i][j+1])-h(m_pf[i][j-1]))*0.5f/param.dy,
	0.0f
	);
  return  v;
}

inline float PhaseField::volume_h(){
  float v=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      v+=h(m_pf[i][j])*param.dx*param.dy;
  }
  return v;
}

inline float PhaseField::volume(){
  float v=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      v+=m_pf[i][j]*param.dx*param.dy;
  }
  return v;
}

inline float PhaseField::surface(){
  float s=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      s+=m_pf[i][j]*(1.0f-m_pf[i][j])*param.dx*param.dy;
  }
  return s;
}

inline float PhaseField::surface_h(){
  float s=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      s+=h(m_pf[i][j])*(1.0f-h(m_pf[i][j]))*param.dx*param.dy;
  }
  return s;
}


void PhaseField::outputtofile(const string fname){
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file "<<fname<<endl;exit(8);}
  for(int j=m_ymax-2;j>=0;j=j-2){
    for(int i=0;i<m_xmax;i=i+2){
      fout<<m_pf[i][j]<<" ";
    }
    fout<<endl;
  }
  fout.close();
}

vector3 PhaseField::normal(const int i,const int j){
  vector3 v=nabla(i,j);
  float vabs=(float)(v.abs());
  if(vabs>0.00001f){
    v/=vabs;
  }
  else
    v.set(0.0f,0.0f,0.0f);
  return v;
}

void PhaseField::intermediate_step(){

  //cout<<"intermediate step ... ";
  float dtau=0.04f;//<=dx*dx/sqrt(D)/4 (Olsson 2005)

  vector3 v(0.0f,0.0f,0.0f);
  VField* n = new VField(m_xmax,m_ymax); //normal vector
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++) n->Set_V(i,j,normal(i,j));
  }
  for(int i=0;i<m_xmax-1;i++){
    n->Set_V(i,0,n->Get_V(i,m_ymax-2));
    n->Set_V(0,i+1,n->Get_V(m_xmax-2,i+1));
    n->Set_V(i+1,m_ymax-1,n->Get_V(i+1,1));
    n->Set_V(m_xmax-1,i,n->Get_V(1,i));
  }
  
  //cout<<m_volume<<" ";

  int count=0;
  float dev=m_tol*dtau+1;
  while(dev>m_tol*dtau){
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]= m_pf[i][j]+dtau*m_D*laplacian8(i,j)
	  -dtau*sqrt(2.0f*m_D) *
	  (
	   (m_pf[i+1][j]*(1.0f-m_pf[i+1][j])*n->Get_V(i+1,j).x
	    -m_pf[i-1][j]*(1.0f-m_pf[i-1][j])*n->Get_V(i-1,j).x)*0.5f/param.dx
	   +(m_pf[i][j+1]*(1.0f-m_pf[i][j+1])*n->Get_V(i,j+1).y
	     -m_pf[i][j-1]*(1.0f-m_pf[i][j-1])*n->Get_V(i,j-1).y)*0.5f/param.dy
	   )
	  ;
      }
    }

    dev=0.0f;
    for(int i=0;i<m_xmax;i++){
      for(int j=0;j<m_ymax;j++) dev+=fabsf(m_pf_temp[i][j]-m_pf[i][j])*param.dx*param.dy;
    }

    update();
    count++;
    //cout<<count<<" "<<dev<<" ";
  }
  delete n;
  m_volume=volume();
  //cout<<m_volume<<" ";
  //cout<<"reshaping count="<<count<<endl;

}



class Lumen;
class ECM;
class Wall;
class Interaction;
class VInteraction;

class Cell: public PhaseField{ //2D
private:
  int   m_imin, m_jmin;       //position in field
  float m_Gx_cell,m_Gy_cell;  //center in field

  float m_max_volume;         //maximum volume of cells
  float m_time;               //cell time
  float m_td;                 //division condition

  vector3 m_pol;               //cell polarity
  vector3 m_pol_tmp;           //cell polarity

  float m_alpha, m_beta, m_gamma, m_eta;
  float m_alpha_V, m_tau_V;

  float m_p_th, m_l_anti;

public:
  Cell();
  Cell(const int xmax_cell, const int ymax_cell);
  Cell(const Cell& cell);

public:
  int Get_imin(), Get_jmin();
  float Get_Gx_cell(), Get_Gy_cell();
  float Get_max_volume();
  float Get_time();
  float Get_td();
  vector3 Get_pol();
  vector3 Get_pol_tmp();

  float Get_alpha();
  float Get_beta();
  float Get_gamma();
  float Get_eta();
  float Get_tau_V();
  float Get_alpha_V();

  void Set_imin(const int imax), Set_jmin(const int jmax);
  void Set_Gx_cell(const float Gx_cell), Set_Gy_cell(const float Gy_cell);
  void Set_max_volume(const float max_volume_cell);
  void Set_pol(const vector3 pol);
  void Set_pol_tmp(const vector3 pol);
  void Set_gamma(const float gamma);
  void Set_eta(const float eta);

  void init_cell(const int num);
  void input_cell(const int num); 
  void output_cell(const int num);
  void update_cell();
  void seed_lumen(float* Pc,Lumen* lumen);
  void CoM();

  void time_evolution(Interaction *phi);
  float f(int& i,int& j,Interaction* field);
  float f(int& i,int& j,Interaction* field,int& pcx);
  void time_evolution(Lumen *lumen,Interaction *phi);
  float f(int& i,int& j,Interaction* field,Lumen* lumen);
  void time_evolution(Lumen *lumen,Interaction *phi,ECM *ecm);
  float f(int& i,int& j,Interaction* field,Lumen* lumen,ECM *ecm);
  void time_evolution(Interaction *phi,ECM *ecm);
  float f(int& i,int& j,Interaction* field,ECM *ecm);
  void time_evolution(Lumen *lumen,Interaction *phi,ECM *ecm,Interaction *pcxs);
  void time_evolution_debug(int t,Interaction *phi,Interaction *pcxs);
  float f(int& i,int& j,Interaction* field,Lumen* lumen,ECM *ecm,Interaction *pcxs);
  void time_evolution(Lumen *lumen,Interaction *phi,Interaction *pcxs);
  float f(int& i,int& j,Interaction* field,Lumen* lumen,int& pcx);
  // void time_evolution(Lumen *lumen,Interaction *phi,Wall *wall,VInteraction* polfield);
  // float f(int& i,int& j,Interaction* field,Lumen* lumen,Wall* wall,VInteraction* polfield);

  void time_evolution_V();
  void add_time();

  void copy_param(Cell*& cell);
  // void division(Cell *&mcell,Lumen* lumen,Interaction* e_eta,int num_cells);
  // //CFG
  // void time_evolution_poles(Interaction *e_eta,float (&r1)[2],float (&r2)[2],int num_cells);
  // float chi(float (&r)[2],float (&r1)[2],float (&r2)[2]);
  // float g(float (&r)[2],float (&r1)[2],float (&r2)[2]);
  // float H(float u);
  // float G(float e_eta);
  // float distance(float (&r1)[2],float (&r2)[2]);
}cell ;

inline Cell::Cell(){
  m_imin=0; m_jmin=0;
  m_Gx_cell=0.0f; m_Gy_cell=0.0f;
  m_pol.set(0.0f,0.0f,0.0f);
  m_pol_tmp.set(0.0f,0.0f,0.0f);
}

inline Cell::Cell(const int xmax_cell, const int ymax_cell):PhaseField(xmax_cell,ymax_cell){}

Cell::Cell(const Cell& cell){

  for (int i=0; i<m_xmax;i++) {
    delete[] m_pf[i];
    delete[] m_pf_temp[i];
  }
  delete[] m_pf;
  delete[] m_pf_temp;

  m_pf = new float*[m_xmax];
  m_pf_temp = new float*[m_xmax];
  for (int i=0;i<m_xmax;i++) {
    m_pf[i] = new float[m_ymax];
    m_pf_temp[i] = new float[m_ymax];
  }
}

inline int Cell::Get_imin(){return m_imin;}
inline int Cell::Get_jmin(){return m_jmin;}
inline float Cell::Get_Gx_cell(){return m_Gx_cell;}
inline float Cell::Get_Gy_cell(){return m_Gy_cell;}
inline float Cell::Get_max_volume(){return m_max_volume;}
inline float Cell::Get_time(){return m_time;}
inline float Cell::Get_td(){return m_td;}
inline vector3 Cell::Get_pol(){return m_pol;}
inline vector3 Cell::Get_pol_tmp(){return m_pol_tmp;}

inline float Cell::Get_alpha(){return m_alpha;}
inline float Cell::Get_beta(){return m_beta;}
inline float Cell::Get_gamma(){return m_gamma;}
inline float Cell::Get_eta(){return m_eta;}
inline float Cell::Get_tau_V(){return m_tau_V;}
inline float Cell::Get_alpha_V(){return m_alpha_V;}

inline void Cell::Set_imin(const int imin){m_imin=imin;}
inline void Cell::Set_jmin(const int jmin){m_jmin=jmin;}
inline void Cell::Set_Gx_cell(const float Gx_cell){m_Gx_cell=Gx_cell;}
inline void Cell::Set_Gy_cell(const float Gy_cell){m_Gy_cell=Gy_cell;}
inline void Cell::Set_max_volume(const float max_volume_cell){m_max_volume=max_volume_cell;}
inline void Cell::Set_pol(const vector3 pol){m_pol=pol;}
inline void Cell::Set_pol_tmp(const vector3 pol){m_pol_tmp=pol;}
inline void Cell::Set_gamma(const float gamma){m_gamma=gamma;}
inline void Cell::Set_eta(const float eta){m_eta=eta;}

inline void Cell::init_cell(const int num){

  m_tau=param.tau;
  m_D  =param.D;
  m_tol=param.tol_u;

  m_alpha=param.alpha;
  m_beta =param.beta;
  m_gamma=param.gamma;
  m_eta  =param.eta;
  m_gamma_curv=param.gamma_curv_u;

  m_tau_V=param.tau_V+(MT_rand()*2.0f-1.0f)*param.noise_tauV*param.tau_V;
  //cout<<"set tauV="<<m_tau_V<<endl;
  m_alpha_V=param.alpha_V;



  m_Gx_cell=param.rx[num]; m_Gy_cell=param.ry[num];
  m_imin=m_Gx_cell/param.dx-m_xmax*0.5f;
  m_jmin=m_Gy_cell/param.dy-m_ymax*0.5f;
  cout<<"center "<<m_Gx_cell<<" "<<m_Gy_cell<<" imin,jmin "<<m_imin<<" "<<m_jmin<<endl;

  float center_x=m_xmax*0.5f*param.dx;
  float center_y=m_ymax*0.5f*param.dy;

  float r0=param.init_r_u[num];
  float r=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      r=sqrtf((param.dx*(i+0.5f)-center_x)*(param.dx*(i+0.5f)-center_x)
	      +(param.dy*(j+0.5f)-center_y)*(param.dy*(j+0.5f)-center_y));
      m_pf[i][j]=(1.0f-tanhf((r-r0)/(sqrtf(2.0f*m_D))))*0.5f;
      m_pf_temp[i][j]=0.0f;
      //cout<<m_pf[i][j]<<endl;
    }
  }

  m_volume=volume_h();
  cout<<"volume="<<m_volume<<endl;
  if(m_alpha_V>0.0f){
    m_max_volume=m_volume*1.1f;
  }
  else{
    m_max_volume=param.V;
  }

  //set cell time
  double rand1=MT_rand();
  double rand2=MT_rand();
  m_td=param.td+(rand1*2.0-1.0)*param.noise_td*param.td;
  m_time=rand2*m_td;
  cout<<num<<" r1="<<rand1<<" r2="<<rand2<<" td="<<m_td<<" time="<<m_time<<endl;

}


inline void Cell::input_cell(const int num){
  m_tau=param.tau;
  m_D  =param.D;

  m_alpha=param.alpha;
  m_beta =param.beta;
  m_gamma=param.gamma;
  m_eta  =param.eta;

  m_tau_V=param.tau_V+(MT_rand()*2.0f-1.0f)*param.noise_tauV*param.tau_V;
  cout<<"set tauV="<<m_tau_V<<endl;

  //m_tau_V  =param.tau_V;
  m_alpha_V=param.alpha_V;

  stringstream ss;
  ss<<num;

  string fname;
  fname="STATE/u_No."+ss.str()+".dat";
  ifstream fin(fname.c_str()); 
  if(!fin.is_open()){cerr<<"ERROR:Could not open input file "<<fname<<endl;exit(8);}
  fin>>m_imin>>m_jmin;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      fin>>m_pf[i][j];
  }
  fin.close();

  m_Gx_cell=(m_imin+m_xmax*0.5f)*param.dx;
  m_Gy_cell=(m_jmin+m_ymax*0.5f)*param.dy;

  m_volume=volume_h();
  if(m_alpha_V>0.0f){
    //m_max_volume=(param.V-param.v_d)*0.5f;
    m_max_volume=m_volume*1.1f;
  }
  else m_max_volume=param.V;
}



inline void Cell::output_cell(const int num){

  //set position
  stringstream ss;
  ss<<num;

  string fname;
  fname=param.Dir+"/u_No."+ss.str()+".dat";
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file "<<fname<<endl;exit(8);}
  fout<<m_imin<<" "<<m_jmin<<endl;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      fout<<fixed<<setprecision(16)<<m_pf[i][j]<<" ";
    fout<<endl;
  }
  fout.close();
}

inline void Cell::update_cell(){
  update();
  m_volume=volume_h();
  m_pol=m_pol_tmp;
}

void Cell::CoM(){

  float sumx=0.0f; float sumy=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      sumx+=h(m_pf[i][j])*param.dx*param.dy*(i+0.5f);
      sumy+=h(m_pf[i][j])*param.dx*param.dy*(j+0.5f);
      //cout<<i<<" "<<j<<" "<<sumx<<" "<<sumy<<" "<<m_pf[i][j]<<" "<<h(m_pf[i][j])<<endl;
    }
  }

  float di=sumx/volume_h()-m_xmax*0.5f;
  float dj=sumy/volume_h()-m_ymax*0.5f;
  //cout<<di<<" "<<dj<<" "<<sumx<<" "<<sumy<<" ";

  int ii,jj;
  //#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    ii=i+(int)di;
    for(int j=0;j<m_ymax;j++){
      jj=j+(int)dj;
      if(ii>0 && ii<m_xmax-1 && jj>0 && jj<m_ymax-1)
	m_pf_temp[i][j]=m_pf[i+(int)di][j+(int)dj];
      else
	m_pf_temp[i][j]=0.0f;
    }
  }
  update();

  m_imin+=(int)di; m_jmin+=(int)dj;
  m_Gx_cell+=(float)(((int)di)*param.dx);
  m_Gy_cell+=(float)(((int)dj)*param.dy);

}


void Cell::time_evolution_V(){
  //saturation
  m_max_volume+=param.dt/m_tau_V*m_alpha_V*(param.V-m_max_volume);
  //linear
  //if(m_max_volume>param.V)
  //  m_max_volume=param.V;
  //else
  //  m_max_volume+=param.dt/m_tau_V*m_alpha_V*(m_volume/m_max_volume);
  // else if(0.8>m_volume/m_max_volume)
  //   m_max_volume+=param.dt/m_tau_V*m_alpha_V;
  //m_max_volume+=param.dt/m_tau_V*m_alpha_V*(m_volume/m_max_volume);
  //exponential
  //m_max_volume+=param.dt/m_tau_V*m_alpha_V*(m_max_volume+param.V);
}

void Cell::add_time(){
  m_time+=param.dt;
}


void Cell::copy_param(Cell*& cell){

  m_tau=cell->Get_tau();
  m_D  =cell->Get_D();
  m_tol=cell->Get_tol();
  m_gamma_curv=cell->Get_gamma_curv();

  m_alpha=cell->Get_alpha();
  m_beta =cell->Get_beta();
  m_gamma=cell->Get_gamma();
  m_eta  =cell->Get_eta();

  m_tau_V=param.tau_V+(MT_rand()*2.0f-1.0f)*param.noise_tauV*param.tau_V;

  //m_tau_V  =cell->Get_tau_V();
  m_alpha_V=cell->Get_alpha_V();

  m_imin=cell->Get_imin();
  m_jmin=cell->Get_jmin();
  m_Gx_cell=cell->Get_Gx_cell();
  m_Gy_cell=cell->Get_Gy_cell();

  m_td=param.td+(MT_rand()*2.0-1.0)*param.noise_td*param.td;

}


class Lumen: public PhaseField{

private:
  float m_rs;
  float m_beta, m_xi, m_xi0;

  float m_eta_s, m_gamma_s; //cell-lumen adhesion

  bool m_leak;              //lumen leakage ON/OFF

  float m_alpha_s;          //local volume conservation
  float m_V_s;              //local volume conservation

public:
  Lumen();
  Lumen(const int xmax, const int ymax);

public:
  float Get_rs();
  float Get_beta();
  float Get_xi();
  float Get_xi0();
  void init_lumen();
  void input_lumen(const int init_t);
  void update_lumen();

  //void time_evolution(Interaction *phi,Wall *wall,VInteraction* polfield);
  void time_evolution(Interaction *phi);
  void time_evolution(Interaction *phi,ECM *ecm);
  void time_evolution(Interaction *phi,ECM *ecm,Interaction *pcx);
  void time_evolution(Interaction *phi,Interaction *u_all,ECM *ecm,Interaction *pcx);
  float f(float phi,float w,vector3 nablahu,vector3 cellpol);
  void update_xi1(Interaction *phi);
  void update_xi2(Interaction *phi);
  void update_xi3(Interaction *phi,ECM *ecm);

  float Get_eta_s();
  float Get_gamma_s();
  float f(float psi,float w,vector3 nablahu,vector3 cellpol,float psilaplacian,int i,int j);

}lumen ;

inline Lumen::Lumen(){}
inline Lumen::Lumen(const int xmax, const int ymax):PhaseField(xmax,ymax){}

inline float Lumen::Get_rs(){return m_rs;}
inline float Lumen::Get_beta(){return m_beta;}
inline float Lumen::Get_xi(){return m_xi;}
inline float Lumen::Get_xi0(){return m_xi0;}
inline void Lumen::init_lumen(){
  m_tau=param.tau_s;
  m_D=param.Ds;
  m_tol=param.tol_s;
  m_rs=param.r_s;
  m_beta=param.beta_s;
  m_xi=param.xi0;
  m_xi0=param.xi0;
  init_field();

  m_eta_s=param.eta_s;
  m_gamma_s=param.gamma_s;

  m_leak=false;

  m_gamma_curv=param.gamma_curv_s;

  m_alpha_s=param.alpha_s;
  m_V_s=param.V_s;
}


inline void Lumen::input_lumen(const int init_t){
  m_tau=param.tau_s;
  m_D=param.Ds;
  m_rs=param.r_s;
  m_beta=param.beta_s;
  m_xi=param.xi0;
  m_xi0=param.xi0;

  m_eta_s=param.eta_s;
  m_gamma_s=param.gamma_s;

  //set position
  stringstream ss;
  ss<<init_t;
  string fname;
  fname="STATE/s_"+ss.str()+".dat";
  ifstream fin(fname.c_str()); 
  if(!fin.is_open()){cerr<<"ERROR:Could not open input file "<<fname<<endl;exit(8);}
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      fin>>m_pf[i][j];
  }
  fin.close();
}

inline float Lumen::Get_eta_s(){return m_eta_s;}
inline float Lumen::Get_gamma_s(){return m_gamma_s;}

inline void Lumen::update_lumen(){
  update();
  m_volume=volume_h();
}

class ECM: public PhaseField{

private:
  float m_alpha;
  float m_beta_cu, m_beta_cs;
  float m_xi_c;

  float m_eta_cu, m_eta_cs, m_gamma;

  float m_Vsys;

public:
  ECM();
  ECM(const int xmax, const int ymax);

public:
  float Get_alpha();
  float Get_beta_cu();
  float Get_beta_cs();
  float Get_xi_c();
  float Get_eta_cu();
  float Get_eta_cs();
  float Get_gamma();
  void init_ecm(Interaction* u_all);
  void input_ecm(const int init_t);
  void update_ecm();

  void time_evolution(Interaction *phi,Lumen *lumen);
  void time_evolution(Interaction *phi);             

  void minuscell(Cell* cell);

  bool check_leakage(Lumen* lumen);

}ecm ;

inline ECM::ECM(){}
inline ECM::ECM(const int xmax, const int ymax):PhaseField(xmax,ymax){}

inline float ECM::Get_alpha(){return m_alpha;}
inline float ECM::Get_beta_cu(){return m_beta_cu;}
inline float ECM::Get_beta_cs(){return m_beta_cs;}
inline float ECM::Get_xi_c(){return m_xi_c;}
inline float ECM::Get_eta_cu(){return m_eta_cu;}
inline float ECM::Get_eta_cs(){return m_eta_cs;}
inline float ECM::Get_gamma(){return m_gamma;}

inline void ECM::update_ecm(){
  update();
  m_volume=volume_h();
}

bool ECM::check_leakage(Lumen* lumen){

  bool leakage=false;
  float vsum=0.0f;

  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      vsum+=m_pf[i][j]*lumen->Get_value(i,j)*param.dx*param.dy;
    }
  }
  if(vsum>param.dx*3.0f) leakage=true;

  return leakage;
}

class Wall: public Field{
private:
  float m_beta, m_eta;

public:
  Wall();
  Wall(const int xmax, const int ymax);

  float Get_beta();
  float Get_eta();
  void Set_param();

  void Set_wall_input();
  void Set_wall_circle();
};

inline Wall::Wall(){}
inline Wall::Wall(const int xmax, const int ymax):Field(xmax,ymax){}

inline float Wall::Get_beta(){return m_beta;}
inline float Wall::Get_eta(){return m_eta;}
inline void Wall::Set_param(){
  m_beta=param.beta_w;
  m_eta =param.eta_w;
}

inline void Wall::Set_wall_input(){
  ifstream ifwall("wall.dat");
  if(!ifwall.is_open()){cerr<<"ERROR:Could not open input file wall data"<<endl;exit(8);}
  for(int j=0;j<m_ymax;j++){
    for(int i=0;i<m_xmax;i++){
      ifwall>>m_f[i][j];
    }
  }
}

inline void Wall::Set_wall_circle(){
  float center_x=m_xmax*0.5f*param.dx;
  float center_y=m_ymax*0.5f*param.dy;
#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      if((param.dx*(i+0.5f)-center_x)*(param.dx*(i+0.5f)-center_x)+(param.dy*(j+0.5f)-center_y)*(param.dy*(j+0.5f)-center_y)<4.0*4.0)
	m_f[i][j]=0.0f;
      else
	m_f[i][j]=1.0f;
    }
  }
}

class Interaction: public Field{

private:
  bool m_e_eta_null;

  float m_volume;

public:
  Interaction();
  Interaction(const int xmax, const int ymax);

public:

  void Set_phi(Cell** cells,const int num_cells);
  void Set_phi2(Cell** cells,const int num_cells,const int num);
  void pluscell(Cell*& cell);
  void minuscell(Cell*& cell);

  bool Get_e_eta_null();
  void Set_e_eta(Cell** cells, int num_cells);
  void reset_e_eta();

  void Set_u_all(Cell** cells,const int num_cells);

  float Get_volume();

  void Set_e_adhe(Cell** cells, int num_cells);
  float cal_emax_x();
  float cal_emax_y();

  void Set_k();
  void init_pf(ECM* ecm,Lumen* lumen);

}interaction ;

inline Interaction::Interaction(){}
inline Interaction::Interaction(const int xmax, const int ymax):Field(xmax,ymax){}

void Interaction::Set_phi(Cell**cells,int num_cells){
  for(int n=0;n<num_cells;n++){
    int imin,jmin;
    imin=cells[n]->Get_imin(); jmin=cells[n]->Get_jmin();
    //cout<<"in Set_phi imin,jmin= "<<imin<<","<<jmin<<endl;
    for(int i=0;i<param.xmax_cell;i++){
      for(int j=0;j<param.ymax_cell;j++){
	m_f[i+imin][j+jmin]+=h(cells[n]->Get_value(i,j));
      }
    }
  }
  m_volume=volume();
}

void Interaction::Set_phi2(Cell**cells,int num_cells,int num0){
  for(int n=0;n<num_cells;n++){
    if(n!=num0){
      int imin,jmin;
      imin=cells[n]->Get_imin(); jmin=cells[n]->Get_jmin();
      for(int i=0;i<param.xmax_cell;i++){
	for(int j=0;j<param.ymax_cell;j++){
	  m_f[i+imin][j+jmin]+=h(cells[n]->Get_value(i,j));
	}
      }
    }
  }
  m_volume=volume();
}


void Interaction::pluscell(Cell*& cell){
  int imin,jmin;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      m_f[i+imin][j+jmin]+=h(cell->Get_value(i,j));
    }
  }
}

void Interaction::minuscell(Cell*& cell){
  int imin,jmin;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      m_f[i+imin][j+jmin]-=h(cell->Get_value(i,j));
    }
  }
}

inline bool Interaction::Get_e_eta_null(){return m_e_eta_null;}

inline void Interaction::Set_e_eta(Cell **cells,int num_cells){

  int imin,jmin,imin1,jmin1,imax,jmax,imax1,jmax1;

  Field* e_eta_m = new Field(param.xmax,param.ymax);
  Field* e_eta_n = new Field(param.xmax,param.ymax);
  e_eta_m->init_field();
  e_eta_n->init_field();

  for(int m=0;m<num_cells;m++){
    imin=cells[m]->Get_imin(); jmin=cells[m]->Get_jmin();
    imax=imin+param.xmax_cell-1;
    jmax=jmin+param.ymax_cell-1;

    e_eta_m->init_field();
    e_eta_n->init_field();

#pragma omp parallel for num_threads(1)
    for(int i=0;i<param.xmax_cell;i++){
      for(int j=0;j<param.ymax_cell;j++){
	e_eta_m->Set_value(i+imin,j+jmin,cells[m]->Get_value(i,j));
      }
    }
    for(int n=0;n<num_cells;n++){
      if(m!=n){
	imin1=cells[n]->Get_imin();
	jmin1=cells[n]->Get_jmin();

	imax1=imin1+param.xmax_cell-1;
	jmax1=jmin1+param.ymax_cell-1;

	if(max(imin,imin1) <= min(imax,imax1) && max(jmin,jmin1) <= min(jmax,jmax1)){
	  for(int i=0;i<param.xmax_cell;i++){
	    for(int j=0;j<param.ymax_cell;j++){
	      e_eta_n->Set_value(i+imin1,j+jmin1,e_eta_n->Get_value(i+imin1,j+jmin1)+cells[n]->Get_value(i,j));
	    }
	  }
	}
      }
    }
    for(int i=1;i<param.xmax-1;i++){
      for(int j=1;j<param.ymax-1;j++){
	m_f[i][j]+=(e_eta_m->nabla_h(i,j)%e_eta_n->nabla_h(i,j))*param.eta/6.0;
      }
    }
  }

  delete e_eta_m;
  delete e_eta_n;

  m_e_eta_null=false;
}

inline void Interaction::Set_e_adhe(Cell **cells,int num_cells){
  int imin,jmin;
  
  Field* e_eta_0 = new Field(m_xmax,m_ymax);
  Field* e_eta_1 = new Field(m_xmax,m_ymax);
  e_eta_0->init_field(); e_eta_1->init_field();
  
  imin=cells[0]->Get_imin(); jmin=cells[0]->Get_jmin();
#pragma omp parallel for num_threads(1)
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      e_eta_0->Set_value(i+imin,j+jmin,cells[0]->Get_value(i,j));
    }
  }
  imin=cells[1]->Get_imin(); jmin=cells[1]->Get_jmin();
#pragma omp parallel for num_threads(1)
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      e_eta_1->Set_value(i+imin,j+jmin,cells[1]->Get_value(i,j));
    }
  }
#pragma omp parallel for num_threads(1)
  for(int i=1;i<param.xmax-1;i++){
    for(int j=1;j<param.ymax-1;j++){
      m_f[i][j]=(e_eta_0->nabla_h(i,j)%e_eta_1->nabla_h(i,j))*param.eta/6.0;
    }
  }

  delete e_eta_0;
  delete e_eta_1;

}

float Interaction::cal_emax_x(){
  float emax_x,sum;
  emax_x=0.0f;
  for(int i=0;i<m_xmax;i++){
    sum=0.0f;
    for(int j=0;j<m_ymax;j++){sum+=m_f[i][j]*param.dx*param.dy;}
    if(sum<emax_x)
      emax_x=sum;
  }
  return emax_x;
}
 
float Interaction::cal_emax_y(){
  float emax_y, sum;
  emax_y=0.0f;
  for(int j=0;j<m_ymax;j++){
    sum=0.0f;
    for(int i=0;i<m_xmax;i++){sum+=m_f[i][j]*param.dx*param.dy;}
    if(sum<emax_y)
      emax_y=sum;
  }
  return emax_y;
}

void Interaction::reset_e_eta(){
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      m_f[i][j]=0.0f;
    }
  }
  m_e_eta_null=true;
}

void Interaction::Set_u_all(Cell**cells,int num_cells){
  for(int n=0;n<num_cells;n++){
  int imin,jmin;
  imin=cells[n]->Get_imin(); jmin=cells[n]->Get_jmin();
#pragma omp parallel for num_threads(1)
    for(int i=0;i<param.xmax_cell;i++){
      for(int j=0;j<param.ymax_cell;j++){
	m_f[i+imin][j+jmin]+=cells[n]->Get_value(i,j);
      }
    }
  }
}



inline float Interaction::Get_volume(){return m_volume;}


void Interaction::Set_k(){
  for(int i=0;i<param.xmax;i++){
    for(int j=0;j<param.ymax;j++) m_f[i][j]=0.0f;//sinf(i/10.0f*M_PI)*sinf(j/10.0f*M_PI);
  }
   m_f[0][0]=1.0f+4.0f*param.D_pol*param.tau_pol/(param.dx*param.dx);
   m_f[0][1]=-param.D_pol*param.tau_pol/(param.dx*param.dx);
   m_f[1][0]=-param.D_pol*param.tau_pol/(param.dx*param.dx);
   m_f[0][param.ymax-1]=-param.D_pol*param.tau_pol/(param.dx*param.dx);
   m_f[param.xmax-1][0]=-param.D_pol*param.tau_pol/(param.dx*param.dx);
}

void Interaction::init_pf(ECM* ecm,Lumen* lumen){
  for(int i=0;i<param.xmax;i++){
    for(int j=0;j<param.ymax;j++){
      m_f[i][j]=(ecm->Get_value(i,j)-lumen->Get_value(i,j)+1.0f)*0.5f;
    }
  }
}




class VInteraction: public VField{

private:
  bool m_null;

public:
  VInteraction();
  VInteraction(const int xmax, const int ymax);

public:
  int Get_xmax(), Get_ymax();
  vector3 Get_vector(const int i, const int j);
  void init_field();

  void Set_polfield(Cell** cells,const int num_cells);
  void plusvec(Cell*& cell);
  void minusvec(Cell*& cell);

}vinteraction ;

inline VInteraction::VInteraction(){}
inline VInteraction::VInteraction(const int xmax, const int ymax):VField(xmax,ymax){}

inline int VInteraction::Get_xmax(){return m_xmax;}
inline int VInteraction::Get_ymax(){return m_ymax;}
inline vector3 VInteraction::Get_vector(const int i, const int j){return m_f[i][j];}

void VInteraction::init_field(){
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
	m_f[i][j].set(0.0f,0.0f,0.0f);
    }
  }
}

void VInteraction::Set_polfield(Cell**cells,int num_cells){
  int imin,jmin;
  for(int n=0;n<num_cells;n++){
    imin=cells[n]->Get_imin(); jmin=cells[n]->Get_jmin();
    for(int i=0;i<param.xmax_cell;i++){
      for(int j=0;j<param.ymax_cell;j++){
	if(cells[n]->Get_value(i,j)>0.1)
	  m_f[i+imin][j+jmin]+=cells[n]->Get_pol();
      }
    }
  }
}

void VInteraction::plusvec(Cell*& cell){
  int imin,jmin;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      m_f[i+imin][j+jmin]+=cell->Get_pol();
    }
  }
}

void VInteraction::minusvec(Cell*& cell){
  int imin,jmin;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();
  for(int i=0;i<param.xmax_cell;i++){
    for(int j=0;j<param.ymax_cell;j++){
      m_f[i+imin][j+jmin]-=cell->Get_pol();
    }
  }
}

//--------------------------------------------------------------------------//
// Cross-reference classes
//--------------------------------------------------------------------------//


inline void Cell::seed_lumen(float* Pc,Lumen* lumen){

  float r=0.0f;
  float seed=0.0f;
#pragma omp parallel for num_threads(1)
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      r=sqrtf((param.dx*(i+m_imin+0.5f)-Pc[0])*(param.dx*(i+m_imin+0.5f)-Pc[0])
	     +(param.dy*(j+m_jmin+0.5f)-Pc[1])*(param.dy*(j+m_jmin+0.5f)-Pc[1]));
      //seed=(1.0f-tanhf((r-lumen->Get_rs())/(2.0f*sqrtf(2.0f*m_D))))*0.5f;
      seed=(1.0f-tanhf((r-lumen->Get_rs())/(sqrtf(2.0f*m_D))))*0.5f;
      if(seed>0.00001){
	lumen->Set_value(i+m_imin,j+m_jmin,lumen->Get_value(i+m_imin,j+m_jmin)+seed);
	if(lumen->Get_value(i+m_imin,j+m_jmin)>1.0f)
	  lumen->Set_value(i+m_imin,j+m_jmin,1.0f);
	m_pf[i][j]-=seed;
	if(m_pf[i][j]<0.0f)
	  m_pf[i][j]=0.0f;
      }
    }
  }
}

void Cell::time_evolution(Interaction *phi){
  if(m_tol>0.00001f){
#pragma omp parallel for num_threads(1)
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*f(i,j,phi))*param.dt/m_tau
                         	  +2.0f*m_gamma_curv*sqrtf(m_D)*laplacian8(i,j)*param.dt/m_tau
	;
      }
    }
  }
  else{
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)
				    +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(i,j,phi))
				    )*param.dt/m_tau;
      }
    }
  }
}




void Cell::time_evolution(Lumen *lumen,Interaction *phi){
//   if(m_tol>0.000001){
// #pragma omp parallel for num_threads(1)
//     for(int i=1;i<m_xmax-1;i++){
//       for(int j=1;j<m_ymax-1;j++){
// 	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*f(i,j,phi,lumen))*param.dt/m_tau
//                                               +2.0f*m_gamma_curv*sqrtf(m_D)*laplacian8(i,j)*param.dt/m_tau
// 	;
//       }
//     }
//   }
//   else{
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++){
      m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)
				  +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(i,j,phi,lumen))
				  )*param.dt/m_tau;
    }
  }
    //  }
}

void Cell::time_evolution(Lumen *lumen,Interaction *phi,ECM *ecm){
//   if(m_tol>0.000001){
// #pragma omp parallel for num_threads(1)
//     for(int i=1;i<m_xmax-1;i++){
//       for(int j=1;j<m_ymax-1;j++){
// 	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*f(i,j,phi,lumen,ecm))*param.dt/m_tau
//                          	  +2.0f*m_gamma_curv*sqrtf(m_D)*laplacian8(i,j)*param.dt/m_tau
// 	;
//       }
//     }
//   }
//   else{
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++){
      m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)
				  +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(i,j,phi,lumen,ecm))
				  )*param.dt/m_tau;
    }
  }
  //}
}

void Cell::time_evolution(Interaction *phi,ECM *ecm){
  if(m_tol>0.000001){
#pragma omp parallel for num_threads(1)
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*f(i,j,phi,ecm))*param.dt/m_tau
                         	  +2.0f*m_gamma_curv*sqrtf(m_D)*laplacian8(i,j)*param.dt/m_tau 
	;
      }
    }
  }
  else{
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)
				    +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(i,j,phi,ecm))
				    )*param.dt/m_tau;
      }
    }
  }
}



float Cell::f(int& i,int& j,Interaction* phi){
  return
    m_alpha*(m_max_volume-m_volume)                          // the constraint of volume
    -m_beta*(phi->Get_value(i+m_imin,j+m_jmin))              // the excluded volume of cells
    +m_eta*phi->laplacian8(i+m_imin,j+m_jmin)                // the adhesion with other cells
    +m_gamma*laplacian8_h(i,j)                               // the adhesion with other cells
    ;
}



float Cell::f(int& i,int& j,Interaction* phi,Lumen* lumen){
  return
    m_alpha*(m_max_volume-m_volume)                          // volume conservation
    -m_beta*(phi->Get_value(i+m_imin,j+m_jmin))              // excluded volume of cells
    -lumen->Get_beta()*h(lumen->Get_value(i+m_imin,j+m_jmin))// excluded volume of lumens

    +m_eta*phi->laplacian8(i+m_imin,j+m_jmin)                // the adhesion with other cells
    +m_gamma*laplacian8_h(i,j)                               // the adhesion with other cells

    //+lumen->Get_eta_s()*lumen->laplacian8_h(i+m_imin,j+m_jmin)             // the adhesion with lumens
    ;
}

float Cell::f(int& i,int& j,Interaction* phi,Lumen* lumen,ECM *ecm){
  return
    m_alpha*(m_max_volume-m_volume)                                // volume conservation

    //excluded volume
    -m_beta*(phi->Get_value(i+m_imin,j+m_jmin))                    //cells
    -lumen->Get_beta()*h(lumen->Get_value(i+m_imin,j+m_jmin))      //lumens
    -ecm->Get_beta_cu()*h(ecm->Get_value(i+m_imin,j+m_jmin))       //medium

    //adhesion
    +m_eta*phi->laplacian8(i+m_imin,j+m_jmin)                      // cells
    +m_gamma*laplacian8_h(i,j)                                     // 
    +lumen->Get_eta_s()*lumen->laplacian8_h(i+m_imin,j+m_jmin)     // cell-lumen
    +ecm->Get_eta_cu()*ecm->laplacian8_h(i+m_imin,j+m_jmin)        // cell-medium
    ;
}


float Cell::f(int& i,int& j,Interaction* phi,ECM *ecm){
  return
    m_alpha*(m_max_volume-m_volume)    // the constraint of volume

    -m_beta*(phi->Get_value(i+m_imin,j+m_jmin))            // the excluded volume of cells

    +m_eta*phi->laplacian8(i+m_imin,j+m_jmin)                // the adhesion with other cells
    +m_gamma*laplacian8_h(i,j)                               //

    +ecm->Get_eta_cu()*ecm->laplacian8_h(i+m_imin,j+m_jmin)             // the adhesion with ecm
    -ecm->Get_beta_cu()*h(ecm->Get_value(i+m_imin,j+m_jmin))
    ;
}



// void Cell::time_evolution(Lumen *lumen,Interaction *phi,Wall *wall,VInteraction *polfield){

// #pragma omp parallel for num_threads(1)
//   for(int i=1;i<m_xmax-1;i++){
//     for(int j=1;j<m_ymax-1;j++){
//       m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)
// 				  +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(i,j,phi,lumen,wall,polfield))
// 				  )*param.dt/m_tau;
//     }
//   }
// }
// float Cell::f(int& i,int& j,Interaction* phi,Lumen* lumen,Wall* wall,VInteraction* polfield){
//   return
//     //f_u M.Nonomura(2012)&M.Akiyama(2018)
//     //param.alpha*(param.V-cell->Get_volume())                 // the constraint of volume
//     // add time evolution of max cell volume 2019----------------------- start
//     m_alpha*(m_max_volume-m_volume)    // the constraint of volume
//     // add time evolution of max cell volume 2019----------------------- end

//     -m_beta*(phi->Get_value(i+m_imin,j+m_jmin))              // the excluded volume of cells
//     //-param.beta*(psi*kdelta(0,0))                            // the excluded volume of cells
//     -lumen->Get_beta()*h(lumen->Get_value(i+m_imin,j+m_jmin))  // the excluded volume of lumens
//     -param.beta_w*h(wall->Get_value(i+m_imin,j+m_jmin))      // the excluded volume of the wall

//     //g_int M.Nonomura(2012)
//     +m_eta*phi->laplacian8(i+m_imin,j+m_jmin)                // the adhesion with other cells
//     +m_gamma*laplacian8_h(i,j)                               // the adhesion with other cells
    
//     //g_s M.Nonomura(2012)
//     +param.eta_w*wall->laplacian8_h(i+m_imin,j+m_jmin)       // the adhesion with the wall

//     // add adhesion between cells and lumens 20200512----------------------- start
//     +lumen->Get_eta_s()*lumen->laplacian8_h(i+m_imin,j+m_jmin)             // the adhesion with lumens
//     // add adhesion between cells and lumens 20200512----------------------- end

//     // add cell polarity 20191218----------------------- start
//     //+param.kappa_n*(cell->nabla_h(i,j)%field->nabla_h(i+imin,j+jmin))*(cell->Get_pol()+polfield->Get_vector(i+imin,j+jmin))
//     //+param.lambda_n*(lumen->nabla_h(i+imin,j+jmin)%cell->Get_pol())
//     // add cell polarity 20191218----------------------- end
//     ;
// }



// void Lumen::time_evolution(Interaction *phi,Wall *wall,VInteraction* polfield){
// #pragma omp parallel for num_threads(1)
//   for(int i=1;i<m_xmax-1;i++){
//     for(int j=1;j<m_ymax-1;j++){
//       m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)+
// 				  m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f+f(phi->Get_value(i,j),
// 									      wall->Get_value(i,j),
// 									      phi->nabla_h(i,j),
// 									      polfield->Get_vector(i,j),
// 									      // add adhesion between cells and lumens 20200512----------------------- start
// 									      phi->laplacian8(i,j),
// 									      i,j
// 									      )
// 									      // add adhesion between cells and lumens 20200512----------------------- end
// 							     )
// 				  )*param.dt/m_tau;
//     }
//   }
// }

//float Lumen::f(float psi,float w,vector3 nablahu,vector3 cellpol){
// add adhesion between cells and lumens 20200512----------------------- start
float Lumen::f(float psi,float w,vector3 nablahu,vector3 cellpol,float psilaplacian,int i,int j){
// add adhesion between cells and lumens 20200512----------------------- end
  return
    m_xi                               // lumen pressure
    -m_beta*psi                        // the excluded volume of cells
    -param.beta_w*h(w)                 // the excluded volume of wall
    // add adhesion between cells and lumens 20200512----------------------- start
    +m_eta_s*psilaplacian              // the adhesion with cells
    +m_gamma_s*laplacian8_h(i,j)       // the adhesion with cells
    // add adhesion between cells and lumens 20200512----------------------- end


    //add cell_polarity 20191216-------------------------start 
    //+param.lambda_n*(nablahu%cellpol)
    //add cell_polarity 20191216-------------------------end
    ;
}

void Lumen::time_evolution(Interaction *phi){
//   if(m_tol>0.000001){
// #pragma omp parallel for num_threads(1)
//     for(int i=1;i<m_xmax-1;i++){
//       for(int j=1;j<m_ymax-1;j++){
// 	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*(+m_xi
// 							       -m_beta*phi->Get_value(i,j)
// 							       // cell-lumen adhesion 20200512
// 							       +m_eta_s*phi->laplacian8(i,j)
// 							       )
// 				    )*param.dt/m_tau
// 	  // //20200808 curvature effects------------------------------------------------ start
//           //                    	  +m_gamma_curv*(m_D*laplacian8(i,j)
// 	  // 					 +m_pf[i][j]*(1-m_pf[i][j])*(m_pf[i][j]-0.5f)
// 	  // 					 )*param.dt/m_tau
// 	  // //20200808 curvature effects------------------------------------------------ end
// 	  ;
//       }
//     }
//   }
//   else{
// #pragma omp parallel for num_threads(1)
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++){
      m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)+
				  m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f
								+m_xi
								-m_beta*phi->Get_value(i,j)
								// cell-lumen adhesion 20200512
								// +m_eta_s*phi->laplacian8(i,j)
								// +m_gamma_s*laplacian8_h(i,j)
								)
				  )*param.dt/m_tau;
    }
  }
    //  }
}

void Lumen::time_evolution(Interaction *phi,ECM *ecm){
//   if(m_tol>0.000001){
// #pragma omp parallel for num_threads(1)
//     for(int i=1;i<m_xmax-1;i++){
//       for(int j=1;j<m_ymax-1;j++){
// 	m_pf_temp[i][j]=m_pf[i][j]
// 	               +(m_pf[i][j]*(1.0f-m_pf[i][j])*(+m_xi
// 						    -m_beta*phi->Get_value(i,j)
// 						    // cell-lumen adhesion 20200512
// 						    +m_eta_s*phi->laplacian8(i,j)
// 						    // ecm 20200527
// 						    -ecm->Get_beta_cs()*h(ecm->Get_value(i,j))
// 						    +ecm->Get_eta_cs()*ecm->laplacian8_h(i,j) 
// 						    )
// 			 )*param.dt/m_tau
// 	  ;
// 	// //20200808 curvature effects------------------------------------------------ start
// 	//              +m_gamma_curv*(m_D*laplacian8(i,j)
// 	// 			      +m_pf[i][j]*(1-m_pf[i][j])*(m_pf[i][j]-0.5f)
// 	// 			      )*param.dt/m_tau
// 	// //20200808 curvature effects------------------------------------------------ end
//       }
//     }
//   }
//   else{
// #pragma omp parallel for num_threads(1)
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++){
      m_pf_temp[i][j]=m_pf[i][j]
	+(m_D*laplacian8(i,j)+
	  m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f
					+m_xi
					-m_beta*phi->Get_value(i,j)
					// cell-lumen adhesion 20200512
					+m_eta_s*phi->laplacian8(i,j)
					+m_gamma_s*laplacian8_h(i,j)
					// ecm 20200527
					-ecm->Get_beta_cs()*h(ecm->Get_value(i,j))
					+ecm->Get_eta_cs()*ecm->laplacian8_h(i,j) 
					)
	  )*param.dt/m_tau;
    }
  }
  //}
}



void Lumen::update_xi1(Interaction *u_all){
  float overlap=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      overlap+=m_pf[i][j]*(1.0f-m_pf[i][j])*(1.0f-u_all->Get_value(i,j))*param.dx*param.dy;
  }
  if(overlap>0)
    m_xi=m_xi0-param.c_s*overlap/surface();
  else
    m_xi=m_xi0;
  //cout<<m_xi0<<" "<<overlap<<" "<<surface()<<" "<<m_volume<<endl;
}

//add feedback control of lumen pressure 2 20200512 ---------------------start
void Lumen::update_xi2(Interaction *u_all){
  float overlap=0.0f;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      overlap+=m_pf[i][j]*(1.0f-m_pf[i][j])*(u_all->Get_value(i,j))*param.dx*param.dy;
  }

  if(overlap>0){
    //m_xi=m_xi0*overlap/m_volume;
    m_xi=m_xi0*overlap/volume();
  }
  else
    m_xi=0;

  //cout<<m_xi0<<" "<<overlap<<" "<<surface()<<" "<<m_volume<<endl;
}
//add feedback control of lumen pressure 2 20200512 ---------------------end

//add feedback control of lumen pressure 3 20200729 ---------------------start
void Lumen::update_xi3(Interaction *u_all,ECM *ecm){
  float a=5.0;
  float overlap=0.0f;
  float v_leak=0.0f;
  float V_leak=5*param.dx*param.dy;
  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++){
      overlap+=m_pf[i][j]*(1.0f-m_pf[i][j])*(u_all->Get_value(i,j))*param.dx*param.dy;
      v_leak+=m_pf[i][j]*(1.0f-m_pf[i][j])*ecm->Get_value(i,j)*param.dx*param.dy;
    }
  }

  //check leakage
  if(v_leak>V_leak && m_leak==false)
    m_leak=true;
  else if(v_leak<param.dx*param.dy && m_leak==true)
    m_leak=false;

  //calculation of lumen pressure
  if(overlap>0){
    if(m_leak==true)
      m_xi=m_xi0*overlap/m_volume*exp(-a*(v_leak/V_leak));
    else
      m_xi=m_xi0*overlap/m_volume;
  }
  else
    m_xi=0;

  //cout<<m_leak<<" "<<overlap<<" "<<v_leak<<" "<<m_volume<<endl;
}
//add feedback control of lumen pressure 3 20200729 ---------------------end


inline void ECM::init_ecm(Interaction* u_all){
  m_tau=param.tau_c;
  m_D=param.Dc;
  m_tol=param.tol_c;
  m_alpha=param.alpha_c;
  m_beta_cu=param.beta_cu;
  m_beta_cs=param.beta_cs;
  m_xi_c=param.xi_c;
  m_eta_cu=param.eta_cu;
  m_eta_cs=param.eta_cs;
  m_gamma=param.gamma_c;

  m_Vsys=param.Vsys;

  m_gamma_curv=param.gamma_curv_c; //20200809

  for(int i=0;i<m_xmax;i++){
    for(int j=0;j<m_ymax;j++)
      m_pf[i][j]=1.0f-u_all->Get_value(i,j);
  }

  m_volume=volume();
}

void ECM::time_evolution(Interaction *phi,Lumen *lumen){

  //float V=m_Vsys-phi->Get_volume()-lumen->Get_volume();
  //cout<<V<<" "<<m_volume;

//   if(m_tol>0.00001f){
// #pragma omp parallel for num_threads(1)
//     for(int i=1;i<m_xmax-1;i++){
//       for(int j=1;j<m_ymax-1;j++){
// 	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*(
// 							       //+m_alpha*(V-m_volume)
// 							       -m_beta_cu*phi->Get_value(i,j)
// 							       -m_beta_cs*h(lumen->Get_value(i,j))
// 							       +m_eta_cu*phi->laplacian8(i,j)
// 							       +m_eta_cs*lumen->laplacian8_h(i,j)
// 							       +m_xi_c
// 							       )
// 				    )*param.dt/m_tau
// 	  // //20200808 curvature effects------------------------------------------------ start
//           //              +m_gamma_curv*(m_D*laplacian8(i,j)
// 	  // 			      +m_pf[i][j]*(1-m_pf[i][j])*(m_pf[i][j]-0.5f)
// 	  // 			      )*param.dt/m_tau
// 	  // //20200808 curvature effects------------------------------------------------ end
// 	  ;
//       }
//     }
//   }
//   else{
// #pragma omp parallel for num_threads(1)
  for(int i=1;i<m_xmax-1;i++){
    for(int j=1;j<m_ymax-1;j++){
      m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)+
				  m_pf[i][j]*(1.0f-m_pf[i][j])*(
								m_pf[i][j]-0.5f
								+m_alpha*(m_Vsys-m_volume)
								-m_beta_cu*phi->Get_value(i,j)
								-m_beta_cs*h(lumen->Get_value(i,j))
								//+m_eta_cu*phi->laplacian8(i,j)
								//+m_eta_cs*lumen->laplacian8_h(i,j)
								//+m_gamma*laplacian8_h(i,j)
								//+m_xi_c
								)
				  )*param.dt/m_tau;
    }
  }
  //}
}



void ECM::time_evolution(Interaction *phi){

  //float V=m_Vsys-phi->Get_volume()-lumen->Get_volume();
  //cout<<V<<" "<<m_volume;

  if(m_tol>0.000001){
#pragma omp parallel for num_threads(1)
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_pf[i][j]*(1.0f-m_pf[i][j])*(
							       //+m_alpha*(V-m_volume)
							       -m_beta_cu*phi->Get_value(i,j)
							       +m_eta_cu*phi->laplacian8(i,j)
							       +m_xi_c
							       )
				    )*param.dt/m_tau
                       +m_gamma_curv*(m_D*laplacian8(i,j)
				      +m_pf[i][j]*(1.0f-m_pf[i][j])*(m_pf[i][j]-0.5f)
				      )*param.dt/m_tau
	  ;
      }
    }
  }
  else{
#pragma omp parallel for num_threads(1)
    for(int i=1;i<m_xmax-1;i++){
      for(int j=1;j<m_ymax-1;j++){
	m_pf_temp[i][j]=m_pf[i][j]+(m_D*laplacian8(i,j)+
				    m_pf[i][j]*(1.0f-m_pf[i][j])*(
							       m_pf[i][j]-0.5f
							       //+m_alpha*(V-m_volume)
							       -m_beta_cu*phi->Get_value(i,j)
							       +m_eta_cu*phi->laplacian8(i,j)
							       +m_gamma*laplacian8_h(i,j)
							       +m_xi_c
							       )
				    )*param.dt/m_tau;
      }
    }
  }
}


//20200624 remove medium after cell division
void ECM::minuscell(Cell* cell){
  int imin,jmin,xmax_cell,ymax_cell;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();
  xmax_cell=cell->Get_xmax(); ymax_cell=cell->Get_ymax();
#pragma omp parallel for num_threads(1)
  for(int i=0;i<xmax_cell;i++){
    for(int j=0;j<ymax_cell;j++){
      m_pf[i+imin][j+jmin]-=cell->Get_value(i,j);
    }
  }
}



class Calculation: public Commonfunction{
public:

  Calculation();
  //old functions
  void   time_evolution_init(Cell *&cells,Interaction *phi);
  float f_bound(float u, Cell *&cells);
  float f(Cell *&cells,float u);
  float gint(Cell*& cell,Interaction* phi,int& i,int& j);
  float chi(float r1_x,float r1_y,float r_x,float r_y,float pm,float theta);

  //CFG
  void time_evolution_poles(Cell *&cells,Interaction *e_eta,float (&r1)[2],float (&r2)[2],int num_cells);
  float chi(float (&r)[2],float (&r1)[2],float (&r2)[2]);
  float g(float (&r)[2],float (&r1)[2],float (&r2)[2]);
  float H(float u);
  float G(float e_eta);
  float distance(float (&r1)[2],float (&r2)[2]);
  void division(Cell *&mcell, Cell *&dcell_2,Lumen* lumen,Interaction* e_eta,int num_cells);

  void time_evolution_cellpolarity(Cell *&cell,VInteraction* polfield,Interaction* field,Lumen* lumen);
} cal ;


inline Calculation::Calculation(){}

void Calculation::time_evolution_init(Cell *&cell, Interaction *field){

  int imin,jmin;
  imin=cell->Get_imin(); jmin=cell->Get_jmin();

  for(int i=1;i<cell->Get_xmax()-1;i++){
    for(int j=1;j<cell->Get_ymax()-1;j++){
      cell->Set_value_temp(i,j,
			    cell->Get_value(i,j)
			    +
			    (param.D*cell->laplacian8(i,j)
			     +f_bound(cell->Get_value(i,j),cell)
			     +gint(cell,field,i,j)
			     )*param.dt/param.tau
			    );
    }
  }
}

inline float Calculation::f_bound(float u,Cell *&cell){
  return u*(1.0f-u)*(u-0.5f+f(cell,u));
}

float Calculation::f(Cell *&cell,float u){
  return
    param.alpha*(param.V-cell->Get_volume());
}


float Calculation::gint(Cell*& cell,Interaction* field,int& i,int& j){

  int imin=cell->Get_imin();
  int jmin=cell->Get_jmin();
  float u=cell->Get_value(i,j);

  return u*(1.0f-u)*
    (
     param.eta*field->laplacian8(i+imin,j+jmin)
     +param.gamma*cell->laplacian8_h(i,j)
     )
    ;
}

//random division
float Calculation::chi(float r1_x,float r1_y,float r_x,float r_y,float pm,float theta){
  float d_rm1=(r_x*0.5f-r1_x)*cosf(theta)+(r_y*0.5f-r1_y)*sinf(theta);
  return 0.5f*(1+pm*tanhf(d_rm1/param.ep_d));
}



void Calculation::time_evolution_poles(Cell*& cell,Interaction* field, float (&r1)[2], float (&r2)[2], int num_cells){

  int DIM = 2;

  int imin=cell->Get_imin();
  int jmin=cell->Get_jmin();

  float r[2]={}; float dr[2]={};
  float F1[2]={}; float F2[2]={};
  dr[0]=param.dx;  dr[1]=param.dy;

  float r1b[2]={}; float r2b[2]={};
  for(int d=0;d<DIM;d++){
    r1b[d]=0.0f; r2b[d]=0.0f;
  }

  float theta;
  float noise1[2]={}; float noise2[2]={};

  int CFGcount=0;
  int CFGcount2=0;
  while(CFGcount<5000 && CFGcount2<50000){
  
    for(int d=0;d<DIM;d++){
      r1b[d]=r1[d]; r2b[d]=r2[d];
      F1[d]=0; F2[d]=0;
    }
    for(int i=0;i<cell->Get_xmax();i++){
      r[0]=(float)(cell->Get_Gx_cell()-cell->Get_xmax()*0.5f*dr[0]+(i+0.5f)*dr[0]);
      for(int j=0;j<cell->Get_ymax();j++){
 	r[1]=(float)(cell->Get_Gy_cell()-cell->Get_ymax()*0.5f*dr[1]+(j+0.5f)*dr[1]);
	if(distance(r1,r)<distance(r2,r)){
 	  for(int d=0;d<DIM;d++){
 	    F1[d]+=G(field->Get_value(i+imin,j+jmin))*H(cell->Get_value(i,j))*(r[d]-r1[d])*dr[d]*dr[d];
 	  }
 	}
 	if(distance(r1,r)>distance(r2,r)){
 	  for(int d=0;d<DIM;d++){
 	    F2[d]+=G(field->Get_value(i+imin,j+jmin))*H(cell->Get_value(i,j))*(r[d]-r2[d])*dr[d]*dr[d];
 	  }
 	}
      }
    }

    // noise
    if(num_cells!=1){
      theta=MT_rand()*M_PI*2.0f;
      // noise1[0]=-cosf(theta)*0.001f; noise1[1]=-sinf(theta)*0.001f;
      // noise2[0]=cosf(theta+M_PI)*0.001f; noise2[1]=sinf(theta+M_PI)*0.001f;
      noise1[0]=-cosf(theta)*param.dx; noise1[1]=-sinf(theta)*param.dy;
      noise2[0]=cosf(theta+M_PI)*param.dx; noise2[1]=sinf(theta+M_PI)*param.dy;
    }

    for(int d=0;d<DIM;d++){
      if(CFGcount>4000 || CFGcount>40000){
 	r1[d]+=(F1[d]-param.sigma*(distance(r1,r2)-param.ls)*((r1[d]-r2[d])/distance(r1,r2)))/param.mu*param.dts;
 	r2[d]+=(F2[d]-param.sigma*(distance(r1,r2)-param.ls)*((r2[d]-r1[d])/distance(r1,r2)))/param.mu*param.dts;
      }
      else{
 	r1[d]+=(F1[d]-param.sigma*(distance(r1,r2)-param.ls)*((r1[d]-r2[d])/distance(r1,r2)))/param.mu*param.dts+noise1[d];
 	r2[d]+=(F2[d]-param.sigma*(distance(r1,r2)-param.ls)*((r2[d]-r1[d])/distance(r1,r2)))/param.mu*param.dts+noise2[d];
      }
    }
    
    if(fabs(r1b[0]-r1[0])<0.000001 && fabs(r1b[1]-r1[1])<0.000001 && fabs(r2b[0]-r2[0])<0.000001 && fabs(r2b[1]-r2[1])<0.000001){
      CFGcount++;
    }
    CFGcount2++;
  }

  if(CFGcount>4999)
    cout<<"division OK"<<endl;
  else if(CFGcount2>49999)
    cout<<"division error"<<endl;
}

//CFG division
float Calculation::chi(float (&r)[2],float (&r1)[2],float (&r2)[2]){
  return 0.5f*(1+tanhf(g(r,r1,r2)/param.ep_d));
}

float Calculation::g(float (&r)[2],float (&r1)[2],float (&r2)[2]){
  return
    (r1[0]-r2[0])/distance(r1,r2)*(r[0]-(r1[0]+r2[0])*0.5f)
    +(r1[1]-r2[1])/distance(r1,r2)*(r[1]-(r1[1]+r2[1])*0.5f)
    ;
}

float Calculation::H(float u){
  return u*(1.0f-u);
}

float Calculation::G(float e_eta){
  return param.rho0-param.rhoe*e_eta;
}

float Calculation::distance(float (&r1)[2],float (&r2)[2]){
  return sqrtf((r1[0]-r2[0])*(r1[0]-r2[0])+(r1[1]-r2[1])*(r1[1]-r2[1]));
}

void Calculation::division(Cell*& mcell,Cell*& dcell_2,Lumen* lumen,Interaction* field,int num_cells){

  cout<<"division"<<endl;
  dcell_2->copy_param(mcell);

  cout<<"mcell "<<mcell->Get_Gx_cell()<<" dcell "<<dcell_2->Get_Gx_cell()<<endl;
  // dcell_2->outputtofile("check_d_0.dat");
  // mcell->outputtofile("check_m_0.dat");

  float G_m[2]={};
  G_m[0]=mcell->Get_Gx_cell(); G_m[1]=mcell->Get_Gy_cell();
  //cout<<"mcell Gx "<<G_m[0]<<" mcell Gy "<<G_m[0]<<endl;

  //set the angle of the division plane
  float theta;
  if(num_cells==1 && param.init_theta>=0.0f)
    theta=(param.init_theta+90.0f)/360.0f*M_PI*2.0f;
  else
    theta=MT_rand()*M_PI*2.0f;

  //20220913 add
  if(num_cells==2 && param.init_theta>=0.0f)
    theta=(param.init_theta)/360.0f*M_PI*2.0f;

  //cerr<<"theta="<<theta/(2.0f*M_PI)*360.0f<<endl;
  cout<<"theta="<<theta/(2.0f*M_PI)*360.0f<<endl;

  int DIM = 2; //2D
  float r[2]={};
  float dr[2]={};
  float r1[2]={};  float r2[2]={};
  dr[0]=param.dx; dr[1]=param.dy;
  r1[0]=G_m[0]-cosf(theta)*dr[0]; r1[1]=G_m[1]-sinf(theta)*dr[1];
  r2[0]=G_m[0]+cosf(theta)*dr[0]; r2[1]=G_m[1]+sinf(theta)*dr[1];
  cout<<"P1xb="<<r1[0]<<" P1yb="<<r1[1]<<" P2xb="<<r2[0]<<" P2yb="<<r2[1]<<endl;
  if(num_cells!=2) //20220913 add
    time_evolution_poles(mcell,field,r1,r2,num_cells);
  cout<<"P1xa="<<r1[0]<<" P1ya="<<r1[1]<<" P2xa="<<r2[0]<<" P2ya="<<r2[1]<<endl;

  float Pc[2]={};
  Pc[0]=(float)(r1[0]+r2[0])*0.5f; Pc[1]=(float)(r1[1]+r2[1])*0.5f;
  cout<<"Pcx="<<Pc[0]<<" Pcy="<<Pc[1]<<endl;
  if(lumen->Get_rs()>0.0000001){
    mcell->seed_lumen(Pc,lumen);
    lumen->Set_volume(lumen->volume_h());
  }
  //dcell_2->outputtofile("check_d_1.dat");
  //mcell->outputtofile("check_m_1.dat");

  float rmin[2]={};
  rmin[0]=G_m[0]-param.xmax_cell*0.5f*dr[0];
  rmin[1]=G_m[1]-param.ymax_cell*0.5f*dr[1];
  //cout<<"xmin="<<rmin[0]<<" ymin="<<rmin[0]<<endl;

  //debag------------------------cfg
  // ofstream foutchig("chi-g.dat"); 
  // if(!foutchig.is_open()){cerr<<"ERROR:Could not open output file lumen data"<<endl;exit(8);}
  //debag------------------------end

  for(int i=0;i<mcell->Get_xmax();i++){
    r[0]=rmin[0]+(i+0.5f)*dr[0];
    for(int j=0;j<mcell->Get_ymax();j++){
      r[1]=rmin[1]+(j+0.5f)*dr[1];

      //random division------------------------
      // dcell_2->Set_value(i,j,mcell->Get_value(i,j)*chi(i,j,param.xmax_cell,param.ymax_cell,1.0f,theta));
      // mcell->Set_value(i,j,mcell->Get_value(i,j)*chi(i,j,param.xmax_cell,param.ymax_cell,-1.0f,theta));
      //random division------------------------end

      //CFG division---------------------------
      dcell_2->Set_value(i,j,mcell->Get_value(i,j)*(1.0f-chi(r,r1,r2)));
      mcell->Set_value(i,j,mcell->Get_value(i,j)*(chi(r,r1,r2)));
      //CFG division---------------------------end

      //debug------------------------CFG
      //cout<<r[0]<<" "<<r[1]<<" "<<chi(r,r1,r2)<<endl;
      //foutchig<<r[0]<<" "<<r[1]<<" "<<chi(r,r1,r2)<<" "<<g(r,r1,r2)<<endl;
      //debug------------------------end
    }
    //debug------------------------CFG
    //foutchig<<endl;
    //debug------------------------end
  }

  //dcell_2->outputtofile("check_d_2.dat");
  //mcell->outputtofile("check_m_2.dat");

  mcell->Set_volume(mcell->volume_h());
  dcell_2->Set_volume(dcell_2->volume_h());
  if(param.alpha_V>0.0f){
    dcell_2->Set_max_volume(mcell->Get_max_volume()*0.5f);
    mcell->Set_max_volume(mcell->Get_max_volume()*0.5f);
  }
  else{
    mcell->Set_max_volume(param.V);
    dcell_2->Set_max_volume(param.V);
  }

  mcell->CoM();
  dcell_2->CoM();

  //add cell_polarity 20191216-------------------------start 
  dcell_2->Set_pol(mcell->Get_pol());
  //add cell_polarity 20191216-------------------------end
}



// void Calculation::time_evolution_cellpolarity(Cell *&cell,VInteraction* polfield,Interaction* field,Lumen* lumen){
//   int imin=cell->Get_imin();
//   int jmin=cell->Get_jmin();
//   vector3 v;
//   v.set(0.0f,0.0f,0.0f);

//   for(int i=1;i<cell->Get_xmax()-1;i++){
//     for(int j=1;j<cell->Get_ymax()-1;j++){
//       // v+=param.kappa_n*h(cell->Get_value(i,j))*field->Get_value(i+imin,j+jmin)*(cell->Get_pol()+polfield->Get_vector(i+imin,j+jmin))
//       // 	-param.lambda_n*h(cell->Get_value(i,j))*lumen->nabla_h(i+imin,j+jmin)
//       v+=param.kappa_n*(cell->nabla_h(i,j)%field->nabla_h(i+imin,j+jmin))*(cell->Get_pol()+polfield->Get_vector(i+imin,j+jmin))
//        	-param.lambda_n*h(cell->Get_value(i,j))*lumen->nabla_h(i+imin,j+jmin)
// 	;
//     }
//   }
//   v=cell->Get_pol()+param.dt*v;
//   cell->Set_pol_tmp(v/((float)v.abs()));
// }


//-----------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------

void timestamp(const int& time){
  string fname;
  fname=param.Dir+"/timestamp.dat";
  ofstream fout(fname.c_str()); 
  if(!fout.is_open()){cerr<<"ERROR:Could not open output file "<<fname<<endl;exit(8);}
  //cout<<"timestamp "<<time<<endl;
  fout<<time<<endl;
  fout.close();
}

int set_init_t(){
  int t;
  ifstream fin("STATE/timestamp.dat"); 
  if(!fin.is_open()){cerr<<"ERROR:Could not open input file "<<endl;exit(8);}
  fin>>t;
  fin.close();
  return t;
}


//-----------------------------------------------------------------------
// main
//-----------------------------------------------------------------------

int main(int argc, char *argv[]){

  string fname;
  string ss;

  cout<<"set parameter start..."<<endl;
  param.paraminput(argv);  //parameter input
  cout<<"Done"<<endl;

  int sign;
  if(param.tau_c!=0.0f && param.r_s!=0.0f)
    sign=3;// cells & lumen & ecm
  else if(param.tau_c==0.0f && param.r_s!=0.0f)
    sign=2;// cells & lumen
  else if(param.tau_c==0.0f && param.r_s==0.0f)
    sign=1;// cells
  else if(param.tau_c!=0.0f && param.r_s==0.0f)
    sign=4;
  cout<<sign<<endl;

  //Mersenne twister---------------------------start
  //int seed=(unsigned)time(NULL);
  //my_srand(seed);
  //my_srand(param.rseed);

  string str_seed=to_string((int)param.rseed);
  str_seed+=to_string((int)param.td);
  str_seed+=to_string(param.xi0*100);
  int seed=stoi(str_seed);
  cout<<"str_seed "<<str_seed<<" outseed "<<seed<<endl;

  my_srand(seed);
  //randseed(seed);
  //Mersenne twister---------------------------end

  int num_cells=param.init_num_cells;
  int max_num_cells=param.MAX_num_cells;
  cout<<"num_cells "<<num_cells<<endl;

  Cell** cells;
  cells = new Cell*[max_num_cells];
  for(int n=0;n<max_num_cells;n++){
    cells[n] = new Cell(param.xmax_cell,param.ymax_cell);
  }
  Lumen* lumen = new Lumen(param.xmax,param.ymax);
  ECM* ecm = new ECM(param.xmax,param.ymax);

  //initial condition
  int init_t;
  cout<<"set initial phase field start...";
  init_t=0;
  // switch(sign){
  // case 3:
  // case 2:
  lumen->init_lumen();
  // case 1:
  // case 4:
  for(int n=0;n<num_cells;n++){
    cout<<" n="<<n<<endl;
    cells[n]->init_cell(n);
  }
  //}
  cout<<" Done"<<endl;

  //read wall
  // cout<<"set wall start...";
  // Wall* wall = new Wall(param.xmax,param.ymax);
  //wall->Set_wall_circle();
  //wall->Set_wall_input();
  // wall->init_field();
  // cout<<"Done"<<endl;

  //interaction field phi & e_eta
  cout<<"set interaction field start...";
  Interaction* phi;
  Interaction* e_eta;
  phi = new Interaction(param.xmax,param.ymax);
  e_eta = new Interaction(param.xmax,param.ymax);

  Interaction* phi2;
  phi2 = new Interaction(param.xmax,param.ymax);

  Interaction* u_adhe;
  u_adhe = new Interaction(param.xmax,param.ymax);

  //output field u_all
  Interaction* u_all;
  u_all = new Interaction(param.xmax,param.ymax);

  bool leakage=false;

  //init field u_all
  u_all->init_field();
  u_all->Set_u_all(cells,num_cells);

  //init field phi
  phi->init_field();
  phi->Set_phi(cells,num_cells);

  //init u adhesion field u_adhe
  u_adhe->init_field();
  u_adhe->Set_e_eta(cells,num_cells);

  // switch(sign){
  // case 3:
  // case 4:
  ecm->init_ecm(u_all);
  ecm->boundary_dirichlet(1.0f);
  //}

  cout<<"Done"<<endl;

  cout<<param.Dir<<endl;

  //output initial condition--------------------------------------start
  cout<<"output t=0...";

  fname=param.Dir+"/u_0.dat";
  u_all->outputtofile(fname);

  fname=param.Dir+"/s_0.dat";
  lumen->outputtofile(fname);

  //fname=param.Dir+"/c_0.dat";
  //ecm->outputtofile(fname);

  cout<<"Done"<<endl;
  //output initial condition--------------------------------------end

  //output file number of cell & center
  fname=param.Dir+"/numofcell.dat";
  ofstream fout_nc(fname.c_str()); 
  if(!fout_nc.is_open()){cerr<<"ERROR:Could not open output file numofcell"<<endl;exit(8);}
  fout_nc<<"0 "<<num_cells<<" ";
  for(int n=0;n<num_cells;n++)
    fout_nc<<cells[n]->Get_Gx_cell()<<" "<<cells[n]->Get_Gy_cell()<<" ";
  for(int n=0;n<max_num_cells-num_cells;n++)
    fout_nc<<"-999 -999 ";
  fout_nc<<endl;

  //output file cell volume
  fname=param.Dir+"/Volume.dat";
  ofstream fout_vol(fname.c_str()); 
  if(!fout_vol.is_open()){cerr<<"ERROR:Could not open output file volume"<<endl;exit(8);}
  fout_vol<<"0 "<<lumen->Get_volume()<<" "<<ecm->Get_volume()<<" ";
  for(int n=0;n<num_cells;n++)
    fout_vol<<cells[n]->Get_volume()<<" "
	    <<cells[n]->Get_max_volume()<<" "
	    <<cells[n]->Get_td()<<" "
	    <<cells[n]->Get_time()<<" ";
  for(int n=0;n<max_num_cells-num_cells;n++) fout_vol<<"0.0 0.0 0.0 0.0 ";
  fout_vol<<endl;


  cout<<"calculation start..."<<endl;
  int out_dt=(int)(param.sampling/param.dt);
  cout<<"out_dt="<<out_dt<<endl;
  for(int t=1+init_t;t<(int)(param.tmax/param.dt)+1;t++){
    //cout<<"t="<<t<<" start"<<endl;

    //time evolution-----------------------------------------start

    for(int n=0;n<num_cells;n++){
      phi2->init_field();
      phi2->Set_phi2(cells,num_cells,n);
      //phi->minuscell(cells[n]);
      cells[n]->time_evolution(lumen,phi2);
      //cells[n]->time_evolution(lumen,phi2,ecm);
      //phi->pluscell(cells[n]);
    }

    lumen->time_evolution(phi);
    //lumen->time_evolution(phi,ecm);
    //lumen->time_evolution(phi,u_all,ecm);

    ecm->time_evolution(phi,lumen);

    //time evolution-----------------------------------------end

    //cout<<"out0"<<endl;

    //update data--------------------------------------------start

    for(int n=0;n<num_cells;n++){
      cells[n]->update_cell();
      cells[n]->CoM();
      if(param.alpha_V>0.0f) cells[n]->time_evolution_V();
      cells[n]->add_time();
    }
    lumen->update_lumen();
    ecm->update_ecm();

    //update field phi
    phi->init_field();
    phi->Set_phi(cells,num_cells);

    //update data--------------------------------------------end

    //cout<<"out1"<<endl;

    //reshaping----------------------------------------------start
    // if((sign==3 || sign==4) && param.tol_c!=0.0f)
    //   ecm->intermediate_step();
    // if(sign!=1 && param.tol_s!=0.0f)
    //   lumen->intermediate_step();
    // if(param.tol_u!=0.0f){
    //   for(int n=0;n<num_cells;n++){
    // 	cells[n]->intermediate_step();
    //   }
    // }
    //reshaping----------------------------------------------end

    //cell division------------------------------------------start
    for(int n=0;n<num_cells;n++){

      //check (volume of cells > volume of division) & (number of cells < maxmum number of cells)
      if(cells[n]->Get_volume()>param.V-param.v_d && num_cells<max_num_cells){

      //check (volume of cells > volume condition for division) 
      //    & (cell time < time condition for division)
      //    & (number of cells < maxmum number of cells)
      // if(cells[n]->Get_volume()>param.V-param.v_d
      // 	 && cells[n]->Get_time()>cells[n]->Get_td()
      // 	 && num_cells<max_num_cells){

	//cal e_eta
	if(e_eta->Get_e_eta_null()==true){
	  e_eta->reset_e_eta();
	  e_eta->Set_e_eta(cells,num_cells);
	  // stringstream ss;
	  // ss<<t;
	  // fname=param.Dir+"/e_eta_"+ss.str()+".dat";
	  // e_eta->outputtofile(fname);
	}

	//cell division
	//cells[num_cells]=division(cells[n],lumen,e_eta,num_cells);
	cal.division(cells[n],cells[num_cells],lumen,e_eta,num_cells);
	//cal.division(cells[n],cells[num_cells],pcxs[n],pcxs[num_cells],e_eta,num_cells);

	//update field phi
	phi->init_field();
	phi->Set_phi(cells,num_cells);

	num_cells++;
	//cerr<<"division end t="<<t<<endl;
	cout<<"division end t="<<t<<endl;
      }
    }

    //reset e_eta
    if(e_eta->Get_e_eta_null()==false){
      e_eta->reset_e_eta();
    }
    //cell division---------------------------------------------------end

    u_all->init_field();
    u_all->Set_u_all(cells,num_cells);

    u_adhe->init_field();
    u_adhe->Set_e_eta(cells,num_cells);


    //output------------------------------------------------------------------start
    if(t%out_dt==0){
      int output_time=t*param.dt;
      timestamp(output_time);

      ss=to_string(output_time);

      //all cells
      fname=param.Dir+"/u_"+ss+".dat";
      u_all->outputtofile(fname);

      //lumen
      fname=param.Dir+"/s_"+ss+".dat";
      lumen->outputtofile(fname);


      //medium
      // fname=param.Dir+"/c_"+ss+".dat";
      // ecm->outputtofile(fname);

      //output volumes
      fout_vol<<output_time<<" "
	      <<lumen->Get_volume()<<" "
	      <<ecm->Get_volume()<<" ";
      for(int n=0;n<num_cells;n++)
	fout_vol<<cells[n]->Get_volume()<<" "
		<<cells[n]->Get_max_volume()<<" "
		<<cells[n]->Get_td()<<" "
		<<cells[n]->Get_time()<<" ";
      for(int n=0;n<max_num_cells-num_cells;n++)
	fout_vol<<"0.0 0.0 0.0 0.0 ";
      fout_vol<<endl;

      //output number of cell & Gx,Gy
      fout_nc<<output_time<<" "<<num_cells<<" ";
      for(int n=0;n<num_cells;n++)
      	fout_nc<<cells[n]->Get_Gx_cell()<<" "<<cells[n]->Get_Gy_cell()<<" ";
      for(int n=0;n<max_num_cells-num_cells;n++)
      	fout_nc<<"-999 -999 ";
      fout_nc<<endl;

      //check leakage-----------------------------------------------------------start
      leakage=ecm->check_leakage(lumen);
      if(leakage) t=(int)(param.tmax/param.dt)+1;
      //check leakage-----------------------------------------------------------end

    }
    //output------------------------------------------------------------------end

  }
  cout<<"Done"<<endl;

  delete[] cells;
  delete lumen;
  delete ecm;
  delete phi;
  delete phi2;
  delete u_all;
  delete e_eta;
  delete u_adhe;

  cout<<"delete end"<<endl;

  fout_vol.close();
  fout_nc.close();
  return 0;
}

