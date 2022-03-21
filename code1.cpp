//solve momentum equation with no pressure in x-axis or F (page33-3.29)
//explicit euler in term t && central finite different interm dimention

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>

using namespace std;

void initialize(double **u, double **u_new, double **v, double **v_new, double **F, double **G, int nx, int ny){

  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-1; j++){
      u[i][j] = 0;
      u_new[i][j] = 0;
      v[i][j] =	0;
      v_new[i][j] = 0;
    }
  }
  
  for(int i=0; i<=nx-1; i++){
    u[i][0] = 1;
  }
  
}

void visualize(double **var, int nx, int ny){

  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny-1; j++){
      cout << var[i][j] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
  
}

void simulation_FG(double **u, double **u_new, double **v, double **v_new, double **F, double **G, int nx, int ny, double Re, double dx, double dy, double dt){
  double RHS;
  double RHS_1;
  double RHS_2;
  double RHS_3;
  for(int i=1; i<=nx-2; i++){
    for(int j =1; j<=nx-2; j++){
      RHS_1 = (((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dx,2.))+((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dy,2.)))/Re;
      RHS_2 = (pow(((u[i][j]+u[i+1][j])/2),2.)-pow(((u[i-1][j]+u[i][j])/2),2.))/dx;
      RHS_3 = (((v[i][j]+v[i+1][j])/2)*((u[i][j]+u[i][j+1])/2)-((v[i][j-1]+v[i+1][j-1])/2)*((u[i][j-1]+u[i][j])/2))/dy;
      RHS = RHS_1 - RHS_2 - RHS_3;
      F[i][j] = u[i][j] + dt*RHS;
    }
  }

  for(int i=1; i<=nx-2; i++){
    F[i][nx-1]=F[i][nx-2];
  }
}



int main(){

  int nx = 11;
  int ny = 11;
  double Re = 300.;
  double dx = 1;
  double dy = 1;
  double dt = 0.01;
  
  
  double **u;
  u = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u[i] = (double *) malloc(ny * sizeof(double));
  }
  double **v;
  v = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    v[i] = (double *) malloc(ny * sizeof(double));
  }
  double **u_new;
  u_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u_new[i] = (double *) malloc(ny * sizeof(double));
  }
  double **v_new;
  v_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    v_new[i] = (double *) malloc(ny * sizeof(double));
  }

  double **F;
  F = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    F[i] = (double *) malloc(ny * sizeof(double));
  }
  double **G;
  G = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    G[i] = (double *) malloc(ny * sizeof(double));
  }

  
  
  initialize(u,u_new,v,v_new,F,G,nx,ny);
  simulation_FG(u, u_new, v, v_new, F, G, nx, ny, Re, dx, dy, dt);
  visualize(F,nx,ny);
  visualize(u,nx,ny);

  
  
  
}
