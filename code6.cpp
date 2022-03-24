//Wrong B.C.

//Refference: momentum equation, continuity equation (page33, equa: 3.29)
//use gauss seidel iteration to solve Pressure
//Use Runge-Kutta to solve for u and v
//find u,v using equa:3.31, P33


#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <string>

using namespace std;

void initialize(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_old, double **P_new, int nx, int ny){

  for(int i=0; i<=nx; i++){
    for(int j=0; j<=ny; j++){
      P[i][j] =	0;
      P_old[i][j] = 0;
      P_new[i][j] = 0;
            
    }
  }
  for(int i=0; i<=nx; i++){
    for(int j=0; j<=ny-1; j++){
      G[i][j] = 0;
      v[i][j] =	0;
      v_new[i][j] = 0;
      
    }
  }
  for(int i=0; i<=nx-1; i++){
    for(int j=0; j<=ny; j++){
      F[i][j] = 0;
      u[i][j] = 0;
      u_new[i][j] = 0;  
    }
  }

  // Boundary Condition of u
  for(int j=1; j<=ny-1; j++){
    u[0][j] = 1;
    //u[i][1] = 0.5;
    
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


//step 1 : find F&G
//Separate Runge-Kutta (Pg34, equa 3.36) into 3 terms: RHS_X1, RHS_X2 and RHS_X3

void simulation_FG(double **u, double **u_new, double **v, double **v_new, double **F, double **G, int nx, int ny, double Re, double dx, double dy, double dt){
  double RHS_F;
  double RHS_F1;
  double RHS_F2;
  double RHS_F3;

  double RHS_G;
  double RHS_G1;
  double RHS_G2;
  double RHS_G3;

  
  for(int i=1; i<=nx-2; i++){
    for(int j =1; j<=ny-1; j++){
      RHS_F1 = (((u[i+1][j]+2*u[i][j]+u[i-1][j])/pow(dx,2.))+((u[i][j+1]+2*u[i][j]+u[i][j-1])/pow(dy,2.)))/Re;
      RHS_F2 = (pow(((u[i][j]+u[i+1][j])/2),2.)-pow(((u[i-1][j]+u[i][j])/2),2.))/dx;
      RHS_F3 = (((v[i][j]+v[i+1][j])/2)*((u[i][j]+u[i][j+1])/2)-((v[i][j-1]+v[i+1][j-1])/2)*((u[i][j-1]+u[i][j])/2))/dy;
      RHS_F = RHS_F1  -RHS_F2- RHS_F3;
      F[i][j] = u[i][j] + dt*RHS_F;
    }
  }
  for(int j =1; j<=ny-1; j++){
   F[0][j]=1;
  }
 

  for(int i=1; i<=nx-1; i++){
    for(int j =1; j<=ny-2; j++){
      RHS_G1 = (((v[i+1][j]+2*v[i][j]+v[i-1][j])/pow(dx,2.))+((v[i][j+1]+2*v[i][j]+v[i][j-1])/pow(dy,2.)))/Re;
      RHS_G3 = (pow(((v[i][j]+v[i][j+1])/2),2.)-pow(((v[i][j-1]+v[i][j])/2),2.))/dy;
      RHS_G2 = (((u[i][j]+u[i][j+1])/2)*((v[i][j]+v[i+1][j])/2)-((u[i-1][j]+u[i-1][j+1])/2)*((v[i-1][j]+v[i][j])/2))/dx;
      RHS_G = RHS_G1 - RHS_G2 -RHS_G3 ;
      G[i][j] = v[i][j] + dt*RHS_G;
    }
  }
  // Neumann Boundary Condition setting for G and F
  for(int j =1; j<=ny-1; j++){
   G[nx-1][j]=G[nx-2][j];
  }
 for(int j =1; j<=ny-2; j++){
   F[nx-1][j]=F[nx-2][j];
  }

}


//step 2 find P[i,j] using gauss seidel iteration
void simulation_p(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_new, int nx, int ny, double Re, double dx, double dy, double dt){
  double aa, bb, cc;
  bool status = true;
  aa = 1/pow(dx,2.);
  bb = 1/pow(dy,2.);
  cc = (pow(dx,2.)*pow(dy,2.))/(2*(pow(dx,2.)+pow(dy,2.)));

  while(true){
    if(status == false){break;}
    for(int i=1; i<=nx-1; i++){
      for(int j =1; j<=ny-1; j++){
	P[i][j]=P_new[i][j];
      }
    }

    //For pressure, use equa 3.38, pg35, to set from i =1 to i = nx-1

    for(int i=1; i<=nx-1; i++){
      for(int j =1; j<=ny-1; j++){
	P_new[i][j]=((((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy)*(-1/dt))+(aa*P[i+1][j]+aa*P_new[i-1][j]+bb*P[i][j+1]+bb*P_new[i][j-1]))*cc;
	status = false;
      }
    }

    for(int i=1; i<=nx-1; i++){
      for(int j =1; j<=ny-1; j++){
	if(abs(P[i][j] - P_new[i][j]) > 10e-50){
	  status = true;
	}
      }
    }
 
  }  
}

//step3 find u, v using Runge-Kutta
void simulation_uv(double **u, double **u_new, double **v, double **v_new, double **F, double **G, double **P, double **P_new, int nx, int ny, double Re, double dx, double dy, double dt){

  for(int i=1; i<=nx-2; i++){
      for(int j =1; j<=ny-1; j++){
	u_new[i][j]=F[i][j]+(dt*(P_new[i+1][j]-P_new[i][j]))/dx;
      }
  }
	
  for(int i=1; i<=nx-1; i++){
      for(int j =1; j<=ny-2; j++){
	v_new[i][j]=G[i][j]+(dt*(P_new[i][j+1]-P_new[i][j]))/dy;
      }
      	
    }
 for(int j =1; j<=ny-1; j++){
   u_new[0][j]=1;
  }
 for(int j =1; j<=ny-1; j++){
   u_new[nx-1][j]=u_new[nx-2][j];
  }
 for(int j =1; j<=ny-2; j++){
   v_new[nx-1][j]=v_new[nx-2][j];
  }

 
}

void update(double **var, double **var_new, int nx, int ny){
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){    
      var[i][j] = var_new[i][j];
    }
  }
}

void paraview(string fileName, double **var, int nx, int ny, double dx, double dy){
  ofstream myfile;
  myfile.open(fileName);
  //------------------------------------------------------------//
    // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

    // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for(int j = 0; j <= ny-1; j++){
    for(int i = 0; i <= nx-1; i++){
      myfile << dx*i << " " << dy*j << " 0\n";
    }
  }
  
  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 0; j <= ny-1; j++){
    for(int i = 0; i <= nx-1; i++){
      myfile << var[i][j] << "\n";
    }
  }
  myfile.close();

}


int main(){

  int nx = 100;
  int ny = 50;
  double Re = 300.;
  double dx = 3;
  double dy = 1;
  double dt = 0.001;
  string fileName;
  
  double **u;
  u = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u[i] = (double *) malloc((ny+1) * sizeof(double));
  }
  double **v;
  v = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    v[i] = (double *) malloc(ny * sizeof(double));
  }
  double **u_new;
  u_new = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    u_new[i] = (double *) malloc((ny+1) * sizeof(double));
  }
  double **v_new;
  v_new = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    v_new[i] = (double *) malloc(ny * sizeof(double));
  }

  double **F;
  F = (double **) malloc (nx * sizeof(double));
  for(int i=0; i<nx; i++){
    F[i] = (double *) malloc((ny+1) * sizeof(double));
  }
  double **G;
  G = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    G[i] = (double *) malloc(ny * sizeof(double));
  }

  double **P;
  P = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    P[i] = (double *) malloc((ny+1) * sizeof(double));
  }
  double **P_old;
  P_old = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    P_old[i] = (double *) malloc((ny+1) * sizeof(double));
  }
  double **P_new;
  P_new = (double **) malloc ((nx+1) * sizeof(double));
  for(int i=0; i<nx+1; i++){
    P_new[i] = (double *) malloc((ny+1) * sizeof(double));
  }

  // visualize(u,nx,ny);
  //visualize(v,nx,ny);
  
  initialize(u,u_new,v,v_new,F,G,P,P_old,P_new,nx,ny);
for (int n = 1; n<=30; n++){  
  simulation_FG(u, u_new, v, v_new, F, G, nx, ny, Re, dx, dy, dt);
  simulation_p(u, u_new, v, v_new, F, G, P, P_new, nx, ny, Re, dx, dy, dt);
  simulation_uv(u, u_new, v, v_new, F, G, P, P_new, nx, ny, Re, dx, dy, dt);
  update(u,u_new,nx,ny+1);
  update(v,v_new,nx+1,ny);
  update(P,P_new,nx+1,ny+1);
  
  cout << "n = " << n << "\n";
 
  //if( n%5 == 0) {
  
  fileName = "phi_" + to_string(n) + ".vtk";
  paraview(fileName, u, nx, ny+1, dx, dy);

  //visualize(F,nx,ny+1);
  
  visualize(u,nx,ny+1);

  //visualize(G,nx+1,ny);
  //visualize(v,nx+1,ny);
  //visualize(P,nx+1,ny+1);
  }
//}
  
  
}
