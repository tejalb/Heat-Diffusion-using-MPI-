#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<time.h>
using namespace std;

inline double sq(double a){
  return a*a;
}

int main(int argc, char *argv[])
{
  clock_t tclock;
  tclock=clock();
  const int siz=atof(argv[1]);
  const double kappa=1;
  double avg=0;

  double** grid = new double*[siz];
  for (int i=0;i<siz;i++){
     grid[i]=new double[siz];
  }

  double** grid_n = new double*[siz];
  for (int i=0;i<siz;i++){
     grid_n[i]=new double[siz];
  }

  double dx=M_PI/siz;
  const double dt=sq(dx)/(8*kappa);
  const double time=0.5*sq(M_PI)/kappa;
  const double nsteps=time/dt;
  //Set boundary conditions - fixed

    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    grid[i][j]=0;
	    }
    }

  for (int i=0;i<siz;i++){
    grid[i][0]=sq(cos(i*M_PI/double(siz)));
    grid[i][siz-1]=sq(sin(i*M_PI/double(siz)));

  }
   
  cout<<"Size is"<<siz<<endl;
  
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    grid_n[i][j]=grid[i][j];
	    }
    }


  for(int tt=0;tt<nsteps;tt++){
	  for(int i=1;i<siz-1;i++){
		  for(int j=1;j<siz-1;j++){	
			  grid_n[i][j]=grid[i][j]+ kappa*dt*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1]-4*grid[i][j])/sq(dx);
		  }
	  }

    for(int i=1;i<siz-1;i++){
	    grid_n[0][i]=grid[0][i]+ kappa*dt*(grid[siz-1][i]+grid[1][i]+grid[0][i-1]+grid[0][i+1]-4*grid[0][i])/sq(dx);
    }

    for(int i=1;i<siz-1;i++){
	    grid_n[siz-1][i]=grid[siz-1][i]+ kappa*dt*(grid[siz-2][i]+grid[0][i]+grid[siz-1][i-1]+grid[siz-1][i+1]-4*grid[siz-1][i])/sq(dx);
    }

    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    grid[i][j]=grid_n[i][j];
	    }
    }

  }
  
 
   for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    avg+=grid[i][j];
	    }
    }
   
   tclock=clock()-tclock;
   cout<<"Total time taken is "<<float(tclock)/CLOCKS_PER_SEC<<endl;
   cout<<"Average temperature is "<<avg/(siz*siz)<<endl;

  char filename[50];
  sprintf(filename,"map_serial_%d.txt",siz);
  ofstream fout(filename);
  
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
	      fout<< i<<" "<<j<<" "<< grid[i][j]<<endl;
	    }fout<<endl;
    }

    fout.close();
}
