#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include<time.h>
using namespace std;

inline double sq(double a){
  return a*a;
}

int main(int argc, char *argv[])
{
  clock_t tclock;
  tclock=clock();
  double avg=0;

    if(argc!=3){
    printf("USAGE: %s <grid size> <nthreads>\n", argv[0]);
    exit(1);
  }

  const int siz=atof(argv[1]);
  const int nthreads=atof(argv[2]);
  cout<<"Size is"<<siz<<endl;
  const double kappa=1;
  //  cout<<"Number of threads="<<omp_get_max_threads();
  //  omp_set_num_threads(8);
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
   
  //cout<<"Size is"<<siz<<endl;

  
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    grid_n[i][j]=grid[i][j];
	    }
    }

    
  for(int tt=0;tt<nsteps;tt++){

#pragma omp parallel for num_threads(nthreads) //for default(none) private(ii)    
      
            for(int ii=1;ii<siz-1;ii++){
	      //       cout<<"Thread number"<<omp_get_thread_num()<<endl;
		  for(int j=1;j<siz-1;j++){	
			  grid_n[ii][j]=grid[ii][j]+ kappa*dt*(grid[ii-1][j]+grid[ii+1][j]+grid[ii][j-1]+grid[ii][j+1]-4*grid[ii][j])/sq(dx);
		  }
	    }
    #pragma omp parallel for num_threads(nthreads)
    for(int i=1;i<siz-1;i++){
	    grid_n[0][i]=grid[0][i]+ kappa*dt*(grid[siz-1][i]+grid[1][i]+grid[0][i-1]+grid[0][i+1]-4*grid[0][i])/sq(dx);
    }
    #pragma omp parallel for num_threads(nthreads)
    for(int i=1;i<siz-1;i++){
	    grid_n[siz-1][i]=grid[siz-1][i]+ kappa*dt*(grid[siz-2][i]+grid[0][i]+grid[siz-1][i-1]+grid[siz-1][i+1]-4*grid[siz-1][i])/sq(dx);
    }

    #pragma omp parallel for num_threads(nthreads)
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    grid[i][j]=grid_n[i][j];
	    }
    }

  }



   tclock=clock()-tclock;
   cout<<"Total time taken is "<<float(tclock)/CLOCKS_PER_SEC<<endl;


  char filename[50];
  sprintf(filename,"map_omp_%d.txt",siz);
  ofstream fout(filename);
  
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
	      fout<< i<<" "<<j<<" "<< grid[i][j]<<endl;
	      //	      avg+=grid[i][j];
	    }fout<<endl;
    }

    fout.close();

       for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){	
		    avg+=grid[i][j];
	    }
	    }
       cout<<"Average temperature is "<<avg/(siz*siz)<<endl;  
}
