#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
using namespace std;

inline double sq(double a){
  return a*a;
}


int main(int argc, char *argv[])
{
  const int siz=atof(argv[1]);
  const double kappa=1;
  int nproc;
  int rank;
  MPI_Status Stat[4];
  MPI_Request request[4];
  /*MPI_Request request2;
  MPI_Request request3;
  MPI_Request request4;*/
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  double** grid = new double*[siz/nproc + 2];
  for (int i=0;i<siz;i++){
    grid[i]=new double[siz];
  }


    for(int i=0;i<siz/nproc + 2;i++){
	    for(int j=0;j<siz;j++){	
		    grid[i][j]=0;
	    }
    }

  for(int i=0;i<siz/nproc;i++){
    //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    grid[i+1][siz-1]=sq(sin((i + rank*siz/nproc)*M_PI/siz));
  }

  for(int i=0;i<siz/nproc;i++){
    //      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      grid[i+1][0]=sq(cos((i + rank*siz/nproc)*M_PI/siz));
  }

  double* ghostf_s=new double[siz];
  double* ghostb_s=new double[siz];
  double* ghostf_r=new double[siz];
  double* ghostb_r=new double[siz];
  
  double** grid_n = new double*[siz/nproc +2];
  for (int i=0;i<siz/nproc +2;i++){
     grid_n[i]=new double[siz];
  }

  double dx=M_PI/siz;
  const double dt=sq(dx)/(8*kappa);
  const double time=0.5*sq(M_PI)/kappa;
  const double nsteps=time/dt;
  //Set boundary conditions - fixed
   
  cout<<"Size is"<<siz<<endl;
  
    for(int i=0;i<siz/nproc +2;i++){
	    for(int j=0;j<siz;j++){	
		    grid_n[i][j]=grid[i][j];
	    }
    }

    


  for(int tt=0;tt<nsteps;tt++){
	  for(int i=1;i<siz/nproc +2-1;i++){
		  for(int j=1;j<siz-1;j++){	
			  grid_n[i][j]=grid[i][j]+ kappa*dt*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1]-4*grid[i][j])/sq(dx);
		  }
	  }

    for(int i=0;i<siz/nproc +2;i++){
	    for(int j=0;j<siz;j++){	
		    grid[i][j]=grid_n[i][j];
	    }
    }



  for(int i=0;i<siz;i++){
    ghostb_s[i]=grid[1][i];
    ghostf_s[i]=grid[siz/nproc][i];
  }

  int rf=(rank+1)%nproc;
  int rb=(rank-1)%nproc;
  if(rb<0) rb=nproc-1;

  int tag1=1;
  int tag2=2;

  //cout<<rank<<" "<<nproc<<" "<<rb<<" "<<rf<<endl;
  cout<<ghostb_r[0]<<" "<<ghostf_s[0]<<endl;
 
  MPI_Isend(&ghostf_s, siz, MPI_DOUBLE, rf , tag1, MPI_COMM_WORLD, &request[0]);
  cout<<"Proc " <<rank<<" sending to  "<<rf<<" with tag "<<tag1<<endl;
  MPI_Isend(&ghostb_s, siz, MPI_DOUBLE, rb , tag2, MPI_COMM_WORLD, &request[1]);
  MPI_Irecv(&ghostf_r, siz, MPI_DOUBLE, rf , tag2, MPI_COMM_WORLD, &request[2]);
  MPI_Irecv(&ghostb_r, siz, MPI_DOUBLE, rb , tag1, MPI_COMM_WORLD, &request[3]);
  cout<<"Proc " <<rank<<" recv from  "<<rb<<" with tag "<<tag1<<endl;

  //  MPI_Wait(&request1,MPI_STATUS_IGNORE);
  //MPI_Wait(&request2,MPI_STATUS_IGNORE);
  //MPI_Wait(&request3,MPI_STATUS_IGNORE);
  //  MPI_Wait(&request4,MPI_STATUS_IGNORE);
  MPI_Waitall(4, request, Stat);
  cout<<"Received"<<endl;
    cout<<ghostb_r[0]<<" "<<ghostf_s[0]<<endl;

  for(int i=0;i<siz;i++){
        cout<<rank<<"in for loop "<<ghostf_r[i]<<" "<<ghostb_r[i]<<endl;
    grid[0][i]= ghostf_r[i];
    grid[siz/nproc+1][i]=ghostb_r[i];
  }

  cout<<"Copied"<<endl;
  }
  
  char filename[50];
  sprintf(filename,"map_mpi%d.txt",rank);
  ofstream fout(filename);
  
    for(int i=0;i<siz/nproc +2;i++){
	    for(int j=0;j<siz;j++){	
	      fout<< i<<" "<<j<<" "<< grid[i][j]<<endl;
	    }fout<<endl;
    }

    fout.close();
    MPI_Finalize();
}
