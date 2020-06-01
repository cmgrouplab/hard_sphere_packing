//generating random packing models of multi-component alloys using 
// adpative shrinking cell MC method 
//The atoms are modeled as hard spheres with different radii

//started on 05/19/19
//author: yang.jiao.2@asu.edu; zhuanghl@asu.edu


//allows the users to 
// (1) specify the total number of atoms
// (2) specify the number of compoenets
// (3) specify the molar fraction/number of each component
// (4) specify the radius of each compnent 

// automatically compute the box length for a default initial density, e.g., 0.1
// use MC + ASC to density the system 



//*********************************************

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iomanip>


#define MAXX 5000  //a large number for random number generator
#define TOL 0.0001 //a small number .... 


//parameters defining the atomic system...

double L_box = 5.0; // the edge length of the cubic box, will be rescaled 
int N_atom; // the total number of atoms
int n_c; //the total number of components 

int* N_c; //the number of atoms for each species [n_c]
double* R_c; //the radius of atoms for each specifies [n_c]
double* V_c; //the volume of the spheres

//the coordinates of the atoms
double* x; //[N_atom]
double* y;
double* z;

//the radius of each atoms
double* R_atom; //[N_atom]


//*******************************************************
//parameters for the ACS MC simulations

int N_stage = 2000; //totoal number of MC stages, a volume compression is applied at the end of each stage
int N_step = 1000; //total number of MC moves per stage for partilce moves; after N_step particle move, the boundary box deform

double delta_L_ratio = 0.1; //L_box is attempted to be re-scaled by (1-delta_L); if not successfull, delta_L is reduced by half repeated until success is achieved

double delta_R_ratio = 0.05; //the ratio for magnitude of particle displacement, meausred in terms of the box length
                             //this is most convenient as coordinates are resclaed w.r.t. box length
//double delta_R; //the magnitude of atom displacement, = delta_R_ratio*R_ave, to be computed later

//for restroing the rejected moves
int MC_index;
double MC_x;
double MC_y;
double MC_z; 

//for restroing rejected deformation
double MC_L_box;


//for input and output
ifstream fin;
ofstream fout;

//******************************************************
//init all the data structures, 
//get parameters from standard run-time input
void get_parameters()
{
     cout<<"Please specify the number of components n_c = :"; 
     cin>>n_c;
     
     N_c = new int[n_c]; //the array for number of atoms of each component
     R_c = new double[n_c]; //the array for the radius of each component
     V_c = new double[n_c];
     
     cout<<"Please specify the number of atoms for each components:"<<endl;
     for(int i=0; i<n_c; i++)
        {
             cout<<"N_c["<<i<<"] = ";
             cin>>N_c[i];
        }
     cout<<endl;
     
     cout<<"Please specify the radius of atoms for each components:"<<endl;
     for(int i=0; i<n_c; i++)
        {
             cout<<"R_c["<<i<<"] = ";
             cin>>R_c[i];
             
             V_c[i] = 4.0*3.1415926*R_c[i]*R_c[i]*R_c[i]/3.0;
             cout<<"V_c[i] = "<<V_c[i]<<endl;
        }
     cout<<endl;
     
     //compute the total number of atoms
     N_atom = 0;
     for(int i=0; i<n_c; i++)
        N_atom += N_c[i];
        
     cout<<"The total of number of atoms N_atom = "<<N_atom<<endl;
     
     //get the arrays for the particle radius and positions
     for(int i=0; i<N_atom; i++)
        {
             x = new double[N_atom];
             y = new double[N_atom];
             z = new double[N_atom];
             
             R_atom = new double[N_atom];
        }
     //these arrays will be initialized when generating the initial packing configurations
     
}


double get_dist(int index1, int index2)
{
       //need to covert the relative coordinates into absolute coordinates
       //need to consider periodic boundary conditions
       double dx = fabs(x[index1] - x[index2]); 
       if(dx >= 0.5) dx = 1.0 - dx;
       
       double dy = fabs(y[index1] - y[index2]); 
       if(dy >= 0.5) dy = 1.0 - dy;
       
       double dz = fabs(z[index1] - z[index2]); 
       if(dz >= 0.5) dz = 1.0 - dz;
       
       //now rescale the dist w.r.t box length 
       double dist = sqrt(dx*dx + dy*dy + dz*dz)*L_box;
       
       return dist;
}


int check_overlap_all()
{
    //compare all pair distances
    int overlap_flag = 0;
    double temp_dist;
    
    for(int i=0; i<N_atom; i++)
       for(int j=0; j<i; j++)
          {
               temp_dist = get_dist(i, j);
               //cout<<"dist "<<i<<" "<<j<<" = "<<temp_dist<<endl;
               
               if(temp_dist < (R_atom[i] + R_atom[j]))
                  {
                      overlap_flag = 1;
                      return 1;
                  }
          }
          
    return overlap_flag;
}

void print_dist()
{
     for(int i=0; i<N_atom; i++)
       for(int j=0; j<i; j++)
          {
               cout<<"dist_"<<i<<"_"<<j<<" = "<<get_dist(i, j)<<endl;
          }
}

double max_dist()
{
     double temp_max = 0.0;  
       
     for(int i=0; i<N_atom; i++)
       for(int j=0; j<i; j++)
          {
               if(temp_max < get_dist(i, j))
                  temp_max = get_dist(i, j);
          }
     
     return temp_max;
}

double min_dist()
{
     double temp_min = 10000000000.0;  
       
     for(int i=0; i<N_atom; i++)
       for(int j=0; j<i; j++)
          {
               if(temp_min > get_dist(i, j))
                  temp_min = get_dist(i, j);
          }
     
     return temp_min;
}

double compute_density()
{
     double sum_V = 0;
     
     for(int i=0; i<n_c; i++)
         sum_V += V_c[i]*N_c[i];
         
     return sum_V/(L_box*L_box*L_box);
}

void init_packing()
{
     //just randomly place particles in the box and rescale if overlap is detected 
     for(int i=0; i<N_atom; i++)
       {
           x[i] = (double)(rand()%MAXX)/(double)MAXX;
           y[i] = (double)(rand()%MAXX)/(double)MAXX;
           z[i] = (double)(rand()%MAXX)/(double)MAXX;
       }
       
     //now we assign the radius to the atoms
     int temp_counter = 0;
     
     for(int i=0; i<n_c; i++)
        for(int j=0; j<N_c[i]; j++)
           {    
                R_atom[temp_counter] = R_c[i];
                
                temp_counter++;
            }   
     if(temp_counter != N_atom)
       {
          cout<<"acounting for assigning particle radius is wrong! Check again!"<<endl;
          cout<<"temp_counter = "<<temp_counter<<endl;
          cout<<"N_atom = "<<N_atom<<endl;
          
          exit(1);
       }
       
     int overlap_flag = check_overlap_all();
     
     while(overlap_flag == 1)
       {
          L_box = 1.05*L_box;
          cout<<"L_box = "<<L_box<<endl;  
          
          overlap_flag = check_overlap_all();              
       }
       
       
     cout<<"packing density phi = "<<compute_density()<<endl;
     
}


void displace_atom(int temp_index)
{
   //save the original position 
   MC_x = x[temp_index];
   MC_y = y[temp_index];
   MC_z = z[temp_index];
   
   //generate new displacement  
   double dx = ((double)(rand()%MAXX)/(double)MAXX - 0.5)*delta_R_ratio;
   double dy = ((double)(rand()%MAXX)/(double)MAXX - 0.5)*delta_R_ratio;
   double dz = ((double)(rand()%MAXX)/(double)MAXX - 0.5)*delta_R_ratio;
     
   //update the new positions  
   x[temp_index] = x[temp_index] + dx;
   if(x[temp_index]>=1.0) x[temp_index] = x[temp_index] - 1.0;
   else if(x[temp_index]<0) x[temp_index] = x[temp_index] + 1.0;  
   
   y[temp_index] = y[temp_index] + dy;
   if(y[temp_index]>=1.0) y[temp_index] = y[temp_index] - 1.0;
   else if(y[temp_index]<0) y[temp_index] = y[temp_index] + 1.0;
   
   z[temp_index] = z[temp_index] + dz;
   if(z[temp_index]>=1.0) z[temp_index] = z[temp_index] - 1.0;
   else if(z[temp_index]<0) z[temp_index] = z[temp_index] + 1.0;

}

void restore_atom(int temp_index)
{
    x[temp_index] = MC_x;
    y[temp_index] = MC_y;
    z[temp_index] = MC_z;
}

void deform_box()
{
    //save the current box size
    MC_L_box = L_box;
    
    double temp_L_ratio = delta_L_ratio;
    
    //now attempt to shrink the box, by (1-temp_L_ratio), if not successful, reduce temp_L_ratio by half
    L_box = (1 - temp_L_ratio)*L_box;
    
    int overlap_flag = check_overlap_all();
    
    while(overlap_flag == 1)
      {
          //restore_box();
          L_box = MC_L_box;
          
          //reduce the rescaling size and try again
          temp_L_ratio = 0.95*temp_L_ratio;
          
          L_box = (1 - temp_L_ratio)*L_box;
          
          overlap_flag = check_overlap_all();
      }
     
}

/*
void restore_box()
{
     L_box = MC_L_box;
}
*/



int MC_move()
{
    //random particle displacement
    //int temp_index;
    int temp_ct = 0;
    
    for(int i=0; i<N_step; i++)
       {
           MC_index = rand()%N_atom;
           
           displace_atom(MC_index);
           
           if(check_overlap_all() == 1)
             {
                restore_atom(MC_index);
                
                temp_ct ++;
             }
       }
       
     cout<<"**********************************"<<endl;
     cout<<"R_acc = "<<(double)(N_step - temp_ct)/(double)N_step<<endl;
     
     //now take care of the boundary
     deform_box();
     
     cout<<"Density = "<<compute_density()<<endl;
     cout<<"**********************************"<<endl;
       
}



void print_packing()
{
     cout<<"printing packing configuration to file ..."<<endl;
     fout.open("packing_atom.txt");
     
     fout<<N_atom<<" //total number of atoms"<<endl;
     fout<<n_c<<"  //number of species"<<endl;
     fout<<L_box<<"  //edge length of cubic box"<<endl;
     
     for(int i=0; i<N_atom; i++)
       fout<<fixed<<setprecision(6)<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\t"<<R_atom[i]<<endl;
       
     fout.close();

}

/*
//rescale the simulation box isotropically to remove any overlap
void rescale_box()
{
}

//void init_system()
//{
//     //use simple cubic lattice for initial configuration
//}
*/


//the main function
int main()
{
    srand(time(NULL));
    
    get_parameters();
    
    init_packing();
    
    //print_dist();
    
    cout<<"min_dist = "<<min_dist()<<endl;
    cout<<"max_dist = "<<max_dist()<<endl;
    cout<<"L_box = "<<L_box<<endl;
    
    //save the density increase....
    fout.open("density_profile.txt");
    
    //now start the MC moves
    for(int i=0; i<N_stage; i++)
    {       
        cout<<"Stage # "<<(i+1)<<endl;
        
        MC_move();
        
        fout<<(i+1)<<"\t"<<compute_density()<<endl;
        
    }
    fout.close();
    
    print_packing();
    
    
    int temp_a;
    cin>>temp_a;
    
    return 0;
}
