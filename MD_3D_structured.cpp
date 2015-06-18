#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <dirent.h>
# include "defs.h"
# include "force_structured.h"

using namespace std;


// random numbers using random_device option with normal distribution
/*void createInitialPosition_N_particles(std::string fileName, int N) {
      std::random_device rd, rd1;
      std::mt19937 genx(rd()),geny(rd1());
      std::ofstream outFile(fileName);
      std::normal_distribution<> d(0,1);
      for(int i=0;i<N;i++) {
      outFile<<d(genx)<<'\t'<<d(geny)<<std::endl;
      }
      outFile.close();
}
*/

// random numbers using rand function
void createInitialPosition_N_particles(std::string fileName,std::string fileName2, int N, double Lx, double Ly, double Lz) {
	double x,y,z, vx, vy, vz;
 	srand (time(NULL)); // initialize random seed
 	std::ofstream outFile(fileName);
 	std::ofstream outFile2(fileName2);
 	for(int i=0;i<N;i++) {
 		x=0.0;//*sigma;//((double) rand() / (RAND_MAX/Lx))-Lx/2.0;  // create particle position from -Lx/2 to Lx/2
		y=Ly/2.5-i*1.5;//abs(((double) rand() / (RAND_MAX/Ly))-Ly/2.0);
		z=0.0;//((double) rand() / (RAND_MAX/Lz))-Lz/2.0;
		vx= ((double) rand()/(RAND_MAX)-0.5);
		vy= ((double) rand()/(RAND_MAX)-0.5);
		vz= ((double) rand()/(RAND_MAX)-0.5);
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
		outFile2<<vx<<'\t'<<vy<<'\t'<<vz<<std::endl;
 	}
 	outFile.close();
}

std::vector<int> radialDistFunc(double XYZ[][3], double Lx,double Ly, double Lz, double dr, int N) {
    std::vector<int> rdf((int) floor(sqrt(pow(Lx/2.0,2.0)+pow(Ly/2.0,2.0)+pow(Lz/2.0,2.0)))/dr,0);
	double r;
	for(int j=0;j<N;j++) {
		r=sqrt(pow(XYZ[j][0],2.0)+pow(XYZ[j][1],2.0)+pow(XYZ[j][2],2.0));
	    rdf[(int) floor(r/dr)]+=1;                        // put each particle in a bin according to its position from origin ie. (0,0)
	}
	return rdf;
}

// forceUpdate fucntion included as force.h header file

std::random_device seed;
std::mt19937 gen{seed()};
std::normal_distribution<> R1(0.0,1.0),R2(0.0,1.0),R3(0.0,1.0);

void brownian( vector<SubData>& particle, vector<vector<mtrx3D>>& Mobility_Tnsr_tt) {

	for(int i=0;i<NrParticles;i++) 
	{
		vctr3D rand(R1(gen), R2(gen), R3(gen));
			for(int j=0;j<NrParticles;j++) 
				{
					particle[i].pos+=(Mobility_Tnsr_tt[i][j]*particle[j].frc)*(dt/(kb*T0));
				}
		particle[i].pos+=(rand*mu_sqrt*kbT_dt);
		particle[i].pos.PBC(box,rbox);

}


}


int main() {
// current date/time based on current system
   time_t now = time(0);
   struct tm *ltm = localtime(&now);
   cout << "start time"<< '\t'<< ltm->tm_hour << ":";
   cout << ltm->tm_min << ":";
   cout << ltm->tm_sec << endl;
         
int if_create_particles = xxcreate, ifrestart=xxrestart;
double tauT=0.1;         
double Temp=T0;
double shear_rate = 0.0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
int Max_Cluster_N=NrParticles;
double simu_time=dt;
int step=0, nSteps=10000, frame=10;
double vel_scale;
int if_Periodic =1;

std::cout<<cellx<<'\t'<<celly<<'\t'<<cellz<<'\t'<<F_g<<std::endl;
double  K_Energy, p_energy=0.0;
vctr3D dR, dr2;
double R, r2;
double dr=0.05; // step size for RDF calculation
// std::vector<int> RDF((int)  floor(sqrt((Lx/2)*(Lx/2)+(Ly/2)*(Ly/2)+(Lz/2)*(Lz/2)))/dr,0), RDF1((int)  floor(sqrt(Lx*Lx+Ly*Ly))/dr,0);

vector<SubData>  particle(NrParticles);
vector<vector<mtrx3D>>  Mobility_Tnsr_tt(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_tr(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_rt(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_rr(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr(NrParticles,vector<mtrx3D>(NrParticles));

// variables for mobility tensor calculation
vctr3D Rij;
double Rij2, Rij2_inv, temp, temp1, temp2, temp3, tau;
//

if(ifrestart)	{
std::string fileName=dataFileName+"/End_positions.dat";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
    }
}	
	fileName=dataFileName+"/Velocities.dat";
	std::ifstream dataFile1(fileName);

	if(!dataFile1.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading Vx, Vy, Vz Velocities"<<std::endl;
    std::string line;
     
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile1,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].vel.comp[0];
        currentLine >> particle[i].vel.comp[1];
        currentLine >> particle[i].vel.comp[2];
    }
}	
} else {

	std::string fileName="../XYZ.dat";
	std::string fileName2="../Initial_Velocities.dat";

if (if_create_particles) {
createInitialPosition_N_particles(fileName,fileName2,NrParticles,Lx,Ly,Lz);
}
//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
std::ifstream dataFile2(fileName2);

if(!dataFile.good() && !dataFile2.good() ) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    std::string line1;


    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];           
    }
}

}		
//delete all files before writing data

// following snippet taken from stakcflow link  http://stackoverflow.com/questions/11007494/how-to-delete-all-files-in-a-folder-but-not-delete-the-folder-c-linux
if (ifrestart) {
dataFileName=dataFileName_new;
}
const char *dataFileNamePointer = dataFileName.c_str();  // covnert the datafilename to a char pointer ans pass it to the snippet below which delete all files in that folder before running the simulation
if (!ifrestart) {
struct dirent *next_file;
DIR *theFolder;
char filepath[256];
theFolder = opendir(dataFileNamePointer);
while (( next_file = readdir(theFolder)) )
	{
    // build the full path for each file in the folder
    sprintf(filepath, "%s/%s",dataFileNamePointer, next_file->d_name);
    remove(filepath);
    }
//
}

/* sort particles into cells */
	

for (int i=0;i<NrParticles;i++) {
	particle[i].radius=sigma/2.0;
	if (!ifrestart) {
		//	particle[i].vel={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
		particle[i].vel=(particle[i].vel)*vel_scale;
	
	}
}

std::ofstream outFile(dataFileName+"/K_energy.dat");
std::ofstream outFile1(dataFileName+"/PE_energy.dat");
std::ofstream outFile2(dataFileName+"/Velocities.dat");
std::ofstream outFile7(dataFileName+"/End_positions.dat");
std::ofstream outFile3(dataFileName+"/Pressure_Tensor.dat");
std::ofstream outFile9(dataFileName+"/particle1.dat");
std::ofstream outFile10(dataFileName+"/particle2.dat");

outFile3<<"simu_time"<<'\t'<<"Pxx"<<'\t'<<"Pyy"<<'\t'<<"Pzz"<<'\t'<<"Pxy"<<'\t'<<"Pxz"<<'\t'<<"Pyz"<<'\t'<<"Mxy"<<'\t'<<"Mxz"<<'\t'<<"Myz"<<'\t'<<"Temp"<<std::endl;
// perfrom MD steps
/*	if (ifrestart) {
	simu_time =10.001;
	new_neighbor = TOPMAP(cellx,celly,cellz,if_Periodic,neighbor,box,(R_cut+R_shell),simu_time*shear_rate*Ly);
	std::tie(Pxx[0],Pyy[0],Pzz[0],Pxy[0],Pxz[0],Pyz[0])= forceUpdate(&p_energy ,XYZ, new_neighbor, cell, cellLength, box, (R_cut+R_shell), N, force, shear_rate, simu_time, ifshear , epsilon, sigma, rs);
}
*/
step = 0;
forceUpdate(particle, &p_energy);
simu_time =dt;
do {
	p_energy=0;	
	
	for (int i=0; i<NrParticles; i++)
	{
		for (int j=0; j<NrParticles; j++)
			{
				Rij=particle[i].pos-particle[j].pos;
				mtrx3D Rij_dyadic_Rij(Rij,Rij);
				vctr3D col1(0.0, -Rij.comp[2], Rij.comp[1]);
				vctr3D col2(Rij.comp[2],0.0,-Rij.comp[0]);
				vctr3D col3(-Rij.comp[1],Rij.comp[0],0.0);
				mtrx3D epsilon_rij(col1 , col2, col3);
				Rij2=Rij.norm2();
				Rij2_inv=1.0/Rij2;
				temp1=1.0/(8.0*M_PI*eta);
				tau = 1.0/(6.0*M_PI*eta*particle[i].radius);
				temp=temp1/(sqrt(Rij2));
				temp2=temp1/(particle[i].radius*particle[i].radius*particle[i].radius);
				temp3=temp/(2.0*Rij2);
				if(i==j) {
				Mobility_Tnsr_tt[i][j]	=	 	Unit_diag * tau;
							
			 } else {
				Mobility_Tnsr_tt[i][j]	=	(	Unit_diag
											+	Rij_dyadic_Rij*Rij2_inv
											+	(Unit_diag*(1.0/3.0)-(Rij_dyadic_Rij)*Rij2_inv)*(particle[i].radius*particle[i].radius+particle[j].radius*particle[j].radius)*Rij2_inv
											)	*	temp;

			 } 
			// 	Mobility_Tnsr[i][j]		=	 Mobility_Tnsr_tt[i][j] -  Mobility_Tnsr_tr[i][j]*Mobility_Tnsr_rr[i][j]*Mobility_Tnsr_rt[i][j] ;

				
			}	
	}
	
for ( int i = 0 ; i < NrParticles; i ++ )
	{
		 particle[i].frc.comp[1]+= F_g;
	}		
 		
	brownian( particle , Mobility_Tnsr_tt)	;
	
 	forceUpdate( particle, &p_energy); 	 	

	
if (step%frame==0) { 
	
     //   std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");
        
	//	outFile5<<NrParticles<<std::endl;
	//	outFile5<<"X Y Z co-ordinates"<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps
	//	for ( int i = 0 ; i < NrParticles; i ++ )
	//		{
					//	outFile5<<'H'<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
						outFile9 <<particle[0].pos.comp[0]<<'\t'<<particle[0].pos.comp[1]<<'\t'<<particle[0].pos.comp[2]<<std::endl;
						outFile10<<particle[1].pos.comp[0]<<'\t'<<particle[1].pos.comp[1]<<'\t'<<particle[1].pos.comp[2]<<std::endl;

//	}
   //   	outFile5<<'\n'<<std::endl;

	//	outFile5.close();
	}

	outFile<<K_Energy<<std::endl;
	outFile1<<p_energy<<std::endl;
	outFile3<<Temp<<'\t'<<std::endl;
	
	step+=1;
	
} while(xxnstep);
	

for (int i=0;i<NrParticles;i++) {
	outFile2<<particle[i].vel.comp[0]<<'\t'<<particle[i].vel.comp[1]<<'\t'<<particle[i].vel.comp[2]<<std::endl;
	outFile7<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;
}
outFile2.close();
outFile.close();
outFile1.close();
outFile3.close();
outFile7.close();

std::ofstream outFile8(dataFileName+"/logfile");
	outFile8<<"NrParticles"<<'\t'<<NrParticles<<std::endl;
	outFile8<<"mass"<<'\t'<<m<<std::endl;
	outFile8<<"kb"<<'\t'<<kb<<std::endl;
	outFile8<<"T0"<<'\t'<<T0<<std::endl;
	outFile8<<"box"<<'\t'<<box.comp[0]<<'\t'<<box.comp[1]<<'\t'<<box.comp[2]<<std::endl;
	outFile8<<"shear rate"<<'\t'<<shear_rate<<std::endl;
	outFile8<<"R_cut"<<'\t'<<r_cut<<std::endl;
	outFile8<<"rs"<<'\t'<<rs<<std::endl;
	outFile8<<"epsilon"<<'\t'<<epsilon<<std::endl;
	outFile8<<"sigma"<<'\t'<<sigma<<std::endl;
outFile8.close();


     // get time now
	now = time(0);
	ltm = localtime(&now);
	cout << "end time"<< '\t'<< ltm->tm_hour << ":";
	cout << ltm->tm_min << ":";
	cout << ltm->tm_sec << endl;
return 0;


}
