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
 		x=((double) rand() / (RAND_MAX/Lx))-Lx/2.0;  // create particle position from -Lx/2 to Lx/2
		y=((double) rand() / (RAND_MAX/Ly))-Ly/2.0;
		z=((double) rand() / (RAND_MAX/Lz))-Lz/2.0;
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
std::normal_distribution<> R1(0,1),R2(0,1),R3(0,1);

void verlet( vector<SubData>& particle, double kbT_dt) {
	

	for(int i=0;i<NrParticles;i++) 
	{
		vctr3D rand(R1(gen), R2(gen), R3(gen));
		particle[i].vel+=(particle[i].frc-particle[i].vel*(mu))*(0.5*dt*inv_mass)+(rand*mu_sqrt*kbT_dt);
		particle[i].pos+=particle[i].vel*dt;
		particle[i].pos.PBC(box,rbox);

	}
}

void verletB(vector<SubData>& particle, double vel_scale, double kbT_dt) {
		

	if(xxthermo) 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				vctr3D rand(R1(gen), R2(gen), R3(gen));
				particle[i].vel+=(particle[i].frc-particle[i].vel*(mu))*(0.5*dt*inv_mass)+(rand*mu_sqrt*kbT_dt);
				particle[i].vel=(particle[i].vel)*vel_scale;
			}
       	} 
	else 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				vctr3D rand(R1(gen), R2(gen), R3(gen));
				particle[i].vel+=(particle[i].frc-particle[i].vel*(mu))*(0.5*dt*inv_mass)+(rand*mu_sqrt*kbT_dt);
			}
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
         
double kb=1.0 , T0=1.0, tauT=0.1;
double Temp=T0;
double shear_rate = 0.0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../xxx",dataFileName_new="../xxxnew" ;
int Max_Cluster_N=NrParticles;
double simu_time=dt;
int step=0, nSteps=10000, frame=100;
double vel_scale;
int if_Periodic =1;
double kbT_dt=sqrt(2.0*kb*Temp*dt);

std::cout<<cellx<<'\t'<<celly<<'\t'<<cellz<<std::endl;
double  K_Energy, p_energy=0;
vctr3D dR, dr2;
double R, r2;
double dr=0.05; // step size for RDF calculation
// std::vector<int> RDF((int)  floor(sqrt((Lx/2)*(Lx/2)+(Ly/2)*(Ly/2)+(Lz/2)*(Lz/2)))/dr,0), RDF1((int)  floor(sqrt(Lx*Lx+Ly*Ly))/dr,0);

vector<SubData>  particle(NrParticles);

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
        Temp+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
				   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
				   + particle[i].vel.comp[2]*particle[i].vel.comp[2]);

    }
		Temp=(Temp)/(1.5*NrParticles*kb);
		vel_scale = sqrt(T0/Temp);
		std::cout<<Temp<<'\t'<<vel_scale<<std::endl;
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
		vel_scale = sqrt(T0/Temp);
		std::cout<<Temp<<'\t'<<vel_scale<<std::endl;

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
	if (!ifrestart) {

for (int i=0;i<NrParticles;i++) {

	//	particle[i].vel={((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5),((double) rand()/(RAND_MAX)-0.5)};
					particle[i].vel=(particle[i].vel)*vel_scale;
	
	}
}

std::ofstream outFile(dataFileName+"/K_energy.dat");
std::ofstream outFile1(dataFileName+"/PE_energy.dat");
std::ofstream outFile2(dataFileName+"/Velocities.dat");
std::ofstream outFile7(dataFileName+"/End_positions.dat");
std::ofstream outFile3(dataFileName+"/Pressure_Tensor.dat");

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
	
	verlet( particle, kbT_dt)	;
	
 	forceUpdate( particle, &p_energy);
	for ( int i = 0 ; i < NrParticles; i ++ )
			{
				particle[i].frc.comp[1]-=09.8;
			}	
	K_Energy=0;
	for ( int i = 0 ; i < NrParticles; i ++ )
			{
						K_Energy+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
									   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
									   + particle[i].vel.comp[2]*particle[i].vel.comp[2]);
	}
	
	Temp=(K_Energy)/(1.5*NrParticles*kb);


if (step%frame==0) { 
	
        std::ofstream outFile5(dataFileName+"/XYZ"+ std::to_string(step/frame) +".xyz");
        
		outFile5<<NrParticles<<std::endl;
		outFile5<<"X Y Z co-ordinates"<<std::endl;
		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps
		for ( int i = 0 ; i < NrParticles; i ++ )
			{
						outFile5<<'H'<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<std::endl;

	}
      	outFile5<<'\n'<<std::endl;

		outFile5.close();
	}

	outFile<<K_Energy<<std::endl;
	outFile1<<p_energy<<std::endl;
	outFile3<<Temp<<'\t'<<std::endl;
	
	step+=1;
	vel_scale = sqrt(1.0+(T0/Temp-1.0)*(dt/tauT));
	kbT_dt=sqrt(2.0*kb*Temp*dt);
	verletB( particle , vel_scale , kbT_dt ) ;

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
