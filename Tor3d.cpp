/***************************************************************************
 *   Copyright (C) 2009 by Ninad Joshi   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*   This code can be used for ion and electron transport in Toroidal      *
 *   Magnetic field with constant value. The magnetic gradient can be      *
 *   added. It uses basic Point Jacobi iterative method to solve Poisson   *
 *   in circular toroidal coordinates. This file can be used as basic code *
 *   for ions in curved magnetic field and can define circular electrodes. *
 *   Any comments and suggestions are welcome. This is part of ongoing     *
 *   project call F8SR, the storage ring for ions.                         *
 *   For convenience of uploading, code has been put into single file.     *
 ***************************************************************************/
     
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>     /* srand, rand */
#include <random>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

std::mt19937 mt_gen(0);		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd() {return rnd_dist(mt_gen);}
/** time step **/
const double dT =  1.0e-10 ;
/**  Boltzmann constant **/
const double kB = 1.38e-23;
const unsigned int maxnumpar = 1000000;
const int maxnumparcell = 10000;
const int timeSteps = 31;
const double pi = 3.14;
const double EPS0 =8.854e-12;
const double R= 1.3;
const double Btor = 0.0125;
const double VLIGHT = 3.0e8;

#define Nr  11
#define Ntheta 6
#define Nzeta   11

using namespace std;
class Volume;


class vectorTor {

public: // Public attributes
  /** x-coordinate */
  double r;
  /** y-coordinate */
  double theta;
  /** z-coordinate */
  double zeta;
  /** norm of vector */
  double norm;
};

class vector3d
{
    public:
        /** coordinate */
        double x,y,z,norm;
        void norma(){ norm=sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));}
        vector3d operator+(vector3d param)
            {vector3d temp;
            temp.x=x+param.x; temp.y=y+param.y;temp.z=z+param.z;
            return(temp);}
        vector3d operator-(vector3d param)
            {vector3d temp;
            temp.x=x-param.x;temp.y=y-param.y;temp.z=z-param.z;
            return(temp);}
        vector3d operator*(double param)
            {vector3d temp;
            temp.x=x*param;temp.y=y*param;temp.z=z*param;
            return(temp);}
        double operator*(vector3d param)
            {double temp;
            temp=x*param.x+y*param.y+z*param.z;
            return(temp);}
        vector3d operator||(vector3d param)
            {vector3d temp;
            temp.x=y*param.z-z*param.y;temp.y=z*param.x-x*param.z;temp.z=x*param.y-y*param.x;
            return(temp);}
};

struct phasespace
{
    vector3d pos;
    vector3d vel;
};

vector3d Tor2Cart(vectorTor X, double R )
{

vector3d Y;
double r= X.r; double theta = X.theta; double zeta= X.zeta;

    Y.x=r*cos(theta);
    Y.y=(R + r*sin(theta))*cos(zeta);
	Y.z=(R + r*sin(theta))*sin(zeta);

return (Y);
}

vectorTor Cart2Tor( vector3d X, double R )
{
vectorTor Y;
double x= X.x;double y= X.y;double z= X.z;

		/** **************    r cordinate******************/
		Y.r=sqrt(x*x +pow(sqrt(y*y+z*z) -R,2.0));

		/** ***********    theta cordinate****************/
		Y.theta=atan2((sqrt(y*y+z*z)-R),x);
		if (Y.theta < 0.0){Y.theta = 2.0*pi+Y.theta;}
		/** **************    zeta cordinate******************/
		/** only  for +ve y */
		Y.zeta=atan2(z,y);
		/*************************/

return(Y);
}

class Particle
{
public:
vector3d pos,vel;
vectorTor pos_t;
vectorTor E;
vectorTor B;
string STATUS;
};

class Specie
{
public:
std :: vector <Particle> particles;
double mass;
double charge;
double spw;
double Energy;
void inject_ions();
void move( double );
int nLost, nFC;
};

void pos2Tor(Specie &p)
{
auto part_it = p.particles.begin();
	while(part_it != p.particles.end())
        {
        Particle &part = *part_it;
        part.pos_t= Cart2Tor( part.pos, R );
        part_it++;
        }
}

vectorTor calcB(vectorTor r)
{
vectorTor B;
    B.r=0.0;
    B.theta=0.0;
    B.zeta=Btor*R/(R+r.r*sin(r.theta));
return B;
}

vector3d torV2cart(vectorTor A, vectorTor r)
{
vector3d B;
    B.x = A.r * cos(r.theta) - A.theta*sin(r.theta);
    B.y = A.r * sin(r.theta)*cos(r.zeta) + A.theta* cos(r.theta)*cos(r.zeta) -A.zeta*sin(r.zeta);
    B.z = A.r * sin(r.theta)*sin(r.zeta) + A.theta* cos(r.theta)*sin(r.zeta) +A.zeta*cos(r.zeta);

return B;

}
void saveIons(Specie &ions)
{
std:: ofstream fout;
fout.open("beam.dat");
auto part_it = ions.particles.begin();
	while(part_it != ions.particles.end())
        {
        Particle &part = *part_it;
        //fout<< part.pos.x<<"\t"<<part.vel.x/part.vel.z<<"\t"<<part.pos.y<<"\t"<<part.vel.y/part.vel.z<<"\n";
        //fout<< part.pos.x*1e3<<"\t"<<atan2(part.vel.x,part.vel.z)*1e3<<"\t"<<(part.pos.y-R)*1e3<<"\t"<<atan2(part.vel.y,part.vel.z)*1e3<<"\n";
        fout<< part.pos.y <<"\t"<<"\t"<<part.pos.z<<"\t"<< part.pos.x<<"\n";
        part_it++;
        }
fout.close();
}



void Specie::inject_ions()
{
/** sample particles **/
    int ninject= 10;
    particles.reserve(100000);
    double rBeam = 30 *1e-3;
    double x,y,xp,yp,radius4d;
    double beamAngle=0.0*pi/180.0;
    double divBeam= 8*1e-3;
    
    double vAbs = sqrt(fabs(2.0*1.6e-19*Energy/mass));
    std :: cout<< vAbs<< " = vel\n";

    for (int i=0;i<ninject;i++)
    {
    Particle p;
    do{
        x = 2.0*rnd() -1.0;
        xp = 2.0*rnd() -1.0;
        y = 2.0*rnd() -1.0;
        yp = 2.0*rnd() -1.0;
        radius4d = sqrt(x*x+xp*xp+y*y+yp*yp);
        }while(radius4d>1.0);

        x = x*rBeam;
        y = y*rBeam;
        xp = xp *divBeam;
        yp = yp *divBeam;
		p.pos.z = 0.005;
		double xnew = x*cos(beamAngle) - xp * sin(beamAngle);
		double xpnew = xp *cos(beamAngle) + x * sin(beamAngle);

        double ynew = y*cos(beamAngle) - yp * sin(beamAngle);
		double ypnew = yp *cos(beamAngle) + y* sin(beamAngle);
		p.pos.x = xnew;
		p.pos.y = ynew;

		p.vel.x=tan(xpnew)*vAbs;
        p.vel.y=tan(ypnew)*vAbs;

        p.vel.z=sqrt(vAbs*vAbs-p.vel.x*p.vel.x - p.vel.y*p.vel.y );

        p.pos.y = p.pos.y +R;

        particles.push_back(p);
    }
}

class Volume
{
public:
	Volume();
	~Volume();

	int n_nodes;

	string*** type;
	double*** rho_i;
	double*** phi;
	double*** phiT;

double*** Er;
        double*** Etheta;
            double*** Ezeta;

	void LoadVolumeMesh();
	double dr, dtheta, dzeta;
    double rMax, zetaMax, zetaMin;
    void setBoundaries();
    void nullRho();
    void smoothRho();
    void solvePotentials();
    void solveEfields();
    void saveFields();
};
/** ********************* **/
Volume::Volume()
{
type = new string**[Nr];
rho_i = new double**[Nr];
phi = new double**[Nr];
phiT = new double**[Nr];
Er = new double**[Nr];
Etheta = new double**[Nr];
Ezeta = new double**[Nr];
for (int i = 0; i < Nr; i++)
	{
		type[i] = new string*[Ntheta];
		rho_i[i] = new double*[Ntheta];
        phi[i] = new double*[Ntheta];
        phiT[i] = new double*[Ntheta];
        Er[i] = new double*[Ntheta];
        Etheta[i] = new double*[Ntheta];
        Ezeta[i] = new double*[Ntheta];
		for (int j = 0; j < Ntheta; j++)
            {
            type[i][j] = new string[Nzeta];
            rho_i[i][j] = new double[Nzeta];
			phi[i][j] = new double[Nzeta];
			phiT[i][j] = new double[Nzeta];
			Er[i][j] = new double[Nzeta];
			Etheta[i][j] = new double[Nzeta];
			Ezeta[i][j] = new double[Nzeta];
			}
	}
}
Volume::~Volume()
{}
void Volume :: LoadVolumeMesh()
{

rMax = 0.100;
zetaMax= 30.0*pi/180.0;
zetaMin = 0.0;
dr = rMax /(Nr-1);
dtheta = 360.0*pi/180.0/Ntheta;
dzeta = zetaMax/(Nzeta-1);
n_nodes = (Nr*Ntheta*Nzeta);
std :: cout<< n_nodes<<" = n nodes\n";

}
/** ********************* **/
void Volume :: setBoundaries()
{

for (int i = 0; i < Nr; i++)
    {
   for (int j = 0; j < Ntheta; j++)
        {
		for (int k = 0; k < Nzeta; k++)
            {
            type[i][j][k] = "Normal";
            if (i ==0)
                {type[i][j][k]="Central";}
			if (i  == Nr-1)
                {type[i][j][k] ="Rbound";   }
            if ( k == 0 || k == Nzeta-1  )
                {type[i][j][k] ="Open";}     
                           
            if ( (k == 1 || k == Nzeta-2) && i ==Nr-2 )
                {type[i][j][k] ="Fixed";}     
           
            }
        }
    }

}
/** ********************* **/
void Volume :: nullRho()
{
for (int i = 0; i < Nr; i++)
    {
    for (int j = 0; j < Ntheta; j++)
        {
		for (int k = 0; k < Nzeta; k++)
            {
				rho_i[i][j][k] = 0.0;
            }}}

}
/** ********************* **/
void Scatter(Specie &p, Volume &v)
{

double charge = p.charge*p.spw;
auto part_it = p.particles.begin();
	while(part_it != p.particles.end())
        {
        Particle &part = *part_it;

        int nr =  trunc(part.pos_t.r/v.dr) ;
        double fr = (part.pos_t.r/v.dr) -nr;
        int ntheta =  trunc(part.pos_t.theta/v.dtheta) ;
        double ftheta = (part.pos_t.theta/v.dtheta)-ntheta;
        int nzeta =  trunc(part.pos_t.zeta/v.dzeta) ;
        double fzeta = (part.pos_t.zeta/v.dzeta)-nzeta;

            v.rho_i[nr][ntheta][nzeta]+= (1. - fr) * (1.0-ftheta) * (1.0-fzeta) * charge;
            v.rho_i[nr+1][ntheta][nzeta]+= (fr) * (1.0-ftheta) * (1.0-fzeta) * charge;
            v.rho_i[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta]+= (1. - fr) * (ftheta) * (1.0-fzeta) * charge;
            v.rho_i[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta] += (fr) * (ftheta) * (1.0-fzeta) * charge;
            v.rho_i[nr][ntheta][nzeta+1] += (1. - fr) * (1.0-ftheta) * (fzeta) * charge;
            v.rho_i[nr+1][ntheta][nzeta+1] += (fr) * (1.0-ftheta) * (fzeta) * charge;
            v.rho_i[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta+1] += (1. - fr) * (ftheta) * (fzeta) * charge;
            v.rho_i[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta+1] += (fr) * (ftheta) * (fzeta) * charge;

        part_it++;
        }

}
/** ********************* **/
/** ********************* **/
void Gather(Specie &p, Volume &v)
{

auto part_it = p.particles.begin();
	while(part_it != p.particles.end())
        {
        Particle &part = *part_it;

        int nr =  trunc(part.pos_t.r/v.dr) ;
        double fr = (part.pos_t.r/v.dr) -nr;
        int ntheta =  trunc(part.pos_t.theta/v.dtheta) ;
        double ftheta = (part.pos_t.theta/v.dtheta)-ntheta;
        int nzeta =  trunc(part.pos_t.zeta/v.dzeta) ;
        double fzeta = (part.pos_t.zeta/v.dzeta)-nzeta;

            part.E.r = (1. - fr) * (1.0-ftheta) * (1.0-fzeta) * v.Er[nr][ntheta][nzeta]
            + (fr) * (1.0-ftheta) * (1.0-fzeta) * v.Er[nr+1][ntheta][nzeta]
            +  (1. - fr) * (ftheta) * (1.0-fzeta) * v.Er[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (fr) * (ftheta) * (1.0-fzeta) * v.Er[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (1. - fr) * (1.0-ftheta) * (fzeta) * v.Er[nr][ntheta][nzeta+1] 
            + (fr) * (1.0-ftheta) * (fzeta) * v.Er[nr+1][ntheta][nzeta+1] 
            + (1. - fr) * (ftheta) * (fzeta) * v.Er[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta+1]
            + (fr) * (ftheta) * (fzeta) * v.Er[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta+1];

            part.E.theta = (1. - fr) * (1.0-ftheta) * (1.0-fzeta) * v.Etheta[nr][ntheta][nzeta]
            + (fr) * (1.0-ftheta) * (1.0-fzeta) * v.Etheta[nr+1][ntheta][nzeta]
            +  (1. - fr) * (ftheta) * (1.0-fzeta) * v.Etheta[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (fr) * (ftheta) * (1.0-fzeta) * v.Etheta[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (1. - fr) * (1.0-ftheta) * (fzeta) * v.Etheta[nr][ntheta][nzeta+1] 
            + (fr) * (1.0-ftheta) * (fzeta) * v.Etheta[nr+1][ntheta][nzeta+1] 
            + (1. - fr) * (ftheta) * (fzeta) * v.Etheta[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta+1]
            + (fr) * (ftheta) * (fzeta) * v.Etheta[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta+1];
            
             part.E.zeta = (1. - fr) * (1.0-ftheta) * (1.0-fzeta) * v.Ezeta[nr][ntheta][nzeta]
            + (fr) * (1.0-ftheta) * (1.0-fzeta) * v.Ezeta[nr+1][ntheta][nzeta]
            +  (1. - fr) * (ftheta) * (1.0-fzeta) * v.Ezeta[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (fr) * (ftheta) * (1.0-fzeta) * v.Ezeta[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta]
            + (1. - fr) * (1.0-ftheta) * (fzeta) * v.Ezeta[nr][ntheta][nzeta+1] 
            + (fr) * (1.0-ftheta) * (fzeta) * v.Ezeta[nr+1][ntheta][nzeta+1] 
            + (1. - fr) * (ftheta) * (fzeta) * v.Ezeta[nr][((ntheta+1)+Ntheta)%Ntheta][nzeta+1]
            + (fr) * (ftheta) * (fzeta) * v.Ezeta[nr+1][((ntheta+1)+Ntheta)%Ntheta][nzeta+1];
            
        part_it++;
        }

}
/** ********************* **/

void Volume :: smoothRho()
{
for (int i = 1; i < Nr; i++)
    {
    for (int j = 0; j < Ntheta; j++)
        {
		for (int k = 0; k < Nzeta; k++)
            {
            double r = i*dr;
            double theta = k * dtheta;
            double dV=(r+dr/2.0)*(R+r*sin(theta))* dr *dtheta*dzeta;
            dV=  pi *(2.0*r*dr +dr*dr) *dtheta *dzeta *(R+r*sin(theta));
            rho_i[i][j][k]= rho_i[i][j][k]/dV;
            }
        }
    }

for (int k = 0; k < Nzeta; k++)
    {
    int i =0;
    double rhoC=0.0;
    for (int j = 0; j < Ntheta; j++)
        {
        rhoC = rhoC + rho_i[i][j][k] ;
        }
    for (int j = 0; j < Ntheta; j++)
        {
        double r = i*dr;
        double theta = k * dtheta;
        double dV = pi *pow(dr/2.0,2.0) *R* dzeta;
	    rho_i[i][j][k] = rhoC/Ntheta/dV;
        }
     }


}

/** ********************* **/
void Volume :: solvePotentials()
{
for (int i = 0; i < Nr; i++)
    for (int j = 0; j < Ntheta; j++)
		for (int k = 0; k < Nzeta; k++)
            {
        phiT [i][j][k]=0.0;
        phi [i][j][k]=0.0;
        }

for (int it=0;it<10000; it++)
    {
    for (int i = 0; i < Nr; i++)
    	for (int j = 0; j < Ntheta; j++)
		for (int k = 0; k < Nzeta; k++)
       		     {
		     phiT [i][j][k]= phi [i][j][k];
        	     }

    for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
                double r = i*dr;
                double theta = j * dtheta;

                double RrsinT = R +r * sin(theta);

                double Aijk = 1.0/ pow(dr , 2.0) + 1.0 / (pow(r*dtheta , 2.0) )
                        + 1.0/ pow(dzeta*RrsinT, 2.0) ;

                double Arpjk = 1.0/ pow( dr , 2.0) +  (R +2.0*r * sin(theta))   / ( 2.0 *r*dr* RrsinT ) ;
                double Armjk = 1.0/ pow( dr , 2.0) -  (R +2.0*r * sin(theta))   / ( 2.0 *r*dr* RrsinT ) ;

                double Aitpk = 1.0/ pow(r *dtheta , 2.0) + cos(theta) / (2.0 *r*dtheta* RrsinT);
                double Aitmk = 1.0/ pow(r *dtheta , 2.0) - cos(theta) / (2.0 *r*dtheta* RrsinT);

                double Aijz = 1.0 / pow(dzeta*RrsinT ,2.0);

                if (type[i][j][k] == "Normal")
                    {
                    phi[i][j][k] = 1.0 /(2.0*Aijk) *
                                (
                                Arpjk * phiT[i+1][j][k]
                            +   Armjk * phiT[i-1][j][k]
                            +   Aitpk * phiT[i][((j+1)+Ntheta)%Ntheta][k]
                            +   Aitmk * phiT[i][((j-1)+Ntheta)%Ntheta][k]
                            +   Aijz * phiT[i][j][k+1]
                            +   Aijz * phiT[i][j][k-1]
                            + rho_i[i][j][k]/EPS0
                                )
                                ;
                    }
        }}}

	double phi1[Nzeta];
	        for (int k = 0; k < Nzeta; k++)
                {phi1[k] =0.0;}
        for (int k = 0; k < Nzeta; k++)
                {
	        for (int j = 0; j < Ntheta; j++)
	            {
			    phi1[k] = phi1[k]+ phiT[1][j][k];
		    }}

	for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
		if (type[i][j][k]  == "Central")
			{
			double theta = dtheta * j;
            double Ak = pi *pow(dr,2.0) /(4.0*R*dzeta) ;
			double Ajk = (Ntheta*pi*R*dzeta) + pi *pow(dr,2.0)/(2.0*R*dzeta)   ;
            double dV= pi *pow(dr/2.0,2.0) * R* dzeta;
			double Qenc = rho_i[i][j][k]*dV ;
			phi[i][j][k] =	1.0/ Ajk *
					(Qenc/EPS0
					 + phi1[k] * (pi *R* dzeta)
					 + Ak* phiT[i][j][k+1]
					 + Ak* phiT[i][j][k-1]
					)
					;
			}
			
	}}}

	for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
                if (type[i][j][k]  == "Open")
                    {
                    if (k ==0)
                        {
                        phi[i][j][k] = phiT[i][j][(Nr*Ntheta)*(k+1)%Nzeta];
                        }
                    if (k ==Nzeta-1)
                        {
                        phi[i][j][k] = phiT[i][j][(Nr*Ntheta)*(k-1)%Nzeta];
                        }
                    }

        }}}
        
    for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
                if (type[i][j][k]  == "Fixed")
                    {
                   phi[i][j][k] = -1000.0;
                    }

        }}}

     double du, dumax = 0.0;
     for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
        du =fabs(phi[i][j][k]-phiT[i][j][k]);

        if (du>=dumax) {dumax=du;}
        }}}

  //      std :: cout<<"error = "<<dumax<<"\t"<<it<<"\n";
        if (dumax<1e-5){break;}
     }

}
/** ************************************ **/

void Volume :: solveEfields()
{
    for (int i = 0; i < Nr; i++)
    	for (int j = 0; j < Ntheta; j++)
		for (int k = 0; k < Nzeta; k++)
       		     {
		     Er [i][j][k]= 0.0;
		     Etheta[i][j][k]=0.0;
		     Ezeta[i][j][k] =0.0;
        	     }
 for (int i = 0; i < Nr; i++)
        {
        for (int j = 0; j < Ntheta; j++)
            {
            for (int k = 0; k < Nzeta; k++)
                {
                double r = i*dr;
                double theta = j * dtheta;
                if (type[i][j][k] == "Normal")
                    {
                    Er[i][j][k] = - (phi[i+1][j][k] - phi[i-1][j][k])/(2.0*dr); 
                    Etheta[i][j][k] = - (phi[i][((j+1)+Ntheta)%Ntheta][k] - 
                                        phi[i][((j-1)+Ntheta)%Ntheta][k])/(2.0*r*dtheta); 
                    Ezeta[i][j][k] = - (phi[i][j][k+1] - 
                                        phi[i][j][k-1])/(2.0*dzeta)/(R+r*sin(theta));                     
                    }
                    
                if (type[i][j][k] == "Central")
                    {
                    Er[i][j][k] = 2.0* Er[i+1][j][k] - Er[i+2][j][k]; 
                    Etheta[i][j][k] = 0.0;
                    Ezeta[i][j][k] = - (phi[i][j][((k+1)+Nzeta)%Nzeta] - 
                                        phi[i][j][((k-1)+Nzeta)%Nzeta])/(2.0*dzeta)/(R+r*sin(theta));                     
                    }
                if (type[i][j][k] == "Rbound")
                    {
                    Er[i][j][k] = phi[i-1][j][k]/dr; 
                    Etheta[i][j][k] = 0.0;
                    Ezeta[i][j][k] = - (phi[i][j][((k+1)+Nzeta)%Nzeta] - 
                                        phi[i][j][((k-1)+Nzeta)%Nzeta])/(2.0*dzeta)/(R+r*sin(theta));                     
                    }

                                        
                }
            }
        }

}


/** ************************************ **/
void Volume :: saveFields()
{
std:: ofstream fout1,fout2,fout3;
vectorTor Vtor;
vector3d V;
fout1.open("rho.dat");
fout2.open("phi.dat");
fout3.open("efield.dat");
int k =1;
  for (int i = 0; i < Nr; i++)
    {
    for (int j = 0; j < Ntheta; j++)
        {
            int index =i + Nr*j + (Nr*Ntheta)* k;
            Vtor.r = i*dr;
            Vtor.theta = j * dtheta;
            Vtor.zeta = k * dzeta;
            V= Tor2Cart(Vtor, R );
            fout1 << V.x<< "\t"<<V.y<<"\t"<<rho_i[i][j][k]<<"\n";
            fout2 << V.x<< "\t"<<V.y<<"\t"<<phi[i][j][k]<<"\n";
            fout3 << V.x<< "\t"<<V.y<<"\t"<<Er[i][j][k]<<"\n";
            }
        }
fout1.close();
fout2.close();
fout3.close();

fout1.open("Bfield_long.dat");
  for (int i = 0; i < Nr; i++)
    {
    double theta = 90.0*pi/180.0;
    for (int k = 0; k < Nzeta; k++)
        {
            Vtor.r = i*dr;
           Vtor.theta = theta;
            Vtor.zeta = k * dzeta;
            V= Tor2Cart(Vtor, R );
            vectorTor B = calcB(Vtor);
            fout1 << V.z<< "\t"<<V.y<<"\t"<<B.zeta<<"\n";
            }
    }
    for (int i = 0; i < Nr; i++)
    {
    double theta = 270.0*pi/180.0;
    for (int k = 0; k < Nzeta; k++)
        {
            Vtor.r = i*dr;
              Vtor.theta = theta;
            Vtor.zeta = k * dzeta;
            V= Tor2Cart(Vtor, R );
            vectorTor B = calcB(Vtor);
            fout1 << V.z<< "\t"<<V.y<<"\t"<<B.zeta<<"\n";
            }
    }
fout1.close();


}
/** ********************* **/
phasespace BorisPush(double dt, vector3d r, vector3d v, vector3d E, vector3d B, 
                        double q, double m )
{
    double hQO2M = 0.5* dt* q /m;
    
    double gamma_t, gammaPrime, sigma, gammaNew, s ;
    vector3d rNew, vNew;

    vector3d u_t, uTpHalf, uPrime, tau, uStar, t, uNew; 

    v.norma();

    gamma_t = sqrt ( 1.0 - pow(v.norm,2.0)/ pow(VLIGHT,2.0));

    u_t  = v * gamma_t ;

    uTpHalf = u_t + (E + v || B ) * hQO2M ;

    rNew = r + uTpHalf * (dt / gamma_t) ;

    uPrime = uTpHalf + E * hQO2M;
    uPrime.norma();

    tau =  B * hQO2M;
    tau.norma();
    
    uStar = uPrime * (tau.norm/VLIGHT) ;
    
    gammaPrime = sqrt ( 1.0 + pow(uPrime.norm,2.0)/ pow(VLIGHT,2.0) ) ;
    
    sigma = pow(gammaPrime, 2.0) - pow(tau.norm, 2.0);
    
    gammaNew = sqrt ( 
                        ( 
                        sigma + sqrt(
                                     pow(sigma, 2.0) +4.0 *(
                                                         pow(tau.norm,2.0)+uStar*uStar
                                                            ) 
                                    ) 
                        ) /2.0 
                     ) ;

    t = tau *(1.0/gammaNew);
    t.norma();
    
    s = 1.0 /(1.0+pow(t.norm,2.0) ) ;
    
    uNew = ( uPrime + t * (uPrime * t) + (uPrime || t) ) * s ;
    
    vNew = uNew * (1.0 / gammaNew) ;
    
    phasespace result = {rNew,vNew};
    
    return result;
}


void Specie :: move(double dt)
{
vectorTor Bt;
vector3d E, B;
phasespace ps;
auto part_it = particles.begin();
	while(part_it != particles.end())
        {

        Particle &part = *part_it;
        
        if (part.STATUS =="Alive")
            {
            Bt = calcB(part.pos_t);
            E = torV2cart(part.E, part.pos_t);
            B = torV2cart(Bt, part.pos_t);
            ps = BorisPush(dt, part.pos, part.vel, E, B, charge, mass );
            part.pos = ps.pos ;
            part.vel = ps.vel ;       
            part_it++;
            }
         else part_it = particles.erase(part_it);	/*outside the mesh*/   
        }

}

void applyDomainConditions(Specie &p, Volume &v)
{
p.nLost =0;
p.nFC =0;
auto part_it = p.particles.begin();
	while(part_it != p.particles.end())
        {
        Particle &part = *part_it;
        if (part.pos_t.r < v.rMax && part.pos_t.zeta > v.zetaMin && part.pos_t.zeta < v.zetaMax)
           {
           part.STATUS = "Alive";
            part_it++;
           }
       if (part.pos_t.r >= v.rMax)
            {
            part.STATUS = "Lost";
            p.nLost = p.nLost+1;
            part_it = p.particles.erase(part_it);
            }
         if (part.pos_t.zeta >= v.zetaMax)
             {
             part.STATUS = "FC";
             p.nFC = p.nFC+1;
             part_it = p.particles.erase(part_it);            
             }
         if (part.pos_t.zeta <= v.zetaMin)
             {
             part.STATUS = "Lost";
             p.nLost = p.nLost+1;
             part_it = p.particles.erase(part_it);
             }                  

        }
    std :: cout << p.nLost<< " =nlost\t"<< p.nFC<< "= nFc\n";    
     
        
}

void InjectMaxwellElectrons(Specie &p, int ne, double sigma)
{
gsl_rng *rng_ptr; // pointer to random number generator (rng)
rng_ptr = gsl_rng_alloc (gsl_rng_taus);
for (int i=0; i< ne; i++ )
    {
        /*new particle*/
		Particle part;
		
		part.pos_t.r = sqrt(rnd()) * 0.080;
        part.pos_t.theta = rnd() * 2.0*pi;
        part.pos_t.zeta = rnd() * 2.0*pi / 12.0;
        part.pos = Tor2Cart(part.pos_t,R );
    part.vel.x = gsl_ran_gaussian_ziggurat(rng_ptr, sigma);
    part.vel.y = gsl_ran_gaussian_ziggurat(rng_ptr, sigma);
    part.vel.z = gsl_ran_gaussian_ziggurat(rng_ptr, sigma);
    
    p.particles.push_back(part);	
    }
gsl_rng_free(rng_ptr);
}

void savePos (Specie &p, string type, int step)
{
std : ofstream fout;
fout.open("epos_"+to_string(10000+step)+".dat");
auto part_it = p.particles.begin();
	while(part_it != p.particles.end())
        {
        Particle &part = *part_it;
        
        fout << part.pos.x << "\t"<<part.pos.y<<"\t" <<part.pos.z<<"\n";
            part_it++;
        }
        
fout.close();

}

int main()
{

    int nProcessors = omp_get_max_threads();

std :: cout << nProcessors << "Nr processors\n";

    Volume volume;

    volume.LoadVolumeMesh();
    volume.setBoundaries();

    Specie ions;
    ions.charge = 1.6e-19;
    ions.mass = 1.66e-27;
    ions.spw = 10000;
    ions.Energy = 10000.0;
    
    Specie electrons;
    electrons.charge = -1.6e-19;
    electrons.mass = 9.1e-31;
    electrons.spw = 10000;
    electrons.Energy = 5.40;
    
    double sigma = sqrt(fabs(electrons.Energy* electrons.charge/electrons.mass));
    
    InjectMaxwellElectrons(electrons, 10000, sigma);
 //   ions.inject_ions();
 //   saveIons(ions);
    for (int ts =0; ts< 1; ts++)
    {
    std:: cout << ts<< "=ts\n";
    if (ts%100 ==0)
        {
        ions.inject_ions();
        pos2Tor(ions);
        applyDomainConditions(ions, volume);
        std :: cout <<ions.particles.size()<<" = size\n"; 
        }
    
    pos2Tor(electrons);
    applyDomainConditions(electrons, volume);
    std :: cout <<electrons.particles.size()<<" = size\n"; 
    InjectMaxwellElectrons(electrons, electrons.nLost+electrons.nFC, sigma);
    
    if(ts%100 ==0)
        {
        volume.nullRho();
        Scatter(ions, volume);
        Scatter(electrons, volume);
        volume.smoothRho();
        volume.solvePotentials();
        volume.solveEfields();
        Gather(ions, volume);
        Gather(electrons, volume);
        ions.move(dT);
        }
     electrons.move(dT);
     if(ts%1000 ==0){
         savePos (electrons, "E", ts);}
    }

//saveIons(ions);
saveIons(electrons);

    volume.saveFields();

    cout << "Hello world!" << endl;
    return 0;
}

