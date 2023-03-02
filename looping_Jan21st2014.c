/***
 * Langevin dynamics of polymer with crowding molecules.
 * Jaeoh Shin (jaeoh.shin@gmail.com), 2013. 11-- 
 * Importan physical quantities:
 * -Time series of the end-to-end distance of polymer chains
 * -Randomly distribute crowding molecuels, gamma_high=2.0: Nov 19th
 * -Add the beding force
 * To do list:
 * turn on or off the bending force-
 * Add sticky ends
 * periodic boundary and remove confinement
 * Jan 16th- Apply moving boundary condition for crowders
 * Jan 21st-R0 is set to 1.5.
 ***/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
//#include "/home/jshin/work/utils/mt19937-64/mt19937-64.c"
#include "/t/jshin/work/utils/mt19937-64/mt19937-64.c"

#define nple  16      // =segments number of one ring polymer*2- it should be even number
#define sple 782      // number of crowding particles
#define tple 798      //=nple+sple
#define config_capture 0 // store equilibration configuration 1 or not
#define one_chain 1 
#define ring_ 0
#define pi 3.141592
//#define pbc 1 //periodic boundary condition
int confinement = 1 ;
int pbc=0;

int roundToint(double x);

// Control the simulation
//#define confinement 1 //1 for confinement 0 for free space
#define bendingpotential 0 //
#define stickyends 1


int init_steps = 2e5 ;//+nple*1e3 ;     //;
int eq_steps = 1e5+sple*5e2 ;
int steps= 2e7+1e5+sple*5e2 ;   //4e6+eq_steps;
int capture= 1e2; //~2e3

/**
 parameters for external potential
**/
double rcy, rcyt=8.5, rcy0, Len ; //cylinder radius and length
double k_z, k_zt=100. ;//the spring constant of external potential along the z-axis spring constant, i.e. ext_u=0.5*k_z*dx^2
double kappa=0.0;
/* parameters for LJ and FENE force **/
double rcut,rcut2; //cutoff distance of LJ potential and its square
double k_fene;     //FENE spring constant
double R0, R02;    //Maximally allowed extension of FENE spring
double eps=1.0, eps_sticky=5.0;

double boxt=16.0;
double box0, box;


/**
* neighbour list and related arrays: Reference, M.P. Allen and D. J. Tildesley (1987)
**/
short int nlist[tple][tple];
double displacement[tple]; //displacement after last nlist update
double rskin, rskin2;
int nlist_updates=0;


/**
 * Update list
 **/
int update_nlist(double *x, double *y, double *z);


/**
 * Check whether updating neighbor list is needed or not
 * */
int check_nlist() ;


/**
 * Force: Lenneard-Jones, FENE, and External- infinte surface or patch
 **/
 
int lj_force(double *x, double *y, double *z, double *fljx, double *fljy, double *fljz);

int fene_force(double *x, double *y, double *z, double *fene_x, double *fene_y, double *fene_z);

int ext_confine(double *x, double *y, double *z, double *extx, double *exty, double *extz);

int bending_force(double *x, double *y, double *z, double *fx_angle, double *fy_angle, double *fz_angle);


/** radius of gyration **/
double rg(double *x, double *y, double *z);
/** end to end distance **/
double ete(double *x, double *y, double *z);


/** Makes two Gaussian random numbers with zero mean and unit variance.
 * Argument is double array with two elements. **/
void gaussrand(double* grnd);

/** 
 * Generation of random number seed, get a current time in micri-second unit
 * */
int gus() {
	struct timeval tp;
	gettimeofday(&tp,0);
	return tp.tv_usec;
}

int main(){
	
	double dt ;
	unsigned long long idum;
	idum=gus();
	init_genrand64(idum);
	time_t start, end ;
	time(&start);
	int simulation_time ;
	int abc ;
	int temp_integer1, temp_integer2, temp_integer3, temp_integer4 ;
	int Nint ;
        int is_overlapped ;	
	double R1, N0;
	double R, xcm, ycm, zcm;
	double xpcm, ypcm, zpcm ;
	double sigma0 ;
	int i, j, k, kk,  ctime, yy;
	FILE *infoFile;   //store general simulation info
	FILE *eqFile;     //equilibrium configuration
	FILE *config;
	FILE *endFile;
	FILE *rgFile;
	  config=fopen("./ini_config.txt","w") ;
	  eqFile=fopen("./eqFile.txt","w");
	  endFile=fopen("./eteFile.txt", "w");
 	  rgFile=fopen("./rg.txt", "w");
	rcut=1.12246205;
	rcut2=rcut*rcut;
	sigma0=1.3 ;

	R0=1.5 ;
	R02=R0*R0;
	k_fene=30.0 ;
	
	rskin=4.0; //Verlet-list skin radius. ~4.0 should be fine. 
	rskin2=rskin*rskin;

	double grnd[2]; // gaussian random numbers
	grnd[0]=0.0;
	grnd[1]=0.0;
	
/**
 system parameters
**/
	double T=1.0,m=1.0; //temperature, mass
	double gamma, gamma_low=1.0, gamma_high=2.0;
	double lj_U;
	
	int need_update_nlist=0;
		
	double x[tple],y[tple],z[tple];
	double dx, dy, dz;
	double vx[tple],vy[tple],vz[tple];
	
	
	double fx[tple], fy[tple], fz[tple] ;
	double fx2[tple], fy2[tple], fz2[tple];
	double fljx[tple],fljy[tple],fljz[tple]; 
	double fene_x[tple], fene_y[tple], fene_z[tple];
	double extx[tple], exty[tple], extz[tple];
	double fx_angle[tple], fy_angle[tple], fz_angle[tple];
	
	
	double drx[tple], dry[tple], drz[tple];
	double drvx[tple], drvy[tple], drvz[tple];
	
	double vmx, vmy, vmz,temp1,temp2,temp3, temp4,temp5 ;
	double fs;
	double radii, ladii;


/** Coefficients for the Li algorithm. The algorithm can be written as:
     x[i+1]=x[i]+cx1*vx[i]+cx2*force(x[i])+drx[i];
     vx[i+1]=cv1*vx[i]+cv2*force(x[i])+cv3*force(x[i+1])+drvx[i];
     The coefficients are constants, which are defined as: (see, e.g., Allen-Tildesley Sec. 9.3)**/
	
/** The random displacements drx and drvx are correlated Gaussian random numbers.
     They can be obtained from zero-mean unit-variance distribution by transformation
     drx[i]=cdr1*grnd[0]+cdr2*grnd[1];
     drvx[i]=cdrv*grnd[0];
     given in, e.g., Allen-Tildesley Sec 9.3 and Appendix G.3.**/

	double c0, c1, c2, cx1, cx2, cv1, cv2, cv3 ;
	double varr, varv, covrv;
	double cdr1, cdr2, cdrv;


/**********       Loop for different chain extension     **********/

  int n_sumsum=0 ;
	yy=0 ;
/** Initialize the position of segment: for ring polymer-straight line then turn back **/
// linear polymer- straight 
// x[0]=0.0;
// y[0]=0.0;
// z[0]=0.0;


 //Initialize the polymer chain
 rcy0=nple*sigma0;
 rcy=rcy0;
 box0=nple*sigma0+4.0;
 box=box0;
  for(j=0;j<nple;j++){
	 x[j]=0.5*sin(0.2*genrand64_real2());
	 y[j]=0.5*cos(0.2*genrand64_real2());
	 z[j]=-rcy+j*sigma0; }
  
  temp3=0.0;
  
  for(j=0;j<nple;j++){
	  temp3=temp3+z[j] ;
	  }
	temp3=temp3/(nple);
  for(j=0;j<nple;j++){
	  z[j]=z[j]-temp3 ;
	  }
	  
	  

/**  
 * Initialize crowding molecules in random position inside the shell
 **/
   box=box0; 
  for(j=nple;j<tple;j++){
    x[j]=(rcy-2.0)*(genrand64_real2()-0.5); //randomply choose the position, it should be inside a sphere!!
    y[j]=(rcy-2.0)*(genrand64_real2()-0.5);
    z[j]=(rcy-2.0)*(genrand64_real2()-0.5);
   is_overlapped=0;

    for (k=0;k<j;k++){
      temp1=(x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k])+(z[j]-z[k])*(z[j]-z[k]) ;
      if( (temp1<1.1*1.1) ) {is_overlapped=1;
	     break;} 
    }
//   if(  x[j]*x[j]+y[j]*y[j]+z[j]*z[j] > (rcy-1.0)*(rcy-1.0) ) {is_overlapped=1; }
if (is_overlapped) {
	j--;
	continue;}
	
  }
  
 for(i=0;i<tple;i++){
   fprintf(config, "%f\t%f\t%f\n", x[i],y[i],z[i]);}
   fclose(config);

/** initialize the velocity **/
	vmx=0., vmy=0.,vmz=0.;
	for(j=0;j<tple;j++){
		vx[j]=(genrand64_real2()-0.5);
		vy[j]=(genrand64_real2()-0.5);
		vz[j]=(genrand64_real2()-0.5);
		vmx=vmx+vx[j];
		vmy=vmy+vy[j];
		vmz=vmz+vz[j];
	}
	
	vmx=vmx/(double)tple;
	vmy=vmy/(double)tple;
	vmz=vmz/(double)tple;
	temp1=0.;
	
	for(j=0;j<tple;j++){
		vx[j]=vx[j]-vmx;
		vy[j]=vy[j]-vmy;
		vz[j]=vz[j]-vmz;
		temp1=temp1 +vx[j]*vx[j] +vy[j]*vy[j] +vz[j]*vz[j];	}

	temp1=temp1/(3.0*(double)tple);
	fs=sqrt(T/temp1);
	temp1=0.;

	for(j=0;j<tple;j++){
		vx[j]=fs*vx[j];
		vy[j]=fs*vy[j];
		vz[j]=fs*vz[j];
		temp1=temp1+vx[j]*vx[j]+vy[j]*vy[j]+vz[j]*vz[j];
	}
	
	// setting the coefficients of Li algorithm for equilibration
	    dt= 0.002;
	    gamma=gamma_low;

		c0=exp(-1.0*gamma*dt);
		c1=(1.0-c0)/(gamma*dt);
		c2=(1.0-c1)/(gamma*dt);
		
		cx1=c1*dt;
		cx2=c2*dt*dt/m;
		
		cv1=c0;
		cv2=(c1-c2)*dt/m;
		cv3=c2*dt/m;
		
		varr=T*dt/(gamma*m)*(2.0-(3.0-4.0*exp(-1.0*gamma*dt)+exp(-2.0*gamma*dt))/(gamma*dt));
		varv=T/m*(1.0-exp(-2.0*gamma*dt));
		covrv=T/(gamma*m)*(1.0-exp(-1.0*gamma*dt))*(1.0-exp(-1.0*gamma*dt))/(sqrt(varr)*sqrt(varv));
		
		cdr1=sqrt(varr)*covrv;
		cdr2=sqrt(varr)*sqrt(1.0-covrv*covrv);
		cdrv=sqrt(varv);


	/** integration loop **/	
	update_nlist(x,y,z);
	lj_force(x, y, z, fljx, fljy, fljz) ;
	fene_force(x, y, z, fene_x, fene_y, fene_z);
	ext_confine(x, y, z, extx, exty,extz);
	bending_force(x, y, z, fx_angle, fy_angle, fz_angle);

     for(i=0;i<=init_steps;i++){
	rcy=rcy0-(rcy0-rcyt)*(double)i/init_steps ;
	box= box0-(box0-boxt)*(double)i/init_steps ;//The cylinder size is continously decreasing
	for(j=0; j<tple;j++){
			gaussrand(grnd);
			drx[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvx[j]=cdrv*grnd[0];
			
			gaussrand(grnd);
			dry[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvy[j]=cdrv*grnd[0];
			
			gaussrand(grnd);
			drz[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvz[j]=cdrv*grnd[0];			
		}
		
		// New value for x and 
		for(j=0;j<tple;j++){
			fx[j]=fljx[j]+fene_x[j]+extx[j]+fx_angle[j];
			fy[j]=fljy[j]+fene_y[j]+exty[j]+fy_angle[j];
			fz[j]=fljz[j]+fene_z[j]+extz[j]+fz_angle[j];
			
			dx=cx1*vx[j]+cx2*fx[j]+drx[j];
			x[j]=x[j]+dx;
			
			dy=cx1*vy[j]+cx2*fy[j]+dry[j];
			y[j]=y[j]+dy;
			
			dz=cx1*vz[j]+cx2*fz[j]+drz[j];
			z[j]=z[j]+dz; 
			displacement[j]=displacement[j]+sqrt(dx*dx+dy*dy+dz*dz);
						}
						
						

//update neighbor list if the displacement is larger than rskin-rcut
			need_update_nlist=check_nlist();
			if(need_update_nlist) {
				update_nlist(x,y,z);
			}
			
	lj_force(x, y, z, fljx, fljy, fljz) ;
	fene_force(x, y, z, fene_x, fene_y, fene_z); 
	ext_confine(x, y, z, extx, exty,extz);
	bending_force(x, y, z, fx_angle, fy_angle, fz_angle);
	
	for(j=0;j<tple;j++){
			fx2[j]=fljx[j]+fene_x[j]+extx[j]+fx_angle[j];
			fy2[j]=fljy[j]+fene_y[j]+exty[j]+fy_angle[j];
			fz2[j]=fljz[j]+fene_z[j]+extz[j]+fz_angle[j];

			vx[j]=cv1*vx[j]+cv2*fx[j]+cv3*fx2[j]+drvx[j];
			vy[j]=cv1*vy[j]+cv2*fy[j]+cv3*fy2[j]+drvy[j];
			vz[j]=cv1*vz[j]+cv2*fz[j]+cv3*fz2[j]+drvz[j];
		}
	}
		    /** initial steps end **/
	

/** 2014 January 16th **/
/** No confinement. Non periodic b.c. for polymer and PBC for crowders   **/
	confinement = 0; 
	pbc=1;
	
	

// setting the coefficients of Li algorithm	
		dt= 0.01 ;
		gamma=gamma_high;
		c0=exp(-1.0*gamma*dt);
		c1=(1.0-c0)/(gamma*dt);
		c2=(1.0-c1)/(gamma*dt);
		
		cx1=c1*dt;
		cx2=c2*dt*dt/m;
		
		cv1=c0;
		cv2=(c1-c2)*dt/m;
		cv3=c2*dt/m;
		
		varr=T*dt/(gamma*m)*(2.0-(3.0-4.0*exp(-1.0*gamma*dt)+exp(-2.0*gamma*dt))/(gamma*dt));
		varv=T/m*(1.0-exp(-2.0*gamma*dt));
		covrv=T/(gamma*m)*(1.0-exp(-1.0*gamma*dt))*(1.0-exp(-1.0*gamma*dt))/(sqrt(varr)*sqrt(varv));
		
		cdr1=sqrt(varr)*covrv;
		cdr2=sqrt(varr)*sqrt(1.0-covrv*covrv);
		cdrv=sqrt(varv);

	update_nlist(x,y,z);
	lj_force(x, y, z, fljx, fljy, fljz) ;
	fene_force(x, y, z, fene_x, fene_y, fene_z); 
	ext_confine(x, y, z, extx, exty,extz);
	bending_force(x, y, z, fx_angle, fy_angle, fz_angle);
	ctime=0;

				/***** integration loop  *****/
  for(i=0;i<=steps;i++){
    
    // Jan 16th 2014-measure cm of polymer and apply periodic bc to crowders
	xpcm=0.0;
	ypcm=0.0;
	zpcm=0.0;
		for(j=0;j<nple;j++){
		  xpcm=xpcm+x[j];
		  ypcm=ypcm+y[j];
		  zpcm=zpcm+z[j];
		}
		
		xpcm=xpcm/(double)nple;
		ypcm=ypcm/(double)nple;
		zpcm=zpcm/(double)nple;
		
		for(j=nple;j<tple;j++){
		  if(pbc){
		   x[j]=x[j]-roundToint((x[j]-xpcm)/box)*(double)box;
		   y[j]=y[j]-roundToint((y[j]-ypcm)/box)*(double)box;
		   z[j]=z[j]-roundToint((z[j]-zpcm)/box)*(double)box;}
				      }

		for(j=0; j<tple;j++){
			gaussrand(grnd);
			drx[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvx[j]=cdrv*grnd[0];
			
			gaussrand(grnd);
			dry[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvy[j]=cdrv*grnd[0];
			
			gaussrand(grnd);
			drz[j]=cdr1*grnd[0]+cdr2*grnd[1];
			drvz[j]=cdrv*grnd[0];			
		}
		
		for(j=0;j<tple;j++){
			fx[j]=fljx[j]+fene_x[j]+extx[j]+fx_angle[j];
			fy[j]=fljy[j]+fene_y[j]+exty[j]+fy_angle[j];
			fz[j]=fljz[j]+fene_z[j]+extz[j]+fz_angle[j];
			
			dx=cx1*vx[j]+cx2*fx[j]+drx[j];
			x[j]=x[j]+dx;
			
			dy=cx1*vy[j]+cx2*fy[j]+dry[j];
			y[j]=y[j]+dy;
			
			dz=cx1*vz[j]+cx2*fz[j]+drz[j];
			z[j]=z[j]+dz; 
			
			displacement[j]=displacement[j]+sqrt(dx*dx+dy*dy+dz*dz);
			 }


			 
	//update neighbor list if displacement is larger than (rskin-rcut)
			need_update_nlist=check_nlist();
			if(need_update_nlist) {
				update_nlist(x,y,z);
			}
			
	lj_force(x, y, z, fljx, fljy, fljz) ;
	fene_force(x, y, z, fene_x, fene_y, fene_z);
	ext_confine(x, y, z, extx, exty,extz);
	bending_force(x, y, z, fx_angle, fy_angle, fz_angle);
	
		for(j=0;j<tple;j++){
			fx2[j]=fljx[j]+fene_x[j]+extx[j]+fx_angle[j];
			fy2[j]=fljy[j]+fene_y[j]+exty[j]+fy_angle[j];
			fz2[j]=fljz[j]+fene_z[j]+extz[j]+fz_angle[j];
			
			vx[j]=cv1*vx[j]+cv2*fx[j]+cv3*fx2[j]+drvx[j];
			vy[j]=cv1*vy[j]+cv2*fy[j]+cv3*fy2[j]+drvy[j];
			vz[j]=cv1*vz[j]+cv2*fz[j]+cv3*fz2[j]+drvz[j];
		}

	


/*** 
 storing data after some equilibration time
***/
	if( (i>eq_steps)&&((i%capture)==0)) {
		ctime++;
		int n_sum = 0 ;
		yy=yy+1 ;
	temp4=sqrt( (x[0]-x[nple-1])*(x[0]-x[nple-1])+(y[0]-y[nple-1])*(y[0]-y[nple-1])+(z[0]-z[nple-1])*(z[0]-z[nple-1]));
	temp5=rg(x,y,z);
		fprintf(endFile, "%d\t%4.3f\n", (int)(i*dt), temp4);
		fprintf(rgFile,  "%d\t%4.3f\n", (int)(i*dt), temp5);

		}
		
/**
 *  Storing chain coordinate and velocity 
**/
	if ( (i>eq_steps+capture)&&(config_capture==1) && (i%1000)==0) {
			fprintf(eqFile, "ITEM: TIMESTEP\n"); //Output format for VMD.
			fprintf(eqFile,"%d\n", (int)i);
			fprintf(eqFile,"ITEM: NUMBER OF ATOMS\n");
			fprintf(eqFile,"%d\n", (int)tple);
			fprintf(eqFile,"ITEM: BOX BOUNDS pp pp pp\n");
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile,"%d %d\n", 0, (int)box);
			fprintf(eqFile, "ITEM: ATOMS id type xu yu zu vx vy vz\n");
			for(j=0;j<tple;j++){
				if(j<nple/2)
					fprintf(eqFile,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f\n", (int)(j+1), 1, x[j], y[j], z[j], vx[j], vy[j], vz[j]);
				else if(j<nple)
					fprintf(eqFile,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f\n", (int)(j+1), 2, x[j], y[j], z[j], vx[j], vy[j], vz[j]);
				else	
				        fprintf(eqFile,"%d %d %.3f %.3f %.3f %.3f %.3f %.3f\n", (int)(j+1), 3, x[j], y[j], z[j], vx[j], vy[j], vz[j]);
								}
													}
													
		}
		//printf("End of equilibration\n");
   				/** Steps loop end. **/


	infoFile = fopen("./info.txt", "a");
 	time(&end);
 	simulation_time = (int)difftime(end, start);
	fprintf(infoFile, "Simulation time= %d minutes %d sec\n", (int)(simulation_time/60),(simulation_time%60));
	fclose(infoFile);
	fclose(eqFile);
	fclose(endFile);
	printf("End of simulation!\n\n");
 	return 0;
}
				      /** End of main ******/




/**
 * Functions are defined below.
 * */


// FENE force
int fene_force(double *x, double *y, double *z, double *fene_x, double *fene_y, double *fene_z){
	int j;
	int aa, bb; //the first and the last segment
	double xdist,ydist,zdist, dist2, fene_factor;
	
	for(j=0;j<tple;j++){
		fene_x[j]=0.0;
		fene_y[j]=0.0;
		fene_z[j]=0.0;
	}

    //FENE force for 1st polymer
    
    if (one_chain){aa=0, bb=nple-1 ;}
    else { aa=0, bb=nple/2-1; }			// the first and last segment of 1st ring polymer
    for(j=aa;j<bb;j++){
		xdist=x[j+1]-x[j];
		ydist=y[j+1]-y[j];
		zdist=z[j+1]-z[j];
		dist2=xdist*xdist+ ydist*ydist+ zdist*zdist;
		fene_factor=k_fene*R02/(R02-dist2);
		fene_x[j]=fene_x[j]+fene_factor*xdist;
		fene_y[j]=fene_y[j]+fene_factor*ydist;
		fene_z[j]=fene_z[j]+fene_factor*zdist;
		
		fene_x[j+1]=fene_x[j+1]-fene_factor*xdist;
		fene_y[j+1]=fene_y[j+1]-fene_factor*ydist;
		fene_z[j+1]=fene_z[j+1]-fene_factor*zdist;
		}
  /** For ring polymer add additional force between two ends segments **/
	if(ring_){xdist=x[bb]-x[aa] ;
		  ydist=y[bb]-y[aa] ;
		  zdist=z[bb]-z[aa] ;
		  dist2=xdist*xdist+ ydist*ydist+ zdist*zdist;
		  
		  fene_factor=k_fene*R02/(R02-dist2);
		  fene_x[aa]=fene_x[aa]+fene_factor*xdist;
		  fene_y[aa]=fene_y[aa]+fene_factor*ydist;
		  fene_z[aa]=fene_z[aa]+fene_factor*zdist;
		  
		  fene_x[bb]=fene_x[bb]-fene_factor*xdist;
		  fene_y[bb]=fene_y[bb]-fene_factor*ydist;
		  fene_z[bb]=fene_z[bb]-fene_factor*zdist;
	}
	return 0;
}



int lj_force(double *x, double *y, double *z, double *fljx, double *fljy, double *fljz){
	int j, k;
	double xdist, ydist, zdist,dist2,dist2i,dist6i, lj_factor, xabs,yabs,zabs,hbox ;
	for(j=0;j<tple;j++){
		fljx[j]=0.0;
		fljy[j]=0.0;
		fljz[j]=0.0; }

	for(j=0;j<tple;j++){
		for(k=j+1 ;k<tple ; k++){
			
		if(nlist[j][k]==1) {   //check the neighbor list
			xdist=x[j]-x[k];
			ydist=y[j]-y[k];
			zdist=z[j]-z[k];
		// Interaction between closed images for crowders
			if(pbc){
				if(j>nple-1){
					xdist=xdist-roundToint(xdist/box)*(double)box;
					ydist=ydist-roundToint(ydist/box)*(double)box;
					zdist=zdist-roundToint(zdist/box)*(double)box;	}
				}

			dist2=xdist*xdist +ydist*ydist +zdist*zdist;
			
			if((j==0 && k==(nple-1)) && stickyends ) {
				if(dist2<9.0){
				dist2i=1./dist2;
				dist6i=dist2i*dist2i*dist2i;
				lj_factor=48.0*eps_sticky*dist2i*dist6i*(dist6i-0.5);
				fljx[j]=fljx[j]+lj_factor*xdist;
				fljy[j]=fljy[j]+lj_factor*ydist;
				fljz[j]=fljz[j]+lj_factor*zdist;
				fljx[k]=fljx[k]-lj_factor*xdist;
				fljy[k]=fljy[k]-lj_factor*ydist;
				fljz[k]=fljz[k]-lj_factor*zdist;
							}
						}
			else {
				if(dist2<rcut2){
					dist2i=1./dist2;
					dist6i=dist2i*dist2i*dist2i;
					lj_factor=48.0*eps*dist2i*dist6i*(dist6i-0.5);
				
					fljx[j]=fljx[j]+lj_factor*xdist;
					fljy[j]=fljy[j]+lj_factor*ydist;
					fljz[j]=fljz[j]+lj_factor*zdist;
					fljx[k]=fljx[k]-lj_factor*xdist;
					fljy[k]=fljy[k]-lj_factor*ydist;
					fljz[k]=fljz[k]-lj_factor*zdist;
						}
				}
			}
		}
	}
	return 0;
}
 
 
/**
 * Confinement force
 **/
int ext_confine(double *x, double *y, double *z, double *extx, double *exty, double *extz){
	int j;
	double f, r;
	double dist, dist2, dist2i, dist6i, lj_factor;
	for(j=0;j<tple;j++){
		extx[j]=0.;
		exty[j]=0.;
		extz[j]=0.;
	}
	
 if(confinement==1) {
		
    for(j=0;j<tple;j++){
	r=rcy-sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]) ;  //the separation distance from sphere wall
	if(r < rcut) {
	  dist2i=1./(r*r) ;
	  dist6i=dist2i*dist2i*dist2i ;
	  lj_factor=48.0*eps*dist2i*dist6i*(dist6i-0.5);
	  extx[j]=-lj_factor*x[j];
	  exty[j]=-lj_factor*y[j];
	  extz[j]=-lj_factor*z[j];}
    }
  }
	return 0;
}


/**
 * Bending force: ref: A Brownian Dynamis program for the simulatioin of linear and circular DNA and other wormlike chain polyelectrolytes
**/
int bending_force(double *x, double *y, double *z, double *fx_angle, double *fy_angle, double *fz_angle)
	{
  int j ;
  double fx_next[nple], fy_next[nple], fz_next[nple] ;
  double fx_prev[nple], fy_prev[nple], fz_prev[nple] ;
  double ux[nple-1], uy[nple-1], uz[nple-1], u[nple-1], usquare  ;
  double ex[nple-1], ey[nple-1], ez[nple-1] ;
  double costheta, theta[nple-1] ;
  double px[nple-1 ], py[nple-1], pz[nple-1], p[nple-1] ; 
  double psqrt  ;

 for(j=0;j<tple;j++){
		fx_angle[j]=0.0;
		fy_angle[j]=0.0;
		fz_angle[j]=0.0;
	}
if(bendingpotential==1) {
  for(j=0 ; j<= nple-1; j++)   {
	fx_angle[j] = 0.0, fy_angle[j] = 0.0, fz_angle[j]= 0.0  ;
	fx_next[j] = 0.0,  fy_next[j] = 0.0,  fz_next[j] = 0.0 ;
	fx_prev[j] = 0.0,  fy_prev[j] = 0.0,  fz_prev[j] = 0.0 ; }

  for(j=0; j<= nple-2; j++) { 
	ux[j] = 0.0, uy[j] = 0.0, uz[j] = 0.0, u[j] = 0.0 ;
	theta[j] = 0.0, px[j]=0.0, py[j]=0.0, pz[j]= 0.0 ;  }

  for(j=0; j<= nple-2; j++) {         // tangential vector
	ux[j] = x[j+1] - x[j];        
	uy[j] = y[j+1] - y[j];
	uz[j] = z[j+1] - z[j];}

  for(j=0; j<= nple-2; j++) {         // tangential vector
	u[j]= ux[j]*ux[j] + uy[j]*uy[j] + uz[j]*uz[j] ;
	u[j]= sqrt(u[j] );

	ex[j] = ux[j]/u[j] ;			//unit tangential vector 
	ey[j] = uy[j]/u[j] ;
	ez[j] = uz[j]/u[j] ; }

  for(j=1; j<=nple-2; j++)  {
	costheta = ex[j]*ex[j-1]+ ey[j]*ey[j-1]+ ez[j]*ez[j-1] ;
	theta [j] =  acos(costheta);}

  for(j=1; j<=nple-2; j++) {
	px[j] = ey[j-1]*ez[j] - ez[j-1]*ey[j] ;
	py[j] = ez[j-1]*ex[j] - ex[j-1]*ez[j] ;
	pz[j] = ex[j-1]*ey[j] - ey[j-1]*ex[j] ; 
	psqrt = px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j] ;
	psqrt =  1.0/sqrt(psqrt);
	px[j] = px[j]*psqrt;  // unit vector
	py[j] = py[j]*psqrt;
	pz[j] = pz[j]*psqrt;}

  for(j=1; j<=nple-2; j++) {
	fx_next[j] = py[j]*ez[j]-pz[j]*ey[j];
	fy_next[j] = pz[j]*ex[j]-px[j]*ez[j];
	fz_next[j] = px[j]*ey[j]-py[j]*ex[j];

	fx_next[j] = kappa*theta[j]*fx_next[j]/u[j];
	fy_next[j] = kappa*theta[j]*fy_next[j]/u[j];
	fz_next[j] = kappa*theta[j]*fz_next[j]/u[j];}

  for(j=1; j<=nple-2 ; j++) {
	fx_prev[j] = py[j]*ez[j-1] - pz[j]*ey[j-1];
	fy_prev[j] = pz[j]*ex[j-1] - px[j]*ez[j-1];
	fz_prev[j] = px[j]*ey[j-1] - py[j]*ex[j-1];

	fx_prev[j] = kappa*theta[j]*fx_prev[j]/u[j-1] ;
	fy_prev[j] = kappa*theta[j]*fy_prev[j]/u[j-1] ;
	fz_prev[j] = kappa*theta[j]*fz_prev[j]/u[j-1] ; }

	/******************************************/
	j=0 ;
	fx_angle[j] =  -fx_prev[j+1];
	fy_angle[j] =  -fy_prev[j+1];
	fz_angle[j] =  -fz_prev[j+1];

	j=1 ;
	fx_angle[j] = fx_prev[j] + fx_next[j] - fx_prev[j+1];
	fy_angle[j] = fy_prev[j] + fy_next[j] - fy_prev[j+1];
	fz_angle[j] = fz_prev[j] + fz_next[j] - fz_prev[j+1];

	j=nple-2 ;
	fx_angle[j] = -fx_next[j-1] +  fx_prev[j] + fx_next[j];
	fy_angle[j] = -fy_next[j-1] +  fy_prev[j] + fy_next[j];
	fz_angle[j] = -fz_next[j-1] +  fz_prev[j] + fz_next[j];

	j=nple-1 ;
	fx_angle[j] = -fx_next[j-1] ;
	fy_angle[j] = -fy_next[j-1] ;
	fz_angle[j] = -fz_next[j-1] ;

  for(j=2; j<=nple-3; j++){
	fx_angle[j] = -fx_next[j-1] +  fx_prev[j] + fx_next[j] - fx_prev[j+1];
	fy_angle[j] = -fy_next[j-1] +  fy_prev[j] + fy_next[j] - fy_prev[j+1];
	fz_angle[j] = -fz_next[j-1] +  fz_prev[j] + fz_next[j] - fz_prev[j+1];
			 }
 }
  
  return 0;
}



//Calculate radius of gyration^2
double rg(double *x, double *y, double *z) {
	int j;
	double xcm, ycm, zcm,inv,rg;
	xcm=0.0 ;ycm=0.0 ;zcm=0.0 ;
	for(j=0;j<nple;j++){
	      xcm=xcm+x[j];
	      ycm=ycm+y[j];
	      zcm=zcm+z[j];
	}
	inv=1.0/(double)nple;
	xcm=xcm*inv;
	ycm=ycm*inv;
	zcm=zcm*inv;
	rg=0.0;
	for(j=0;j<nple;j++){
		rg=rg+(x[j]-xcm)*(x[j]-xcm)+(y[j]-ycm)*(y[j]-ycm)+(z[j]-zcm)*(z[j]-zcm);
	}
	rg=rg*inv;
	return rg;
}

//
//double ete(double *x, double *y, double *z){
//}




/** 
 * Neighbor list
 **/
int update_nlist(double *x, double *y, double *z){
	int j,k;
	double xdist, ydist, zdist, dist2;
	
	for(j=0;j<tple-1;j++){
		for(k=j+1;k<tple; k++){
			xdist=x[j]-x[k];
			ydist=y[j]-y[k];
			zdist=z[j]-z[k];
						// Interaction between closed images for crowders
			if(pbc){
				if(j>nple-1){
					xdist=xdist-roundToint(xdist/box)*(double)box;
					ydist=ydist-roundToint(ydist/box)*(double)box;
					zdist=zdist-roundToint(zdist/box)*(double)box;	}
				}
				
			dist2=xdist*xdist +ydist*ydist +zdist*zdist;
			if(dist2<rskin2){
				nlist[j][k]=1;
				nlist[k][j]=1;
			}else{
				nlist[j][k]=0;
				nlist[k][j]=0;
			}
		}
	}
	for(j=0;j<tple;j++){
		displacement[j]=0.0;
	}
	nlist_updates++;
	return 0;
}



	
	
	

/**
 * Check neighbor list
 * */
int check_nlist(){
	int j;
	
	double maxdisp=0.0;
	double max2disp=0.0;
	
	for(j=0;j<tple;j++){ //Give two largest displacement
		if(displacement[j]>maxdisp){
			max2disp=maxdisp;
			maxdisp=displacement[j];
		}else if(displacement[j]>max2disp){
			max2disp=displacement[j];
		}else {
			//do nothing
		}
	}
	if(maxdisp+max2disp>rskin-rcut){
		return 1;
	}else {
		return 0;
	}
}

/**
 * Gives two gaussian random numbers with zero mean and unit variance.
 * Argument is double array with two elements.
 * Uses the polar form of Box-Muller transform to convert two random
 * variables of [0,1)-uniform distribution into gaussian random numbers.
 */
void gaussrand(double* grnd){
	double r1, r2,y1,y2, w;
		do{
			r1=2.0*genrand64_real2()-1.0;
			r2=2.0*genrand64_real2()-1.0;
			w=r1*r1+r2*r2;
		}while(w >= 1.0);

	w=sqrt((-2.0*log(w))/w);
	y1=r1*w;
	y2=r2*w;
	grnd[0]=y1;
	grnd[1]=y2;
	}
	
	// Define roundToint
int roundToint(double x){
	if(x>=0) return (int)(x+0.5);
	return (int)(x-0.5);
}
