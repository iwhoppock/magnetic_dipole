#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <sys/time.h>

#define time_steps 2000000
#define vAc 3e-4   	//Alfvén velocity over speed of light
#define PI 4*atan(1)
#define dt 0.01

#define mass (1.0)
#define charge (1.0)

static inline double WTime(void){
    //timing the simulation for efficiency
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

void Boris(double x[3], double v[3], double E[3], double B[3]){
    //Boris Rotation for non-relativistic particle tracing -- "It really is state of the art" - D. Verscharen
    double v_minus[3];
    double v_prime[3];
    double v_plus[3];
    double t[3];
    double s[3];
    double t_mag2;
    
    //t vector
    t[0] = charge / mass * B[0] * 0.5 * dt;
    t[1] = charge / mass * B[1] * 0.5 * dt;
    t[2] = charge / mass * B[2] * 0.5 * dt;
    
    //|t|^2
    t_mag2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    //s vector
    s[0] = 2.0 * t[0] / (1.0 + t_mag2);
    s[1] = 2.0 * t[1] / (1.0 + t_mag2);
    s[2] = 2.0 * t[2] / (1.0 + t_mag2);

    //v minus
    v_minus[0] = v[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
    v_minus[1] = v[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
    v_minus[2] = v[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;
    
    //v prime
    v_prime[0] = v_minus[0] + ( v_minus[1] * t[2] - v_minus[2] * t[1]);
    v_prime[1] = v_minus[1] + (-v_minus[0] * t[2] + v_minus[2] * t[0]);
    v_prime[2] = v_minus[2] + ( v_minus[0] * t[1] - v_minus[1] * t[0]);
    
    //v plus:
    v_plus[0] = v_minus[0] + ( v_prime[1] * s[2] - v_prime[2] * s[1]);
    v_plus[1] = v_minus[1] + (-v_prime[0] * s[2] + v_prime[2] * s[0]);
    v_plus[2] = v_minus[2] + ( v_prime[0] * s[1] - v_prime[1] * s[0]);
    
    //final v_n+1/2
    v[0] = v_plus[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
    v[1] = v_plus[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
    v[2] = v_plus[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;

    //pusher
    x[0] += v[0] * dt;
    x[1] += v[1] * dt;
    x[2] += v[2] * dt;    
}

void dipole(double B[3], double x[3]){
    double position_mag = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    const double m = 1000.; //magnetic moment
    B[0] = (3. * m * x[0] * x[2]) / (position_mag * position_mag * position_mag * position_mag * position_mag);
    B[1] = (3. * m * x[1] * x[2]) / (position_mag * position_mag * position_mag * position_mag * position_mag);
    B[2] = ((m) / (position_mag * position_mag * position_mag)) * (( (3. * x[2] * x[2]) / (position_mag * position_mag) ) - 1.);
}


int main(){
    double E[3] = {0.,0.,0.};
    double B[3];
    double x[3] = {10.,0.,0.};
    double v[3] = {0.,0.1,0.05};
      
    //Openfile
	FILE *fout1 = NULL;
	fout1 = fopen("dipole.csv", "w");

	//TIMEING
    double time1 = WTime();

    //TIME LOOP:
    for (int nt = 0; nt < time_steps; nt++){
        dipole(B, x);
        Boris(x, v, E, B);
        if (nt % 2 == 0){
            fprintf(fout1,"%g, %g, %g, %g, %g, %g, %g\n",nt*dt,x[0],x[1],x[2],v[0],v[1],v[2]); 
        }
    }

    fclose(fout1); 

    double time2 = WTime();
    printf("TIME = %g\n",time2-time1);

    return 0;
} 