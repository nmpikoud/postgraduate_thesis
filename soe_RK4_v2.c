/* This program implements the classical fourth-order Runge-Kutta method for
a system of ordinary differential equations using a constant step.
The user must do three things.
First, he must define the multide of the equations as a directive.
Second, at the first lines of the main program he must specify the initial conditions
and a name for the file where the data will be saved
Third, at the end of the program, inside the function, he must define the first
derivative of each dependant variable as they are given by his equations. */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#define n 4 // Multitude of DE

void der(double der2[],double x,double y[]);
clock_t t1,t2;

int main(){
	
	t1=clock();
	double t;
	
	int i,c=0/*Multitude of time steps*/;
	double x0,x1=0;
	double y0[n],y1[n],y[n],k1[n],k2[n],k3[n],k4[n];
	
	// Define initial, final and step value for the independent variable
	double xi=0,xf=5E+8,h=10;
	// Define the length of time interval at the end of the numerical solution, which will be saved
	double a=1E+7;
	// The multitude of time steps to skip saving
	int p1=20;
	p1=p1+1;
	int p2=p1;
	// Initial conditions
	y0[0]=0.000004070;//unstable
	y0[1]=0.000086;
	y0[2]=0.00000009653;
	y0[3]=0.0000031989;
	// Open a txt file where the data will be saved 
	FILE *file=fopen("orbit.txt","w");
	if (file == NULL){
		printf("Error opening file!\n");
		exit(1);
	}
	
	x0=xi;
	while(x1<=xf+h){
		// Evaluate k1s
		der(k1,x0,y0);
		// Evaluate k2s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k1[i]*h/2;
		}
		der(k2,x0+h/2,y);
		// Evaluate k3s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k2[i]*h/2;
		}
		der(k3,x0+h/2,y);
		// Evaluate k4s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k3[i]*h;
		}
		der(k4,x0+h,y);
		// Evaluate y1s and x1
		for(i=0;i<n;i++){
			y1[i]=y0[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
		}
		x1=x0+h;
		c=c+1;
		
		if(x1>(xf-a)){
			if((p2%p1)==0){
				p2=p1;
				// Print results in the txt file
				fprintf(file,"%.15lf ",x1);
				for(i=0;i<n;i++){
					fprintf(file,"%.15lf ",y1[i]);
				}
				fprintf(file,"\n");
				p2=p2+1;
			}
			else{
				p2=p2+1;
			}
		}
		
		// Make new values, old
		x0=x1;
		for(i=0;i<n;i++){
			y0[i]=y1[i];
		}
	}
	
	fclose(file);// Close the txt file
	
	printf("The multitude of steps is %d",c);
	
	t2=clock();
	t=(double)(t2-t1)/CLOCKS_PER_SEC;
	printf("\nThe time needed for the program to run is %lf seconds",t);
	
	return 0;
}

// In this function the derivatives of the dependant variables are defined
void der(double der2[],double x,double y[]){
	// x -> independent variable
	// y[0],y[1],y[2],....,y[n-1] -> dependant variables
	/*
	 * y[0]=x1
	 * y[1]=x2
	 * y[2]=y1
	 * y[3]=y2
	*/
	//****Numerical values of the coefficients from the paper without dissipation*****
	double A=7.6983E-3,B=-177261.1194,C=-7.897728E-3,D=-5.021339E-3,E=8.113136E-7;
	double F=2.198007E-8,I=-4.897480E-9,K=6.9785637E-6,R=-5.3880293E-6,S=-1.2618983E-6;
	double eyilon=1;//Magnitude of perturbation
	der2[0]=2*A*y[2]+4*B*y[2]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*C*y[2]+E*y[3]-K*y[3]-2*R*y[2]+eyilon*1.4888230310773993E-7*y[0];
	der2[1]=2*A*y[3]+4*B*y[3]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*D*y[3]+E*y[2]-K*y[2]-2*S*y[3]-eyilon*2.205045054221079E-6*y[1];
	der2[2]=-2*A*y[0]-4*B*y[0]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*C*y[0]-E*y[1]-F-K*y[1]-2*R*y[0]-eyilon*1.4888230310773993E-7*y[2];
	der2[3]=-2*A*y[1]-4*B*y[1]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*D*y[1]-E*y[0]-I-K*y[0]-2*S*y[1]-eyilon*2.205045054221079E-6*y[3];
}
