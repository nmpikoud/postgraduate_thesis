/*This program implements the classical fourth-order Runge-Kutta method for
a system of ordinary differential equations using a constant step.
The user must do three things.
First, he must define the multide of the equations as a directive.
Second, at the first lines of the main program he must specify the initial conditions.
Third, at the end of the program, inside the function, he must define the first
derivative of each dependant variable as they are given by his equations.*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define n 4 // Multitude of DE

double *der(double,double y[]);

int main(){
	int i;
	int c1=0;//Multitude of time steps
	double x0,x1=0;
	double *pointer;
	double y0[n],y1[n],y[n],k1[n],k2[n],k3[n],k4[n];
	
	// Define initial, final and step value for the independent variable
	double xi=0,xf=10E+9,h=10;
	// Define the length of the time interval at the end of the numerical solution, which will be saved
	double a=1E+7;
	//The multitude of time steps to skip saving
	int p1=15;
	p1=p1+1;
	int p2=p1;
	
	// Define initial values for the dependant variables
	y0[0]=3.97089E-6;
	y0[1]=0.0001;
	y0[2]=-0.0000005;
	y0[3]=0.00005;
	
	// Open a txt file where the data will be saved 
	FILE *file=fopen("orbit.txt","w");
	if (file == NULL){
		printf("Error opening file!\n");
		exit(1);
	}
	
	//Print the initial conditions
	fprintf(file,"%.15lf ",xi);
		for(i=0;i<n;i++){
			fprintf(file,"%.15lf ",y0[i]);
		}
		fprintf(file,"\n");
	
	x0=xi;
	while(x1<=xf+h){
		// Evaluate k1s
		pointer=der(x0,y0);
		for(i=0;i<n;i++){
			k1[i]=*(pointer+i);
		}
		// Evaluate k2s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k1[i]*h/2;
		}
		pointer=der(x0+h/2,y);
		for(i=0;i<n;i++){
			k2[i]=*(pointer+i);
		}
		// Evaluate k3s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k2[i]*h/2;
		}
		pointer=der(x0+h/2,y);
		for(i=0;i<n;i++){
			k3[i]=*(pointer+i);
		}
		// Evaluate k4s
		for(i=0;i<n;i++){
			y[i]=y0[i]+k3[i]*h;
		}
		pointer=der(x0+h,y);
		for(i=0;i<n;i++){
			k4[i]=*(pointer+i);
		}
		// Evaluate y1s and x1
		for(i=0;i<n;i++){
			y1[i]=y0[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
		}
		x1=x0+h;
		c1=c1+1;
		
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
	
	printf("%d",c1);
	
	return 0;
}

// In this function the derivatives of the dependant variables are defined
double *der(double x,double y[]){
	// x -> independent variable
	// y[0],y[1],y[2],....,y[n-1] -> dependant variables
	static double der[n];
	//****Numerical values of the coefficients from the paper without dissipation*****
	double A=7.6983E-3,B=-177261.1194,C=-7.897728E-3,D=-5.021339E-3,E=8.113136E-7;
	double F=2.198007E-8,I=-4.897480E-9,K=6.9785637E-6,R=-5.3880293E-6,S=-1.2618983E-6;
	double b1=1.5217E-7,b2=1.92794E-7;
	double mi1=1293.47,mi2=1293.47;
	double b15=pow(b1,5),b25=pow(b2,5),mi12=pow(mi1,2),mi22=pow(mi2,2);
	double GAMMA1=pow(1.088523E-5,4),GAMMA2=pow(1.735596E-4,4);
	der[0]=2*A*y[2]+4*B*y[2]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*C*y[2]+E*y[3]-K*y[3]-2*R*y[2]-b15*mi12*y[0]/8/GAMMA1;
	der[1]=2*A*y[3]+4*B*y[3]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*D*y[3]+E*y[2]-K*y[2]-2*S*y[3]-b25*mi22*y[1]/8/GAMMA2;
	der[2]=-2*A*y[0]-4*B*y[0]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*C*y[0]-E*y[1]-F-K*y[1]-2*R*y[0]-b15*mi12*y[2]/8/GAMMA1;
	der[3]=-2*A*y[1]-4*B*y[1]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*D*y[1]-E*y[0]-I-K*y[0]-2*S*y[1]-b25*mi22*y[3]/8/GAMMA2;
	
	return der;
}
