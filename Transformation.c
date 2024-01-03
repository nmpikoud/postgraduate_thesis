#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int main(){
	int n=4;
	double gamma1=1.088523E-5,gamma2=1.735596E-4;
	double b1=1.5217E-7,b2=1.92794E-6,mi1=1293.47,mi2=1293.47;
	
	//Open a file where the first variables will be saved
	FILE *file1=fopen("first_variables.txt","w");
	if (file1 == NULL){
    	printf("Error opening file!\n");
    	exit(1);
	}
	
	//Open the file where the time evolution of the fourth variables is stored
	FILE *file;
	file=fopen("orbit1.txt","r");
	
	double y[n];
	double z[n];
	double t;
	
	//Loop over the entries of the file where the data are stored
	while(fscanf(file,"%lf",&t)==1){
	
		fscanf(file,"%lf %lf %lf %lf",&y[0],&y[1],&y[2],&y[3]);
		
		//From xi and yi to Ii and ói
		//z1[0]=I1,z1[1]=I2,z1[2]=ó1,z1[3]=ó2
		z[0]=(y[0]*y[0]+y[2]*y[2])/2;
		z[1]=(y[1]*y[1]+y[3]*y[3])/2;
		// atan(x) returns the arc tangent of x in radians in the interval [-pi/2,+pi/2]
		z[2]=atan(y[2]/y[0]);
		z[3]=atan(y[3]/y[1]);
		
		//From Ii to Li and Gi
		//y[0]=L1,y[1]=L2,y[2]=G1,y[3]=G2
		y[0]=gamma1-z[0]-z[1];
		y[1]=gamma2+2*z[0]+2*z[1];
		y[2]=gamma1-2*z[0]-z[1];
		y[3]=gamma2+2*z[0]+z[1];
		
		//From Li and Gi to ai and ei
		z[0]=y[0]*y[0]/b1/b1/mi1;
		z[1]=y[1]*y[1]/b2/b2/mi2;
		z[2]=sqrt(1-y[2]*y[2]/y[0]/y[0]);
		z[3]=sqrt(1-y[3]*y[3]/y[1]/y[1]);
		fprintf(file1,"%lf %lf %lf %lf %lf\n",t,z[0],z[1],z[2],z[3]);
	}
	
	fclose(file);
	fclose(file1);
	
	return 0;
}
