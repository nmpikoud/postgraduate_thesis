#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#define n 4 // Multitude of DE

void der(double z[],double x,double y[]);
clock_t t1,t2;

int main(){
	
	t1=clock();
	double t;
	
	int i,j;
	int moic;// Variable that holds the multitude of initial conditions
	double x0,x1;
	double y0[n],y1[n],y[n],k1[n],k2[n],k3[n],k4[n];
	
	// Define initial, final and step value for the independent variable
	double xi=0,xf=2E+2,h=0.01;
	// Define the length of time interval at the end of the numerical solution, which will be saved
	double a=2E+2;
	// The multitude of time steps to skip saving
	int p1=0;
	p1=p1+1;
	int p2=p1;
	
	/* Open the file where the multitude of initial conditions is saved
	   read it and store it*/
	FILE *file1;
	file1=fopen("multitude_of_initial_conditions.txt","r");
	
	fscanf(file1,"%d",&moic);
	
	fclose(file1);
	printf("%d\n",moic);
	//
	
	double initial_conditions[moic][n];// Matrix that holds the initial conditions
	
	/* Open the file where the initial conditions are saved,
	   read them and store them in an two dimensional array */
	FILE *file2;
	file2=fopen("initial_conditions.txt","r");
	
	for(i=0;i<moic;i++){
		fscanf(file2,"%*d");
		for(j=0;j<n;j++){
			fscanf(file2," %lf",&initial_conditions[i][j]);
		}
	}
	
	fclose(file2);
	
	for(i=0;i<moic;i++){printf("%.15lf %.15lf %.15lf %.15lf\n",initial_conditions[i][0],initial_conditions[i][1],initial_conditions[i][2],initial_conditions[i][3]);}
	
	char filename_format[]="fourth_variables%d.txt";
	char filename[sizeof(filename_format) + 3];
	
		for(j=1;j<=moic;j++){
			snprintf(filename,sizeof(filename),filename_format,j);
			
			// Open a txt file where the data will be saved
			FILE *file=fopen(filename,"w");
			if (file == NULL){
				printf("Error opening file!\n");
				exit(1);
			}
	
			for(i=0;i<n;i++){
				y0[i]=initial_conditions[j-1][i];
				printf("%.15lf ",y0[i]);
			}
			printf("\n\n");
			
			getchar();
			
			x0=xi;
			x1=0;
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
				
				if(x1>(xf-a)){
					if(p2%p1==0){
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
	
		}
	
	t2=clock();
	t=(double)(t2-t1)/CLOCKS_PER_SEC;
	printf("The time needed for the program to run is %lf seconds.\n",t);
	
	return 0;
}

// In this function the derivatives of the dependant variables are defined
void der(double z[],double x,double y[]){
	// x -> independent variable
	// y[0],y[1],y[2],....,y[n-1] -> dependant variables
	/*
	 * y[0]=x1
	 * y[1]=x2
	 * y[2]=y1
	 * y[3]=y2
	*/
	// ****Numerical values of the coefficients from the paper without dissipation*****
	double A=7.6983E-3,B=-177261.1194,C=-7.897728E-3,D=-5.021339E-3,E=8.113136E-7;
	double F=2.198007E-8,I=-4.897480E-9,K=6.9785637E-6,R=-5.3880293E-6,S=-1.2618983E-6;
	double eyilon=1;// Magnitude of perturbation
	z[0]=2*A*y[2]+4*B*y[2]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*C*y[2]+E*y[3]-K*y[3]-2*R*y[2]+eyilon*1.4888230310773993E-7*y[0];
	z[1]=2*A*y[3]+4*B*y[3]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])+2*D*y[3]+E*y[2]-K*y[2]-2*S*y[3]-eyilon*2.205045054221079E-6*y[1];
	z[2]=-2*A*y[0]-4*B*y[0]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*C*y[0]-E*y[1]-F-K*y[1]-2*R*y[0]-eyilon*1.4888230310773993E-7*y[2];
	z[3]=-2*A*y[1]-4*B*y[1]*(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]+y[3]*y[3])-2*D*y[1]-E*y[0]-I-K*y[0]-2*S*y[1]-eyilon*2.205045054221079E-6*y[3];
}
