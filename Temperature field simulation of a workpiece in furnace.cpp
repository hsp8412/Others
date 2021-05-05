/*
Author: Sipeng He && Wenhao Sun
Date: May 23, 2019
Program: Simulation of a workpiece's temperature field inside a furnace, using difference methods
Features: 
1.Let user pick the time length of the simulation
2.Simulation data is stored as C:\\temperatureField.dat C:\\keyPoint.dat
*/


#include <stdio.h>
/*Path of result data files*/
#define result1 "C:\\temperatureField.dat"
#define result2 "C:\\keyPoint.dat"
int main()
{
	/*Define the variables*/
	int R1=50,R2=100,N=80;
    int i,j=1;
    double DR=0.05,DZ=0.05,DS=0.005,TF=650,T0=20;
	double P=7.8,CP=0.5,L=0.3,H=0.2,h=0.03,S,S1;
	double FR,FZ,F1,F2,F3,F4;
	double RR1=DR*R1,RR2=DR*R2;
	/*Define the arrays*/
    double T1[200][200],T2[200][200];
	/*Define the file pointers*/
	FILE* fp1;
	FILE* fp2;
	/*Open(create) the data documents*/
    fp1=fopen(result1,"w+");
	fp2=fopen(result2,"w+");
	/*Calculation of variables*/
	FR=L*DS/(CP*P*DR*DR);
	FZ=L*DS/(CP*P*DZ*DZ);
	F1=(8*RR1+4*DR)/(4*RR1+DR);
	F2=8*RR1/(4*RR1+DR);
	F3=(8*RR2-4*DR)/(4*RR2-DR);
	F4=8*RR2/(4*RR2-DR);
	/*Initialization*/
	for(i=R1;i<=R2;i++)
		for(j=1;j<=N;j++)
		{T1[i][j]=T0;
		T2[i][j]=T0;}
	/*Input the duration of simulation*/
	printf("Input the time:");
	scanf("%lf",&S1);
	printf("S1=%lf\n",S1);
	printf("Press any key to continue...\n");
	getchar();
	getchar();
	/*Write the x-coordinates for keyPoints.dat*/
    fprintf(fp2,"S T[R1][1] T[(R1+R2)/2][1] T[R2][1] T[R1][N] T[(R1+R2)/2][N] T[R2][N]\n");
	/*Start the loop*/
	for(S=0.00;S<=S1+DS;S=S+DS)
	{
		/*Write the time of each line in keyPoints.dat*/
		fprintf(fp2,"%.3lf ",S);
		/*Write the temperature data in temperatureField.dat*/
		fprintf(fp2,"%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",T1[R1][1],T1[(R1+R2)/2][1],T1[R2][1],T1[R1][N],T1[(R1+R2)/2][N],T1[R2][N]);
		/*Write the time in temperatureField.dat, printing on the screen*/
		printf("S=%.3lf\n",S);
		fprintf(fp1,"S=%.3lf\n",S);
		fprintf(fp1,"%8.3lf ",S);
		printf("%8.3lf ",S);
		/*Write the x-coodinators in temperatureField.dat, printing on the screen*/
        for(i=R1;i<=R2;i++)
		{fprintf(fp1,"%8.2lf ",i*DR);
		printf("%8.2lf ",i*DR);}
		fprintf(fp1,"\n");
		printf("\n");
		/*Write the y-coordinator and data of each line, printing on the screen*/
		for(i=1;i<=N;i++)
		{
			fprintf(fp1,"%8.2lf",i*DZ);/*Write and Print y-coordinator*/
		    printf("%8.2lf",i*DZ);
			for(j=R1;j<=R2;j++)/*Write and Print each line of data*/
			{
			     printf("%8.2lf ",T2[j][i]);
			     fprintf(fp1,"%8.2lf ",T2[j][i]);
			}
		    fprintf(fp1,"\n");
		    printf("\n");
		}
		/*Differential equation*/
           T2[R1][1]=T1[R1][1]+F1*FR*(T1[R1+1][1]-T1[R1][1])+F2*FR*H*DR*(TF-T1[R1][1])/L+2*FZ*(T1[R1][2]-T1[R1][1]);/*The bottom left corner*/
		   T2[R2][1]=T1[R2][1]+F3*FR*(T1[R2-1][1]-T1[R2][1])+F4*FR*H*DR*(TF-T1[R2][1])/L+2*FZ*(T1[R2][2]-T1[R2][1]);/*The bottom right corner*/
           T2[R1][N]=T1[R1][N]+F1*FR*(T1[R1+1][N]-T1[R1][N])+F2*FR*H*DR*(TF-T1[R1][N])/L+2*FZ*(T1[R1][N-1]-T1[R1][N])+2*FZ*H*DZ*(TF-T1[R1][N])/L;/*The top left corner*/
           T2[R2][N]=T1[R2][N]+F3*FR*(T1[R2-1][N]-T1[R2][N])+F4*FR*H*DR*(TF-T1[R2][N])/L+2*FZ*(T1[R2][N-1]-T1[R2][N])+2*FZ*H*DZ*(TF-T1[R2][N])/L;/*The top right corner*/
		for(i=R1+1;i<R2;i++)
			for(j=2;j<N;j++)
        T2[i][j]=T1[i][j]+FR*(T1[i+1][j]+T1[i-1][j]-2*T1[i][j])+(DR/(2*DR*i))*FR*(T1[i+1][j]-T1[i-1][j])+FZ*(T1[i][j+1]+T1[i][j-1]-2*T1[i][j]);/*inner point*/
		for(j=2;j<N;j++)
        T2[R2][j]=T1[R2][j]+F3*FR*(T1[R2-1][j]-T1[R2][j])+F4*FR*(H*DR/L)*(TF-T1[R2][j])+FZ*(T1[R2][j+1]-2*T1[R2][j]+T1[R2][j-1]);/*The right boundary*/
		for(i=R1+1;i<R2;i++)
           T2[i][N]=T1[i][N]+FR*(T1[i+1][N]+T1[i-1][N]-2*T1[i][N])+(DR/(2*DR*i))*FR*(T1[i+1][N]-T1[i-1][N])+2*FZ*(T1[i][N-1]-T1[i][N])+(2*H*DZ/L)*FZ*(TF-T1[i][N]);/*The upper boundary*/
        for(i=R1+1;i<R2;i++)
           T2[i][1]=T1[i][1]+FR*(T1[i+1][1]+T1[i-1][1]-2*T1[i][1])+(DR/(2*DR*i))*FR*(T1[i+1][1]-T1[i-1][1])+2*FZ*(T1[i][2]-T1[i][1]);/*The lower boundary*/
        for(j=2;j<N;j++)
           T2[R1][j]=T1[R1][j]+F1*FR*(T1[R1+1][j]-T1[R1][j])+F2*FR*(H*DR/L)*(TF-T1[R1][j])+FZ*(T1[R1][j+1]-2*T1[R1][j]+T1[R1][j-1]);/*The left boundary*/
		   /*Operation of Arrays*/
		   for(i=R1;i<=R2;i++)
			for(j=1;j<=N;j++)
				T1[i][j]=T2[i][j];
}
	   /*Close the files*/
       fclose(fp1);
	   fclose(fp2);
	   getchar();
	   return 0;
	}
