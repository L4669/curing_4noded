#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double autoclave_input_cycle(double time)
{
	int NP = 9, n = 0; //no of points available from the curing cycle
	double *x, *y, T;
	x = (double *)malloc(sizeof(double)*NP);
	y = (double *)malloc(sizeof(double)*NP);

	//input the cycle here or read from a file
	x[0] = 0*60; x[1] = 10*60; x[2] = 65*60; x[3] = 73*60; x[4] = 115*60;
	x[5] = 135*60; x[6] = 140*60; x[7] = 177*60; x[8] = 195*60;

	y[0] = 25; y[1] = 80; y[2] = 80; y[3] = 85; y[4] = 87; //90 & 89
	y[5] = 86; y[6] = 118; y[7] = 135; y[8] = 135;

	for(n = 1; n < NP; n++)
	{
		if(time <= x[n])
		{
			T =  y[n-1] + (y[n] - y[n-1])*(time - x[n-1])/(x[n] - x[n-1]);
			break;
		}
	}

	free(x); free(y);

	return T;
}

double autoclave_input_cycle_2(double time)
{
	int NP = 7, n = 0; //no of points available from the curing cycle
	double *x, *y, T;
	x = (double *)malloc(sizeof(double)*NP);
	y = (double *)malloc(sizeof(double)*NP);

	//input the cycle here or read from a file
	x[0] = 0*60; x[1] = 50*60; x[2] = 160*60; x[3] = 170*60; x[4] = 230*60;
	x[5] = 250*60; x[6] = 274*60;

	y[0] = 20; y[1] = 85; y[2] = 85; y[3] = 94; y[4] = 100;
	y[5] = 110; y[6] = 113;

	for(n = 1; n < NP; n++)
	{
		if(time <= x[n])
		{
			T =  y[n-1] + (y[n] - y[n-1])*(time - x[n-1])/(x[n] - x[n-1]);
			break;
		}
	}

	free(x); free(y);

	return T;
}


