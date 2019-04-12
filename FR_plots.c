#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
	double Imax;
	double H;
	double sigmaA;
	double sigmaJ;
	double muA;
	double muJ;
	double T;
	double sb;
	double sm;
} Parameters;

double functionalResponse(double x, Parameters *parms);
double fitness(double R, Parameters *parms, double *fractionAdults);
void linspace(double min, double max, int points, double *c);
double trapz(double x[], double y[], int n);
void readData(FILE *parmFile, double *R, double *FA, int *ct_ptr);
void fitnessFromData(double *Rvec, double *FA, Parameters *parms, double *fitness, int count);
void printFitnessResults(FILE *file, double *R_vec, double *fitness_vec, int count);

int main() {
	double *data1_Rup = (double*)malloc(500*sizeof(double));
	double *data1_FAup = (double*)malloc(500*sizeof(double));
	double *data1_Rdown = (double*)malloc(500*sizeof(double));
	double *data1_FAdown = (double*)malloc(500*sizeof(double));
	double *data2_Rup = (double*)malloc(500*sizeof(double));
	double *data2_FAup = (double*)malloc(500*sizeof(double));
	double *data2_Rdown = (double*)malloc(500*sizeof(double));
	double *data2_FAdown = (double*)malloc(500*sizeof(double));
	
	Parameters *parms1 = (Parameters*)malloc(sizeof(Parameters));
	Parameters *parms2 = (Parameters*)malloc(sizeof(Parameters));
	double fa1;
	double fa2;
	double *fractionAdults1 = &fa1;
	double *fractionAdults2 = &fa2;
	parms1->Imax = 1.0;
	parms1->H = 3.0;
	parms1->sigmaA = 0.35;
	parms1->sigmaJ = 0.5;
	parms1->muA = 0.015;
	parms1->muJ = 0.015;
	parms1->T = 0.1;
	parms1->sb = 5.0e-6;
	parms1->sm = 1e-4;
	
	parms2->Imax = 1.0;
	parms2->H = 3.0;
	parms2->sigmaA = 0.5;
	parms2->sigmaJ = 0.5;
	parms2->muA = 0.25;
	parms2->muJ = 0.015;
	parms2->T = 0.1;
	parms2->sb = 5.0E-6;
	parms2->sm = 1e-4;
	
	int count = 0;
	int *ct_ptr = &count; 
	double fitness_vec[10000];
	
	
	FILE *dataFile1_up = fopen("CASE18_sp1_FA_resources.dat","r");
	readData(dataFile1_up, data1_Rup, data1_FAup, ct_ptr);
	fclose(dataFile1_up);
	fitnessFromData(data1_Rup, data1_FAup, parms1, fitness_vec, count);
	FILE *file1 = fopen("CASE18_1UR.dat","w");
	printFitnessResults(file1, data1_Rup, fitness_vec, count);
	fclose(file1);
	/*
	count = 0;
	FILE *dataFile1_down = fopen("CASE13_sp1_down.dat","r");
	readData(dataFile1_down, data1_Rdown, data1_FAdown, ct_ptr);
	fclose(dataFile1_down);
	fitnessFromData(data1_Rdown, data1_FAdown, parms1, fitness_vec, count);
	FILE *file2 = fopen("CASE13_1DR.dat","w");
	printFitnessResults(file2, data1_Rdown, fitness_vec, count);
	fclose(file2);
	*/
	
	count = 0;
	FILE *dataFile2_up = fopen("CASE18_sp2_FA_resources.dat","r");
	readData(dataFile2_up, data2_Rup, data2_FAup, ct_ptr);
	fclose(dataFile2_up);
	fitnessFromData(data2_Rup, data2_FAup, parms2, fitness_vec, count);
	FILE *file3 = fopen("CASE18_2UR.dat","w");
	printFitnessResults(file3, data2_Rup, fitness_vec, count);
	fclose(file3);
	
	/*
	count = 0;
	FILE *dataFile2_down = fopen("CASE13_sp2_down.dat","r");
	readData(dataFile2_down, data2_Rdown, data2_FAdown, ct_ptr);
	fclose(dataFile2_down);
	fitnessFromData(data2_Rdown, data2_FAdown, parms2, fitness_vec, count);
	FILE *file4 = fopen("CASE13_2DR.dat","w");
	printFitnessResults(file4, data2_Rdown, fitness_vec, count);
	fclose(file4);
	 */
	/*
	FILE *file = fopen("FR.dat","w");
	
	int n = 100000;
	double R_vec[n];
	linspace(0.0, 7.0, n, R_vec);
	double y1[n];
	double y2[n];
	double R1_eq = 10.0;
	double R2_eq = 10.0;
	
	for (int i = n-1; i > -1; i--) {
		y1[i] = fitness(R_vec[i], parms1, fractionAdults1);
		if (fabs(y1[i]) < 1e-5 && isnormal(y1[i])) {
			R1_eq = R_vec[i];
		}
		y2[i] = fitness(R_vec[i], parms2, fractionAdults2);
		if (fabs(y2[i]) < 1e-5 && isnormal(y2[i])) {
			R2_eq = R_vec[i];
		}
		if (i%10 == 0) {
			fprintf(file,"%f %f %f\n",R_vec[i], y1[i], y2[i]);
			//fprintf(file,"%f %f %f\n",R_vec[i], fa1, fa2);
		}
	}
	fitness(R1_eq, parms1, fractionAdults1);
	printf("sp1 R1 eq: %f\n", R1_eq);
	printf("sp1 fraction Ad at eq: %f\n", fa1);
	fitness(R2_eq, parms2, fractionAdults2);
	printf("sp2 R1 eq: %f\n", R2_eq);
	printf("sp2 fraction Ad at eq: %f\n", fa2);
	*/
	free(parms1);
	free(parms2);
	free(data1_Rup);
	free(data2_Rup);
	free(data1_FAup);
	free(data2_FAup);
	free(data1_Rdown);
	free(data2_Rdown);
	free(data1_FAdown);
	free(data2_FAdown);
	//fclose(file);
	return 0;
};


void linspace(double min, double max, int points, double *c) {
	double h = (max-min)/(points-1);
	for (int i = 0; i < points; i++) {
		*(c+i) = min + h*i;
	}
}

double functionalResponse(double x, Parameters *parms) {
	double H = parms->H;
	double Imax = parms->Imax;
	return Imax*(x/(H+x));
}

double fitness(double R, Parameters *parms, double *fractionAdults) {
	double H = parms->H;
	double Imax = parms->Imax;
	double sigmaA = parms->sigmaA;
	double sigmaJ = parms->sigmaJ;
	double muA = parms->muA;
	double muJ = parms->muJ;
	double sm = parms->sm;
	double sb = parms->sb;
	double T = parms->T;
	double numerator;
	int n = 1000;
	double s_vec[n];
	linspace(sb,sm,n,s_vec);
	double birth, vgrowth, deathA, deathJ, K, fa;
	if (R > (H*T)/((sigmaA*Imax)-T)) {
		birth = (sm/sb)*(Imax*sigmaA*(R/(H+R))-T);
		deathA = muA;
	} else {
		birth = 0.0;
		deathA = muA - (Imax*sigmaA*(R/(H+R))-T);
	}
	if (R > (H*T)/((sigmaJ*Imax)-T)) {
		vgrowth = Imax*sigmaJ*(R/(H+R))-T;
		deathJ = muJ;
	} else {
		vgrowth = 0.0;
		deathJ = muJ - vgrowth;
	}
	double growth;
	double y_vec[n];
	for (int i = 0; i < n; i++) {
		y_vec[i] = pow((s_vec[i]),-((muJ/vgrowth)+1.0));
	}
	K = (birth/vgrowth)*trapz(s_vec,y_vec,n)*pow(sb,(muJ/vgrowth));
	numerator = (birth/(muA*sm*vgrowth))*pow((sm/sb),-muJ/vgrowth);
	//printf("%f\n",K);
	fa = numerator/(1.0+K);
	if (isnormal(K)) {
	*fractionAdults = fa;
	} else {
		//set it to previous normal value:
		fa = *fractionAdults;
	}
	if (!isnormal(birth)) {
		birth = 0.0;
	}
	//printf("%f\n",fa);
	return birth*fa - deathA*fa - deathJ*(1.0-fa);
}

//this needs to return an array for sure.
void fitnessFromData(double *Rvec, double *FA, Parameters *parms, double *fitness_vec, int count) {
	double H = parms->H;
	double Imax = parms->Imax;
	double sigmaA = parms->sigmaA;
	double sigmaJ = parms->sigmaJ;
	double muA = parms->muA;
	double muJ = parms->muJ;
	double sm = parms->sm;
	double sb = parms->sb;
	double T = parms->T;
	int n = 1000;
	double s_vec[n];
	linspace(5e-6,1e-4,n,s_vec);
	double birth, vgrowth, deathA, deathJ, K, fa, R;
	
	for (int i = 0; i < count; i++) {
		R = Rvec[i];
		fa = FA[i];
		if (R > (H*T)/((sigmaA*Imax)-T)) {
			birth = (sm/sb)*(Imax*sigmaA*(R/(H+R))-T);
			deathA = muA;
		} else {
			birth = 0.0;
			deathA = muA - (Imax*sigmaA*(R/(H+R))-T);
		}
		if (R > (H*T)/((sigmaJ*Imax)-T)) {
			vgrowth = Imax*sigmaJ*(R/(H+R))-T;
			deathJ = muJ;
		} else {
			vgrowth = 0.0;
			deathJ = muJ - vgrowth;
		}
		if (!isnormal(birth)) {
			birth = 0.0;
		}
		fitness_vec[i] = birth*fa - deathA*fa - deathJ*(1.0-fa);
	}
}

double trapz(double x[], double y[], int n) {
	double z = 0.0;
	double triangle = 0.0;
	double rectangle = 0.0;
	for (int i = 1; i < n; i++) {
		rectangle = y[i-1]*(x[i]-x[i-1]);
		triangle = 0.5*(x[i]-x[i-1])*(y[i]-y[i-1]);
		z = z + triangle + rectangle;
	}
	return z;
}

void readData(FILE *parmFile, double *R, double *FA, int *count) {
	int i = 0;
	char *throw_away;
	int array_size = 100;
	char char_array[array_size];
	char c;
	int indexR = 0;
	int indexFA = 0;
	
	for (int k = 0; k < array_size; k++) {
					char_array[k] = '\0';
				}
	
	while ((c = fgetc(parmFile))!=EOF) {
		if (i > 0) {
			switch(c) {
			case '-': char_array[i] = c; i++; break;
			case '.': char_array[i] = c; i++; break;
			case '1': char_array[i] = c; i++; break;
			case '2': char_array[i] = c; i++; break;
			case '3': char_array[i] = c; i++; break;
			case '4': char_array[i] = c; i++; break;
			case '5': char_array[i] = c; i++; break;
			case '6': char_array[i] = c; i++; break;
			case '7': char_array[i] = c; i++; break;
			case '8': char_array[i] = c; i++; break;
			case '9': char_array[i] = c; i++; break;
			case '0': char_array[i] = c; i++; break;
			case 'e': char_array[i] = c; i++; break;
			case 'E': char_array[i] = c; i++; break;
			case '+': char_array[i] = c; i++; break;
			default: i = 0; break;
		}
		} else {
			switch(c) {
			case '-': char_array[i] = c; i++; break;
			case '.': char_array[i] = c; i++; break;
			case '1': char_array[i] = c; i++; break;
			case '2': char_array[i] = c; i++; break;
			case '3': char_array[i] = c; i++; break;
			case '4': char_array[i] = c; i++; break;
			case '5': char_array[i] = c; i++; break;
			case '6': char_array[i] = c; i++; break;
			case '7': char_array[i] = c; i++; break;
			case '8': char_array[i] = c; i++; break;
			case '9': char_array[i] = c; i++; break;
			case '0': char_array[i] = c; i++; break;
			default: i = 0; break;
		}
		}
		
		if (c == ' ' || c == '\n') {
			//upon space check if there was a number in 2nd element of character array (not first because of e's)
			if (char_array[0]!='\0') {
				//convert the array contents into a double and store the number in an array:
				//fprintf(dest,"%e\n",strtod(char_array,&throw_away));
				i = 0;
				if (c == '\n') {
					FA[indexFA] = strtod(char_array,&throw_away);
					indexFA++;
					(*count)++;
				}
				if (c == ' ') {
					R[indexR] = strtod(char_array,&throw_away);
					indexR++;
				}
				//reset the elements of the char_array to null:
				for (int k = 0; k < array_size; k++) {
					char_array[k] = '\0';
				}
			}
		}
	}
}

void printFitnessResults(FILE *file, double *R_vec, double *fitness_vec, int count) {
	for (int i = 0; i < count; i++) {
		fprintf(file, "%f %f\n", R_vec[i], fitness_vec[i]);
	}
}