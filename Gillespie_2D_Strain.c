#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dSFMT.h"    //random number generator, compile with dSFMT.c (check page Random_Numbers)

double det_bal_break(double u_a, double u_b, double e0, double mu, double rhokbt, double prf, double expon, double shft){
    return exp(-(e0+0.5*mu*rhokbt*(u_a*u_a + u_b*u_b)));//-prf*exp(-expon*(u-shft)*(u-shft)-diff);
}

double arrhenius(double u_a, double u_b, double k0, double kl2){
    return k0*exp(-kl2*(u_a*u_a + u_b*u_b));
}

void save_stress(int stress_cnt, double ts[], double stressa[], double stressb[], int repeat, char *folder, size_t flen){
    printf("ARGH stress\n");
	char filename[255];
    char fl[flen+15];
    strcpy(fl,folder);
    strcat(fl,"/stress_%d.dat\0");
	sprintf(filename, fl, repeat);
	FILE *f = fopen(filename, "ab+");
    int j;
    for(j=0; j<stress_cnt; j++){
        fprintf(f, "%0.15le %le %le\n", ts[j], stressa[j], stressb[j]);
    }
	fclose(f);
}

void save_X(double t, double Ua[], double Ub[], long int X[], long int X_T, int numN2, double du, char *folder){
    double nwrite[numN2+1];
    int j;
    for(j=0; j<numN2+1; j++){
        nwrite[j] = ((double)X[j])/(du*du*((double) X_T));
    }
	char filename[255];
	strcpy(filename,folder);
    strcat(filename,"/N_Final.dat\0");
	FILE *f = fopen(filename, "ab+");
    if(f==NULL){
        printf("Ugh\n");
    }
    if(nwrite[0] > 0){fprintf(f, "%f %f %le %le\n", -2.0, -2.0, nwrite[0]*du*du, t);}
    for(j=1; j<numN2+1; j++){
        fprintf(f, "%f %f %le %le\n", Ua[j-1], Ub[j-1], nwrite[j], t);
    }
	fclose(f);
}

dsfmt_t  dsfmt;

void save_param(char *folder, double e0, double mu, double rho, double prf, double expon, double shft, long int X_T, 
                double k_a, double k_d, int numN, double qu_range, double du){

	char filename[255];
    strcpy(filename,folder);
    strcat(filename,"/parameters.dat\0");

	FILE *f = fopen(filename, "ab+");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }


    fprintf(f, "Number of molecules: %ld \nNumber of bins: %d \nRange: %le \ndu: %le\n", X_T, numN, qu_range, du);
    fprintf(f, "mu: %f \ne0: %le \nrho: %le \nk_a^0: %le\nk_d: %le\n", mu, e0, rho, k_a, k_d);
    fprintf(f, "Prefactor: %f \nExponent: %le \nShift: %le\n", prf, expon, shft);
	fclose(f);
}

double start_phi_b = 0.75;
int numN = 50;

double qu_range = 0.25;

long int X_T = 10000000;

double k_a = 4E6*125E-6;

double t_max = 50000.0;

double rho = 530E18;
double kbT = 300*1.381E-23;
double mu = 500;
double e0 = -9;

double k_cross = 3E-4;
double l_cross = 40E-9;

int main(int argc, char **argv){

    int rpt;
    double rhokbt = 1/(rho*kbT);
    //double kl2 = k_cross * l_cross * l_cross/kbT;
    double k_d = k_a / det_bal_break(0, 0, e0, mu, rhokbt, 0, 0, 0);
    for(rpt=0;rpt<1;rpt++){

        unsigned long int cycles = 0;
        double t = 0;

        char *arg1, *arg2, *arg3;

        double prf = strtod(argv[1], &arg1);
        double expon = strtod(argv[2], &arg2);
        double shft = strtod(argv[3], &arg3);
        char *folder = argv[4];

        size_t flen = strlen(folder);

        printf("%le %le %le %s\n", prf, expon, shft, folder);

        int seed;
        seed=time(NULL);
        dsfmt_init_gen_rand(&dsfmt, seed);

        int numN2 = numN*numN;

        double Ua[numN2];
        double Ub[numN2];
        double qu[numN];

        double n_u_q[numN2];

        long int X[numN2+1];
        long int V[numN2*2]; // For each reaction, number is index of rxn, +ve is gain
        double C[numN2*2];

        double phibm1 = 1;

        int i;
        double du = 2*qu_range/(numN-1);
        // Set up values of U and Q in a grid
        for(i=0; i<numN; i++){
            qu[i] = -qu_range+du*i;
        }
        for(i=0; i<numN2; i++){
            Ua[i] = qu[i%numN];
            Ub[i] = qu[(int)floor(i/numN)];
            phibm1 += det_bal_break(Ua[i], Ub[i], e0, mu, rhokbt, prf, expon, shft)*du*du;
        }
        for(i=0; i<numN2; i++){
            n_u_q[i] = det_bal_break(Ua[i], Ub[i], e0, mu, rhokbt, prf, expon, shft)/(phibm1);//(1-start_phi_b)/numN; // Initialise values of n(u,q)
        }

        // Initialise X
        X[0] = X_T;
        for(i=1; i<numN2+1; i++){
            X[i] = (long int) floor(X_T*n_u_q[i-1]*du*du);
            X[0] -= X[i];
        }

        // Initialise V & C

        for(i=0; i<numN2; i++){
            V[i*2] = (long int)(i+1);
            V[(i*2)+1] = (long int) -(i+1);
            C[(i*2)+1] = k_d;
            C[i*2] = C[(i*2)+1]*det_bal_break(Ua[i], Ub[i], e0, mu, rhokbt, prf, expon, shft);
        }

        double a[numN2*2];
        double cs, a0;
        double r1, r2;
        double r2a0, tau;
        long int cs_cnt;

        int ar_len = 50000;
        double stressa[ar_len];
        double stressb[ar_len];
        double ts[ar_len];
        int stress_cnt = 0;
        double str_sum_a, str_sum_b;

        save_param(folder, e0, mu, rho, prf, expon, shft, X_T, k_a, k_d, numN, qu_range, du);

        clock_t start_time = clock();
        double elapsed_time;
        int rxn_cnt;

        while(t < t_max){

            if(cycles%10000000==0){
                printf("%f\n", t);
                printf("Unbound num: %d\n", X[0]);
                elapsed_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
                printf("Elapsed time %lf\n", elapsed_time);
                start_time = clock();
            }
            if(cycles%10000000==0){save_X(t, Ua, Ub, X, X_T, numN2, du, folder);}
            if(cycles>1250&&cycles%25000==0){
                if(stress_cnt>=ar_len){
                    save_stress(stress_cnt, ts, stressa, stressb, rpt, folder, flen);
                    stress_cnt = 0;
                }
                else{
                    str_sum_a = 0;
                    str_sum_b = 0;
                    for(i=0; i<numN2; i++){
                        str_sum_a += sqrt((double) X_T)*mu*(Ua[i])*((double) X[i+1])/((double) X_T);
                        str_sum_b += sqrt((double) X_T)*mu*(Ub[i])*((double) X[i+1])/((double) X_T);
                    }
                    ts[stress_cnt] = t;
                    stressa[stress_cnt] = str_sum_a;
                    stressb[stress_cnt] = str_sum_b;
                    stress_cnt++;
                }
            }


            a0 = 0;

            for(i=0; i<numN2; i++){ // Calculate a=c*x and sum of all a
                a[i*2] = C[i*2]*((double) X[0]);
                a0 += a[i*2];
                a[(i*2)+1] = C[(i*2)+1]*((double) X[i+1]);
                a0 += a[(i*2)+1];
            }

            r1 = dsfmt_genrand_open_open(&dsfmt);
            r2 = dsfmt_genrand_open_open(&dsfmt);

            tau = log(1/r1)/a0;
            t += tau;

            //if(cycles%100000==0){printf("%f %f %ld\n", a[0], C[0], X[1]);}

            r2a0 = r2*a0;
            cs = 0.0;
            cs_cnt = 0;

            while(cs<r2a0){
                cs += a[cs_cnt];
                cs_cnt++;
            }
            cs_cnt--; // Correction

            if(V[cs_cnt] < 0){
                X[0]++;
                X[-V[cs_cnt]]--;
            }
            else{
                X[0]--;
                X[V[cs_cnt]]++;
            }

            //if(cycles%50000==0){printf(" %ld", V[cs_cnt]);}

            cycles++;

        }
        save_stress(stress_cnt, ts, stressa, stressb, rpt, folder, flen); // Final catch of stresses
        save_X(t, Ua, Ub, X, X_T, numN2, du, folder);
    }

    return 0;
}
