#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dSFMT.h"    //random number generator, compile with dSFMT.c (check page Random_Numbers)

double det_bal_break(double u, double q, double e0, double mu, double D, 
            double chi, double rhokbt, double prf, double expon, double shft){

    return exp(-(e0+0.5*rhokbt*(mu*u*u + 2*D*u*q + chi*q*q)))-prf*exp(-expon*(u-shft)*(u-shft)-diff);
}

double arrhenius(double u_a, double u_b, double k0, double kl2){
    return k0*exp(-0.5*kl2*(u_a*u_a + u_b*u_b));
}

void save_stress(int stress_cnt, double ts[], double stress[], int repeat, char *folder, size_t flen){
    printf("ARGH stress\n");
	char filename[255];
    char fl[flen+15];
    strcpy(fl,folder);
    strcat(fl,"/stress_%d.dat\0");
	sprintf(filename, fl, repeat);
	FILE *f = fopen(filename, "ab+");
    int j;
    for(j=0; j<stress_cnt; j++){
        fprintf(f, "%0.15le %le\n", ts[j], stress[j]);
    }
	fclose(f);
}

void save_X(double t, double U[], double Q[], long int X[], long int X_T, int numN2, double du, char *folder){
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
        fprintf(f, "%f %f %le %le\n", U[j-1], Q[j-1], nwrite[j], t);
    }
	fclose(f);
}

/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
/* RAND_MAX = 2147483647 */
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}


dsfmt_t  dsfmt; // Must be declared before poisson_single()

int poisson_Single(double lambda){
    double exp_lambda = exp(-lambda);
    double randUni; //uniform variable
    double prodUni; //product of uniform variables
    int randPoisson; //Poisson variable
    
    //initialize variables
    randPoisson = -1;
    prodUni = 1;
    do
    {
        randUni = dsfmt_genrand_open_open(&dsfmt); //generate uniform variable
        prodUni = prodUni * randUni; //update product
        randPoisson++; //increase Poisson variable

    } while (prodUni > exp_lambda);
    return randPoisson;
}

void save_param(char *folder, double e0, double mu, double D, double chi,
                double rho, double prf, double expon, double shft, long int X_T, 
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
    fprintf(f, "mu: %f \nD: %f \nchi: %f \e0: %le \nrho: %le \nk_a^0: %le\nk_d: %le\n", mu, D, chi, e0, rho, k_a, k_d);
    fprintf(f, "Prefactor: %f \nExponent: %le \nShift: %le\n", prf, expon, shft);
	fclose(f);
}

int numN = 50;

double qu_range = 0.25;

long int X_T = 100000000;

double k_a = 4E6*125E-6;

double t_max = 50000.0;

double rho = 530E18;
double kbT = 300*1.381E-23;
double mu = 500;
double D = 50;
double chi = 1000;
double e0 = -9;

double k_cross = 3E-4;
double l_cross = 40E-9;

int main(int argc, char **argv){

    int rpt;
    double rhokbt = 1/(rho*kbT);
    //double kl2 = k_cross * l_cross * l_cross/kbT;
    double k_d = k_a / det_bal_break(0, 0, e0, mu, D, chi, rhokbt, 0, 0, 0);
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

        double U[numN2];
        double Q[numN2];
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
            U[i] = qu[i%numN];
            Q[i] = qu[(int)floor(i/numN)];
            phibm1 += det_bal_break(U[i], Q[i], e0, mu, D, chi, rhokbt, prf, expon, shft)*du*du;
        }
        for(i=0; i<numN2; i++){
            n_u_q[i] = det_bal_break(U[i], Q[i], e0, mu, D, chi, rhokbt, prf, expon, shft)/(phibm1); // Initialise values of n(u,q)
        }

        // Initialise X
        X[0] = X_T;
        for(i=1; i<numN2+1; i++){
            X[i] = (long int) (floor(X_T*n_u_q[i-1]*du*du));
            X[0] -= X[i];
        }

        // Initialise V & C

        for(i=0; i<numN2; i++){
            V[i*2] = (long int)(i+1);
            V[(i*2)+1] = (long int) -(i+1);
            C[(i*2)+1] = k_d;
            C[i*2] = C[(i*2)+1]*det_bal_break(U[i], Q[i], e0, mu, D, chi, rhokbt, prf, expon, shft);
        }

        double a[numN2*2];
        double cs, a0;
        double r1, r2;
        double r2a0, tau;
        long int cs_cnt;

        int ar_len = 50000;
        double stress[ar_len];
        double ts[ar_len];
        int stress_cnt = 0;
        double str_sum;

        int j;

        double tau1, tau2;

        save_param(folder, e0, mu, D, chi, rho, prf, expon, shft, X_T, k_a, k_d, numN, qu_range, du);

        clock_t start_time = clock();
        double elapsed_time;
        int rxn_cnt;

        while(t < t_max){
            //elapsed_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
            //printf("Elapsed time %lf\n", elapsed_time);
            //start_time = clock();

            if(cycles%10000000==0){
                printf("%f\n", t);
                printf("Unbound num: %ld\n", X[0]);
                elapsed_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
                printf("Elapsed time %lf\n", elapsed_time);
                start_time = clock();
            }
            if(cycles%10000000==0){save_X(t, U, Q, X, X_T, numN2, du, folder);}
            if(cycles>1250&&cycles%5000==0){
                if(stress_cnt>=ar_len){
                    save_stress(stress_cnt, ts, stress, rpt, folder, flen);
                    stress_cnt = 0;
                }
                else{
                    str_sum = 0;
                    for(i=0; i<numN2; i++){
                        str_sum +=  sqrt((double) X_T)*(mu*U[i]+D*Q[i])*((double) X[i+1])/((double) X_T);
                    }
                    ts[stress_cnt] = t;
                    stress[stress_cnt] = str_sum;
                    stress_cnt++;
                }
            }

            /* Gillespie */

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

            r2a0 = r2*a0;
            cs = 0.0;
            cs_cnt = 0;

            while(cs<r2a0){
                cs += a[cs_cnt];
                cs_cnt++;
            }
            cs_cnt--; // correction

            if(V[cs_cnt] < 0){
                X[0]++;
                X[-V[cs_cnt]]--;
            }
            else{
                X[0]--;
                X[V[cs_cnt]]++;
            }

                cycles++;
        }
        save_stress(stress_cnt, ts, stress, rpt, folder, flen); // Final catch of stresses
        save_X(t, U, Q, X, X_T, numN2, du, folder);
    }

    return 0;
}
