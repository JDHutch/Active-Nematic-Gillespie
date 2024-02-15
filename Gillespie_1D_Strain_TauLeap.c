#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dSFMT.h"    //random number generator, compile with dSFMT.c (check page Random_Numbers)

double det_bal_break(double u, double e0, double mu, double rhokbt, double prf, double expon, double shft){
    return exp(-(e0+0.5*mu*u*u*rhokbt))-prf*exp(-expon*(u-shft)*(u-shft));
}

double arrhenius(double u, double k0, double kl2){
    return k0*exp(-0.5*kl2*u*u);
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
    printf("ARGH stress\n");
}

void save_X(double t, double U[], long int X[], long int X_T, int numN, double du, char *folder){
    double nwrite[numN+1];
    int j;
    for(j=0; j<numN+1; j++){
        nwrite[j] = ((double)X[j])/(du*((double) X_T));
    }
	char filename[255];
    //char fl[flen+15];
    strcpy(filename,folder);
    strcat(filename,"/N_Final.dat\0");
	//sprintf(filename, fl);
	FILE *f = fopen(filename, "ab+");
    if(nwrite[0] > 0){fprintf(f, "%f %le %le\n", -2.0, nwrite[0]*du, t);}
    for(j=1; j<numN+1; j++){
        fprintf(f, "%f %le %le\n", U[j-1], nwrite[j], t);
    }
	fclose(f);
}

/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
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


dsfmt_t  dsfmt;

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


double start_phi_b = 0.75;
int numN = 50;

double qu_range = 0.25;

long int X_T = 100000000;

double k_a = 4E6*125E-6;

double t_max = 50000000.0;

double rho = 530E18;
double kbT = 300*1.381E-23;
double mu = 500;
double e0 = -9;

double k_cross = 3E-4;
double l_cross = 40E-9;

int N_crit = 10;
double e_control = 0.03;
double zeta = 10.0;

int main(int argc, char **argv){

    int rpt;
    double rhokbt = 1/(rho*kbT);
    double kl2 = k_cross * l_cross * l_cross/kbT;
    double X_T_root = sqrt((double) X_T);
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

        double U[numN];
        double qu[numN];

        double n_u_q[numN];

        long int X[numN+1];
        long int V[numN*2]; // For each reaction, number is index of rxn, +ve is gain
        double C[numN*2];

        int crit[numN*2];
        int noncrit[numN*2];

        double mu_i[numN+1];
        double s2_i[numN+1];

        int n_crit;
        int n_noncrit;

        double phibm1 = 1;

        int i;
        double du = 2*qu_range/(numN-1);
        // Set up values of U and Q in a grid
        for(i=0; i<numN; i++){
            qu[i] = -qu_range+du*i;
        }
        for(i=0; i<numN; i++){
            U[i] = qu[i%numN];
            phibm1 += det_bal_break(U[i], e0, mu, rhokbt, prf, expon, shft)*du;
        }
        for(i=0; i<numN; i++){
            n_u_q[i] = det_bal_break(U[i], e0, mu, rhokbt, prf, expon, shft)/(phibm1);//(1-start_phi_b)/numN; // Initialise values of n(u,q)
        }

        // Initialise X
        X[0] = X_T;
        for(i=1; i<numN+1; i++){
            X[i] = (long int) (floor(X_T*n_u_q[i-1]*du));
            X[0] -= X[i];
        }

        //X[0] = floor(X[0]*du);
        //X_T = floor(X_T*du);

        // Initialise V & C

        for(i=0; i<numN; i++){
            V[i*2] = (long int)(i+1);
            V[(i*2)+1] = (long int) -(i+1);
            C[i*2] = arrhenius(U[i], k_a, kl2);
            C[(i*2)+1] = C[i*2]/det_bal_break(U[i], e0, mu, rhokbt, prf, expon, shft);
        }

        double a[numN*2];
        double cs, a0;
        double r1, r2;
        double r2a0, tau;
        long int cs_cnt;

        int ar_len = 50000;
        double stress[ar_len];
        double ts[ar_len];
        int stress_cnt = 0;
        double str_sum;

        int rj;
        int rct;
        double tau_prime;
        double tau_dprime;
        double sum_mean; // For our system two denominators mu and sig^2 are equal
        double max_xi;
        double tau_tentative;
        double a0c;
        double k[numN*2];
        long int X_temp[numN+1];
        int x_neg = 0;

        int j;

        while(t < t_max){

            if(cycles%1000000==0){printf("%f\n", t);}
            if(cycles%1000000==0){save_X(t, U, X, X_T, numN, du, folder);}
            if(cycles>1250&&cycles%500==0){
                if(stress_cnt>=ar_len){
                     save_stress(stress_cnt, ts, stress, rpt, folder, flen);
                    stress_cnt = 0;
                }
                else{
                    str_sum = 0;
                    for(i=0; i<numN; i++){
                        str_sum += X_T_root*mu*(U[i])*((double) X[i+1])/((double) X_T);
                    }
                    ts[stress_cnt] = t;
                    stress[stress_cnt] = str_sum;
                    stress_cnt++;
                }
            }

            /* Gillespie */

            a0 = 0;
            mu_i[0] = 0; // Reset divisors for tprime
            s2_i[0] = 0;

            for(i=0; i<numN; i++){ // Calculate a=c*x and sum of all a
                a[i*2] = C[i*2]*((double) X[0]);
                a0 += a[i*2];
                a[(i*2)+1] = C[(i*2)+1]*((double) X[i+1]);
                a0 += a[(i*2)+1];

                mu_i[i+1] = 0;
                s2_i[i+1] = 0;
            }
            


            /* Identify critial reactions */
            n_crit = 0;
            n_noncrit = 0;


            for(i=0; i<numN*2; i++){
                rj = V[i];
                if(rj < 0){
                    rct = X[-rj];
                }
                else{
                    rct = X[0];
                }
                if(rct != 0){
                    if(rct < N_crit){
                        crit[n_crit] = i;
                        n_crit++;
                    }
                    else{
                        noncrit[n_noncrit] = i;
                        n_noncrit++;
                        if(rj < 0){
                            mu_i[0] += a[i];
                            mu_i[-rj] -= a[i];
                            s2_i[-rj] += a[i];
                        }
                        else{
                            mu_i[0] -= a[i];
                            mu_i[rj] += a[i];
                            s2_i[rj] += a[i];
                        }
                        s2_i[0] += a[i];
                    }
                }
            }

            /* Find max tau */
            tau_prime = 1000.0;
            for(i=0; i<numN+1; i++){ // numN+1
                max_xi = ((double)X[i])*e_control > 1.0 ? ((double)X[i])*e_control : 1.0;
                if(mu_i[i]<0){mu_i[i] = -mu_i[i];}
                tau_tentative = max_xi/mu_i[i] < (max_xi*max_xi)/s2_i[i] ? max_xi/mu_i[i] : (max_xi*max_xi)/s2_i[i];
                tau_prime = tau_tentative < tau_prime ? tau_tentative : tau_prime;
            }

            for(j=0;j<n_noncrit;)




            /*for(j=0; j<n_noncrit; j++){
                i = noncrit[j];
                max_xi = ((double)X[i+1])*e_control > 1.0 ? ((double)X[i+1])*e_control : 1.0;
                sum_mean = a[(noncrit[j]*2)] + a[(noncrit[j]*2)+1];
                tau_tentative = max_xi/sum_mean < (max_xi*max_xi)/sum_mean ? max_xi/sum_mean : (max_xi*max_xi)/sum_mean;
                tau_prime = tau_tentative < tau_prime ? tau_tentative : tau_prime;
            }*/

            /* Should we do normal Gillespie? */

            //printf("ARgh\n");

            if(tau_prime < zeta/a0){

                //printf("What's the fucking point\n");

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

            else{

                /* Calculate tau double prime */
                // Time that a critical reaction would take place in
                a0c = 0;
                x_neg = 0;

                if(n_crit > 0){
                    for(i=0; i<n_crit; i++){ // Calculate a=c*x and sum of all a
                        a0c += a[crit[i]];
                        //if(crit[i]%2==0){a0c += a[crit[i]];}
                        //else{a0c += a[crit[i]+1];}
                    }
                    r1 = dsfmt_genrand_open_open(&dsfmt);
                    tau_dprime = log(1/r1)/a0c;
                }
                else{
                    tau_dprime = 10000;
                }

                do{
                    x_neg = 0;

                    tau = tau_prime < tau_dprime ? tau_prime : tau_dprime; // Leap is smaller of two 

                    // k determines how many firings of each reaction in tau leap
                    for(i=0; i<numN*2; i++){
                            k[i] = 0;
                    }

                    // Do Gilespie for one critical reaction if leap large enough for 1 critical to occur
                    if(tau_prime > tau_dprime){ 
                        r2 = dsfmt_genrand_open_open(&dsfmt);
                        r2a0 = r2*a0c;
                        cs = 0.0;
                        cs_cnt = 0;

                        while(cs<r2a0){
                            cs += a[crit[cs_cnt]];
                            cs_cnt++;
                        }
                        cs_cnt--;

                        k[crit[cs_cnt]] = 1;

                    }

                    // For non critical, assign random poisson distributed number of firings with mean/variance a(X)tau
                    for(i=0; i<n_noncrit; i++){
                            k[noncrit[i]] = poisson_Single(tau*a[noncrit[i]]);
                    }

                    int assign_order[n_noncrit];
                    for(i=0; i<n_noncrit; i++){assign_order[i]=i;}

                    shuffle(assign_order, n_noncrit);

                    // Create temporary populations to check for negatives -> rejection of step
                    for(i=0; i<numN+1; i++){
                        X_temp[i] = X[i];
                    }

                    for(j=0; j<n_noncrit; j++){
                        i = noncrit[assign_order[j]];
                        if(V[i] < 0){
                            X_temp[0] += k[i];
                            X_temp[-V[i]] -= k[i];

                            if(X_temp[-V[i]] < 0){
                                x_neg = 1;
                                break;
                            }
                        }
                        else{
                            X_temp[0] -= k[i];
                            X_temp[V[i]] += k[i];

                            if(X_temp[0] < 0){
                                x_neg = 1;
                                break;
                            }
                        }
                    }

                    // Population gone negative -> tau leap too large
                    if(x_neg){
                        tau_prime /= 2;
                    }

                }while(x_neg);

                // Affirm new populations
                for(i=0; i<numN*2; i++){
                    X[i] = X_temp[i];
                }

                t += tau;
                cycles++;
            }
        }
         save_stress(stress_cnt, ts, stress, rpt, folder, flen); // Final catch of stresses
        save_X(t, U, X, X_T, numN, du, folder);
    }

    return 0;
}
