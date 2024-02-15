#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dSFMT.h"    //random number generator, compile with dSFMT.c (check page Random_Numbers)

double det_bal_break(double u_a, double u_b, double q_a, double q_b, double e0, double mu,
        double D, double chi, double rhokbt, double prf, double expon, double shft){

                return exp(-(e0+rhokbt*(mu*(u_a*u_a + u_b*u_b) + 
                2*D*(u_a*q_a + u_b*q_b) + chi*(q_a*q_a + q_b*q_b))));//-prf*exp(-expon*(u-shft)*(u-shft)-diff);
}

double arrhenius(double u_a, double u_b, double k0, double kl2){
    return k0*exp(-0.5*kl2*(u_a*u_a + u_b*u_b));
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

void save_X(double t, double Ua[], double Ub[], double Qa[], double Qb[], long int X[], 
            long int X_T, int numN4, double du, char *folder, double C[]){
    double* nwrite = malloc((numN4+1)*sizeof(double));
    int j;
    for(j=0; j<numN4+1; j++){
        nwrite[j] = ((double)X[j]);///(du*du*du*du*((double) X_T));
    }
	char filename[255];
	strcpy(filename,folder);
    strcat(filename,"/N_Final.dat\0");
	FILE *f = fopen(filename, "ab+");
    if(f==NULL){
        printf("Ugh\n");
    }
    if(nwrite[0] > 0){fprintf(f, "%f %f %f %f %le %le %f %f\n", -2.0, -2.0, -2.0, -2.0, nwrite[0], t, 0.0, 0.0);} //*du*du*du*du
    for(j=1; j<numN4+1; j++){
        fprintf(f, "%f %f %f %f %le %le %f %f\n", Ua[j-1], Ub[j-1], Qa[j-1], Qb[j-1], nwrite[j], t, C[2*(j-1)], C[2*(j-1)+1]);
    }
	fclose(f);
    free(nwrite);
}

void save_n(double t, double Ua[], double Ub[], double Qa[], double Qb[], double n[], 
            long int X_T, int numN4, double du, char *folder){
    double* nwrite = malloc((numN4+1)*sizeof(double));
	char filename[255];
	strcpy(filename,folder);
    strcat(filename,"/N_Initial.dat\0");
	FILE *f = fopen(filename, "ab+");
    if(f==NULL){
        printf("Ugh\n");
    }
    int j;
    for(j=0; j<numN4; j++){
        fprintf(f, "%f %f %f %f %f\n", Ua[j-1], Ub[j-1], Qa[j-1], Qb[j-1], n[j]);
    }
	fclose(f);
    free(nwrite);
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
    fprintf(f, "mu: %f \nD: %f \nchi: %f \ne0: %le \nrho: %le \nk_a^0: %le\nk_d: %le\n", mu, D, chi, e0, rho, k_a, k_d);
    fprintf(f, "Prefactor: %f \nExponent: %le \nShift: %le\n", prf, expon, shft);
	fclose(f);
}

int numN = 20;

double qu_range = 0.4;

long int X_T = 1E9;

double k_a = 4E6*125E-6;

double t_max = 100000.0;

double rho = 530E18;
double kbT = 300*1.381E-23;
double mu = 1000;
double D = 50;
double chi = 1000;
double e0 = -9;

double k_cross = 3E-4;
double l_cross = 40E-9;

int N_crit = 10;
double e_control = 0.08;
double zeta = 10.0;

int stress_sample = 5;

int main(int argc, char **argv){

    int rpt;
    double rhokbt = 1/(rho*kbT);
    //double kl2 = k_cross * l_cross * l_cross/kbT;
    double k_d = k_a / det_bal_break(0, 0, 0, 0, e0, mu, D, chi, rhokbt, 0, 0, 0);
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

        int numN2 = numN*numN;
        int numN3 = numN*numN*numN;
        int numN4 = numN*numN*numN*numN;


        double* Ua =  malloc(numN4*sizeof(double));
        double* Ub =  malloc(numN4*sizeof(double));
        double* Qa =  malloc(numN4*sizeof(double));
        double* Qb =  malloc(numN4*sizeof(double));
        double qu[numN];

        double* n_u_q = malloc(numN4*sizeof(double));
        double* n_test = malloc(numN4*sizeof(double));

        long int* X = malloc((numN4+1)*sizeof(long int));
        long int* V = malloc(numN4*2*sizeof(long int)); // For each reaction, number is index of rxn, +ve is gain
        double* C = malloc(numN4*2*sizeof(double)); 

        long int* crit = malloc((numN4+1)*sizeof(long int));
        long int* noncrit = malloc((numN4+1)*sizeof(long int));

        double* mu_i = malloc((numN4+1)*sizeof(double));
        double* s2_i = malloc((numN4+1)*sizeof(double));

        int n_crit;
        int n_noncrit;
        
        double phibm1 = 1;

        int i;
        double du = 2*qu_range/(numN-1);
        double du4 = du*du*du*du;
        // Set up values of U and Q in a grid
        for(i=0; i<numN; i++){
            qu[i] = -qu_range+du*i;
        }
        for(i=0; i<numN4; i++){
            Ua[i] = qu[(int)floor(i/numN3)];
            Ub[i] = qu[(int)(floor(i/numN2)-floor(i/numN3)*numN)];
            Qa[i] = qu[(int)(floor(i/numN)-floor(i/numN2)*numN)];
            Qb[i] = qu[(int)(i-floor(i/numN)*numN)];
            phibm1 += det_bal_break(Ua[i], Ub[i], Qa[i], Qb[i], e0, mu, D, chi, rhokbt, prf, expon, shft)*du4;
        }
        double phib = 0;
        for(i=0; i<numN4; i++){
            n_u_q[i] = det_bal_break(Ua[i], Ub[i], Qa[i], Qb[i], e0, mu, D, chi, rhokbt, prf, expon, shft)/(phibm1);//(1-start_phi_b)/numN; // Initialise values of n(u,q)
            phib += n_u_q[i]*du4;
        }

        // Initialise X
        X[0] = X_T;
        for(i=1; i<numN4+1; i++){
            X[i] = (long int) (floor(X_T*n_u_q[i-1]*du4));
            X[0] -= X[i];
        }
        save_n(t, Ua, Ub, Qa, Qb, n_u_q, X_T, numN4, du, folder);
        free(n_u_q);

        //X[0] = floor(X[0]*du);
        //X_T = floor(X_T*du);

        // Initialise V & C

        for(i=0; i<numN4; i++){
            V[i*2] = (long int)(i+1);
            V[(i*2)+1] = (long int) -(i+1);
            C[(i*2)+1] = k_d;
            C[i*2] = C[(i*2)+1]*det_bal_break(Ua[i], Ub[i], Qa[i], Qb[i], e0, mu, D, chi, rhokbt, prf, expon, shft);
        }
        
        double* a = malloc(numN4*2*sizeof(double)); 
        double cs, a0;
        double r1, r2;
        double r2a0, tau;
        long int cs_cnt, crit_cnt;

        int ar_len = 5000;
        double stressa[ar_len];
        double stressb[ar_len];
        double ts[ar_len];
        int stress_cnt = 0;
        double str_sum_a, str_sum_b;

        int rj;
        int rct;
        double tau_prime;
        double tau_dprime;
        //double sum_mean; // For our system two denominators mu and sig^2 are equal
        double max_xi;
        //double tau_tentative;
        double a0c;
        double* k = malloc(numN4*2*sizeof(double));
        long int* X_temp = malloc((numN4+1)*sizeof(long int));
        long int* X_old = malloc((numN4+1)*sizeof(long int));
        int x_neg = 0;

        int j;

        double tau1, tau2;

        save_param(folder, e0, mu, D, chi, rho, prf, expon, shft, X_T, k_a, k_d, numN, qu_range, du);

        clock_t start_time = clock();
        double elapsed_time;
        int rxn_cnt;

        int start = 0;

        while(t < t_max){

            if(cycles == 50000){start=1;t=0;}

            if(cycles%10000==0&&start){
                printf("%f\n", t);
                printf("Unbound num: %ld\n", X[0]);
                elapsed_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
                printf("Elapsed time %lf\n", elapsed_time);
                start_time = clock();
            }
            if(cycles%100000==0){save_X(t, Ua, Ub, Qa, Qb, X, X_T, numN4, du, folder, C);}
            if(start&&cycles%stress_sample==0){
                if(stress_cnt>=ar_len){
                    save_stress(stress_cnt, ts, stressa, stressb, rpt, folder, flen);
                    stress_cnt = 0;
                }
                else{
                    str_sum_a = 0;
                    str_sum_b = 0;
                    for(i=0; i<numN4; i++){
                        str_sum_a +=  X_T_root*(mu*(Ua[i])+D*Qa[i])*((double) X[i+1])/((double) X_T);
                        str_sum_b +=  X_T_root*(mu*(Ub[i])+D*Qb[i])*((double) X[i+1])/((double) X_T);
                    }
                    ts[stress_cnt] = t;
                    stressa[stress_cnt] = str_sum_a;
                    stressb[stress_cnt] = str_sum_b;
                    stress_cnt++;
                }
            }
            
            /* Gillespie */

            a0 = 0;
            mu_i[0] = 0; // Reset divisors for tprime
            s2_i[0] = 0;

            for(i=0; i<numN4; i++){ // Calculate a=c*x and sum of all a
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


            for(i=0; i<numN4*2; i++){
                rj = V[i];
                // Find population
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
            
            //elapsed_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
            //printf("Elapsed time %lf\n", elapsed_time);
            //start_time = clock();

            /* Find max tau */
            tau_prime = 1000.0;
            /*for(i=0; i<numN2+1; i++){ // numN+1
                max_xi = ((double)X[i])*e_control > 1.0 ? ((double)X[i])*e_control : 1.0;
                if(mu_i[i]<0){mu_i[i] = -mu_i[i];}
                tau_tentative = max_xi/mu_i[i] < (max_xi*max_xi)/s2_i[i] ? max_xi/mu_i[i] : (max_xi*max_xi)/s2_i[i];
                tau_prime = tau_tentative < tau_prime ? tau_tentative : tau_prime;
            }*/

             max_xi = ((double)X[0])*e_control > 1.0 ? ((double)X[0])*e_control : 1.0;
            if(mu_i[0]<0){mu_i[0] = -mu_i[0];} // Abs
            tau1 = max_xi/(mu_i[0]);
            tau2 = (max_xi*max_xi)/s2_i[0];
            tau_prime = tau1 < tau2 ? tau1 : tau2;


            /* Should we do normal Gillespie? */

            //printf("Gil? zeta: %f tauP: %f\n", 1000*zeta/a0, 1000*tau_prime);

            if(tau_prime < zeta/a0){

                printf("What's the fucking point\n");

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
                    }
                    r1 = dsfmt_genrand_open_open(&dsfmt);
                    tau_dprime = log(1/r1)/a0c;
                }
                else{
                    tau_dprime = 10000000;
                }
                
                do{
                    x_neg = 0;

                    tau = tau_prime < tau_dprime ? tau_prime : tau_dprime; // Leap is smaller of two 

                    // k determines how many firings of each reaction in tau leap
                    for(i=0; i<numN4*2; i++){
                            k[i] = 0;
                    }

                    // Create temporary populations to check for negatives -> rejection of step
                    for(i=0; i<numN4+1; i++){
                        X_temp[i] = X[i];
                        if(X[i] < 0){printf("Error: Negative population for %d.\n Old population: %ld New population: %ld\n", i, X_old[i], X[i]);exit(1);}
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

                        crit_cnt = crit[cs_cnt];

                        if(V[crit_cnt] < 0){
                            X_temp[0]++;
                            X_temp[-V[crit_cnt]]--;
                            if(X_temp[-V[crit_cnt]] < 0){printf("Error: Negative population for %d.\n", -V[crit_cnt]);}
                        }
                        else{
                            X_temp[0]--;
                            X_temp[V[crit_cnt]]++;
                        }
                    }
                    rxn_cnt = 0;
                    // For non critical, assign random poisson distributed number of firings with mean/variance a(X)tau
                    for(i=0; i<n_noncrit; i++){
                            k[noncrit[i]] = poisson_Single(tau*a[noncrit[i]]);
                            rxn_cnt += k[noncrit[i]];
                    }

                    //printf("Number of reactions: %d\n", rxn_cnt);
                    
                    int assign_order[n_noncrit];
                    for(i=0; i<n_noncrit; i++){assign_order[i]=i;}

                    shuffle(assign_order, n_noncrit);
                   
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
                for(i=0; i<numN4+1; i++){
                    X_old[i] = X[i];
                    X[i] = X_temp[i];
                }


                if(start){t += tau;}
                cycles++;
            }
        }
        save_stress(stress_cnt, ts, stressa, stressb, rpt, folder, flen); // Final catch of stresses
        save_X(t, Ua, Ub, Qa, Qb, X, X_T, numN4, du, folder, C);

        free(Ua);free(Ub);free(Qa);free(Qb);free(X);free(V);free(C);free(a);
        free(crit);free(noncrit);free(mu_i);free(s2_i);free(X_temp);free(k);
    }
    return 0;
}
