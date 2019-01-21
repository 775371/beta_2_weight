/*
 * split.Rule = CT
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

static double *sums, *wtsums, *treatment_effect, *treatments_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;
 /*categorical var*/
static double *y_, *z_ , *yz_ ,  *yy_ , *zz_ , *k_ , *kz_ ,  *ky_, *kk_ ;

int
CTinit(int n, double *y[], int maxcat, char **error,
        int *size, int who, double *wt, double *treatment, 
        int bucketnum, int bucketMax, double *train_to_est_ratio)
{ 
    if (who == 1 && maxcat > 0) {
        
        graycode_init0(maxcat);
        countn = (int *) ALLOC(2 * maxcat, sizeof(int));
        tsplit = countn + maxcat;
        treatment_effect = (double *) ALLOC(9 * maxcat, sizeof(double));
        wts = treatment_effect + maxcat;
        trs = wts + maxcat;
        sums = trs + maxcat;
        wtsums = sums + maxcat;
        trsums = wtsums + maxcat;
        wtsqrsums = trsums + maxcat;
        trsqrsums = wtsqrsums + maxcat;
	treatments_effect = trsqrsums + maxcat;
	    
         y_ = (double *) ALLOC(9 * maxcat, sizeof(double));
        z_ = y_ + maxcat;
        yz_ = z_ + maxcat;
        yy_ = yz_ + maxcat;
        zz_ = yy_ + maxcat;
        k_ = zz_ + maxcat;
        kz_ = k_ + maxcat;
        ky_ = kz_ + maxcat;
        kk_ = ky_ + maxcat;
    }
    *size = 1;
    *train_to_est_ratio = n * 1.0 / ct.NumHonest;
   
    return 0;
       
}

void
CTss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, 
     double *risk, double *wt, double *treatment, double *treatments, double max_y,
     double alpha, double eta, double train_to_est_ratio)
	
{   Rprintf("CTss in CT.c start\n");
    int i;
    double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
    double ttreat = 0.;
    double effect;
    double effects;
 
    double tr_var, con_var;
    double con_sqr_sum = 0., tr_sqr_sum = 0.;
    double var_beta = 0., beta1_sqr_sum = 0.; /* var */
    double  y_sum = 0., z_sum = 0.;
    double yz_sum = 0.,  yy_sum = 0., zz_sum = 0.;
    
    double k_sum =0. ; /* two beta*/
    double kz_sum = 0.,  ky_sum = 0., kk_sum = 0.;
    
    double  beta_1 = 0., beta_0 = 0., beta_2=0.;    
    double beta2_sqr_sum = 0.; /* var */  
   
 
    for (i = 0; i < n; i++) {
        temp1 += *y[i] * wt[i] * treatment[i];
        temp0 += *y[i] * wt[i] * (1 - treatment[i]);
        twt += wt[i];
        ttreat += wt[i] * treatment[i];
        tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
        con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
        
        y_sum += treatment[i];
        z_sum += *y[i];   
        yz_sum += *y[i] * treatment[i];
       
        yy_sum += treatment[i] * treatment[i];
        zz_sum += *y[i] * *y[i];
        k_sum+= treatments[i];
        kk_sum += treatments[i] * treatments[i];
        ky_sum+= treatments[i] * treatment[i];
        kz_sum+= *y[i] * treatments[i];
            
    }

   
    //effect = temp1 / ttreat - temp0 / (twt - ttreat);        
    tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
    con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));
   
   
     /* Y= beta_0 + beta_1 T_1+beta_2 T_2 */
    beta_1 = (
            (twt* yz_sum *twt* kk_sum - twt* yz_sum * k_sum * k_sum - y_sum * z_sum *twt* kk_sum + y_sum * z_sum * k_sum * k_sum)
            -(twt* kz_sum *twt* ky_sum-twt* kz_sum * y_sum * k_sum - z_sum * k_sum *twt* ky_sum + z_sum * k_sum * k_sum * y_sum)) 
            / ( (twt * yy_sum - y_sum * y_sum) * (twt* kk_sum - k_sum * k_sum) - (twt * ky_sum - yy_sum * kk_sum)); 
        
    beta_2 = ((twt* kz_sum *twt* yy_sum-twt* kz_sum * y_sum * y_sum- z_sum * k_sum *twt*yy_sum + z_sum * k_sum * y_sum * y_sum)
              -(twt* yz_sum *twt* ky_sum -twt* yz_sum * y_sum *k_sum - z_sum * y_sum *twt* ky_sum + z_sum * y_sum * y_sum * k_sum)) 
            / ((twt* yy_sum - y_sum * y_sum)*(twt* kk_sum - k_sum * k_sum)-(twt*ky_sum-yy_sum*kk_sum) ); 
        
    beta_0 = (z_sum - beta_1 * y_sum -beta_2 * k_sum) / twt;
        
    effect = beta_1;
    effects=beta_2;
    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;
    var_beta = eta*(beta1_sqr_sum /twt- beta_1 * beta_1 / (twt* twt)) + (1-eta)*(beta2_sqr_sum /twt- beta_2 * beta_2 / (twt* twt));
    
    *tr_mean= temp1 / ttreat;
    *con_mean= temp0 / (twt - ttreat);
    *value = effect;

    //*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect + 
    //(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));
    *risk = 4 * twt * max_y * max_y - alpha * twt * (eta*effect*effect+(1-eta)*effects*effects) 
            + (1 - alpha) * (1 + train_to_est_ratio) * twt * ( var_beta);
    Rprintf("twt in CTss in CT.c %d.\n", twt);   
 }


void CT(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, 
        int *csplit, double myrisk, double *wt, double *treatment, double *treatments,  int minsize, double alpha, double eta,
        double train_to_est_ratio)
{   
	Rprintf("CT in CT.c start\n");
    int i, j;
    double temp;
    double temps;
 
 
    double left_sum, right_sum;
    double left_tr_sum, right_tr_sum;
    double left_tr, right_tr;
    double left_wt, right_wt;
    int left_n, right_n;
    double best;
    int direction = LEFT;
    int where = 0;
    double node_effect, left_effect, right_effect;
    double left_temp, right_temp;
    double left_temps, right_temps;
 
    int min_node_size = minsize;
    
    double tr_var, con_var;
    double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
    double left_tr_var, left_con_var, right_tr_var, right_con_var;

    right_wt = 0.;
    right_tr = 0.;
    right_sum = 0.;
    right_tr_sum = 0.;
    right_sqr_sum = 0.;
    right_tr_sqr_sum = 0.;
    right_n = n;
	
    double   right_y_sum = 0., right_z_sum = 0.;
    double  left_y_sum = 0., left_z_sum = 0.;
    double right_yz_sum = 0.,  right_yy_sum = 0., right_zz_sum = 0.;
    double left_yz_sum = 0.,  left_yy_sum = 0., left_zz_sum = 0.;
    double  beta_1 = 0., beta_0 = 0.; 
    double   beta1_sqr_sum = 0.,  var_beta = 0.; /* beta*/
    double   beta2_sqr_sum = 0.; /* beta*/ 
    double left_k_sum =0. ; /* two beta*/
    double left_kz_sum = 0.,  left_ky_sum = 0., left_kk_sum = 0.;
    double right_k_sum =0. ; /* two beta*/
    double right_kz_sum = 0.,  right_ky_sum = 0., right_kk_sum = 0.;
    double  beta_2=0.;     
   
    for (i = 0; i < n; i++) {
        right_wt += wt[i];
        right_tr += wt[i] * treatment[i];
        right_sum += *y[i] * wt[i];
        right_tr_sum += *y[i] * wt[i] * treatment[i];
        right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
        right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
      
       
        right_y_sum += treatment[i];
        right_z_sum += *y[i];
        right_yz_sum += *y[i] * treatment[i];
       
        right_yy_sum += treatment[i] * treatment[i];
        right_zz_sum += *y[i] * *y[i];
        right_k_sum+= treatments[i];
        right_kk_sum += treatments[i] * treatments[i];
        right_ky_sum+= treatments[i] * treatment[i];
        right_kz_sum+= *y[i] * treatments[i];
    }

    beta_1 = ((right_wt * right_yz_sum *right_wt* right_yy_sum- right_wt * right_yz_sum * right_y_sum * right_y_sum-right_y_sum * right_z_sum *right_wt *right_kk_sum + right_y_sum * right_z_sum * right_k_sum * right_k_sum)
              -(right_wt * right_kz_sum * right_wt * right_ky_sum-right_wt * right_kz_sum * right_y_sum *right_k_sum - right_z_sum * right_k_sum * right_wt * right_ky_sum + right_z_sum * right_k_sum * right_k_sum * right_y_sum)) 
            / ((right_wt * right_yy_sum - right_y_sum * right_y_sum)*(right_wt * right_kk_sum - right_k_sum * right_k_sum) 
               -(right_wt*right_ky_sum-right_yy_sum*right_kk_sum)); 
        
    beta_2 =  ((right_wt * right_kz_sum *right_wt* right_kk_sum- right_wt * right_kz_sum * right_y_sum * right_y_sum- right_z_sum * right_k_sum *right_wt *right_yy_sum + right_z_sum * right_k_sum * right_y_sum * right_y_sum)
              -(right_wt * right_yz_sum * right_wt * right_ky_sum-right_wt * right_yz_sum * right_y_sum *right_k_sum - right_z_sum * right_y_sum * right_wt * right_ky_sum + right_z_sum * right_y_sum * right_y_sum * right_k_sum)) 
            / ((right_wt * right_yy_sum - right_y_sum * right_y_sum)*(right_wt * right_kk_sum - right_k_sum * right_k_sum)
               - (right_wt*right_ky_sum-right_yy_sum*right_kk_sum));
        
    beta_0 = (right_z_sum - beta_1 * right_y_sum -beta_2 * right_k_sum) / right_wt;
        
    temp = beta_1;
    temps=beta_2;
    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;
    var_beta = eta*(beta1_sqr_sum / right_wt - beta_1 * beta_1 / (right_wt * right_wt) )+ 
            (1-eta)*(beta2_sqr_sum / right_wt - beta_2 * beta_2 / (right_wt * right_wt));
    
    Rprintf("beta_1 in CTss in CT.c %d.\n", beta_1); 
    Rprintf("beta_2 in CTss in CT.c %d.\n", beta_2); 
  
        
   /* beta_1 = (right_n * right_yz_sum - right_z_sum * right_y_sum) / (right_n * right_yy_sum - right_y_sum * right_y_sum);
    beta_0 = (right_z_sum - beta_1 * right_y_sum) / right_n;
    temp = beta_1;
    beta_sqr_sum = beta_1 * beta_1 ;
    var_beta = beta_sqr_sum / n - beta_1 * beta_1 / (n * n);*/
    //temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
    tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
    con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
        - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
        / ((right_wt - right_tr) * (right_wt - right_tr));
   /* node_effect = alpha * temp * temp * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
        * right_wt * (tr_var / right_tr  + con_var / (right_wt - right_tr));*/
   
    node_effect = alpha * (eta*temp*temp+(1-eta)*temps*temps)  * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
        * right_wt * (var_beta);
 
 
 
    Rprintf("nclass in CT in CT.c %d.\n", nclass); 
 
 
 
    if (nclass == 0) {
        /* continuous predictor */
        left_wt = 0;
        left_tr = 0;
        left_n = 0;
        left_sum = 0;
        left_tr_sum = 0;
        left_sqr_sum = 0;
        left_tr_sqr_sum = 0;
        best = 0;
        
        for (i = 0; right_n > edge; i++) {
            
         //   Rprintf("The wt[i] in function CT in CT.c is %d\n",(int)wt[i]);
        //    Rprintf("The treatment[i] in function CT in CT.c is %d\n",(int)treatment[i]);
         //   Rprintf("The y[i] in function CT in CT.c is %d\n",(int)y[i]);
       
            left_wt += wt[i];
            right_wt -= wt[i];
            left_tr += wt[i] * treatment[i];
            right_tr -= wt[i] * treatment[i];
            left_n++;
            right_n--;
            temp = *y[i] * wt[i] * treatment[i];
            left_tr_sum += temp;
            right_tr_sum -= temp;
            left_sum += *y[i] * wt[i];
            right_sum -= *y[i] * wt[i];
            temp = (*y[i]) *  (*y[i]) * wt[i];
            left_sqr_sum += temp;
            right_sqr_sum -= temp;
            temp = (*y[i]) * (*y[i]) * wt[i] * treatment[i];
            left_tr_sqr_sum += temp;
            right_tr_sqr_sum -= temp;
                
           
            left_y_sum += treatment[i];
            right_y_sum -= treatment[i];
            left_z_sum += *y[i];
            right_z_sum -= *y[i];
            left_yz_sum += *y[i] * treatment[i];
            right_yz_sum -= *y[i] * treatment[i];
           
            left_yy_sum += treatment[i] * treatment[i];
            right_yy_sum -= treatment[i] * treatment[i];
            left_zz_sum += *y[i] * *y[i];
            right_zz_sum -= *y[i] * *y[i];
              /* add treatments */  
             left_k_sum += treatments[i];
             right_k_sum -= treatments[i];
             left_ky_sum += *y[i] * treatments[i];
             right_ky_sum -= *y[i] * treatments[i];
           
            left_kk_sum += treatments[i] * treatments[i];
            right_kk_sum -= treatments[i] * treatments[i];
            left_kz_sum += treatments[i] * *y[i];
            right_kz_sum -= treatments[i] * *y[i]; 
                
            
           /* if (x[i + 1] != x[i] && left_n >= edge &&
                (int) left_tr >= min_node_size &&
                (int) left_wt - (int) left_tr >= min_node_size &&
                (int) right_tr >= min_node_size &&
                (int) right_wt - (int) right_tr >= min_node_size) {*/                      
                                            
    
                    
            if (x[i + 1] != x[i] &&
                (int) left_wt >= min_node_size &&
                (int) right_wt  >= min_node_size) {
                    
                    
                    
     
     beta_1 = ((left_wt * left_yz_sum *left_wt* left_yy_sum- left_wt * left_yz_sum * left_y_sum * left_y_sum-left_y_sum * left_z_sum *left_wt *left_kk_sum + left_y_sum * left_z_sum * left_k_sum * left_k_sum)
              -(left_wt * left_kz_sum * left_wt * left_ky_sum-left_wt * left_kz_sum * left_y_sum *left_k_sum - left_z_sum * left_k_sum * left_wt * left_ky_sum + left_z_sum * left_k_sum * left_k_sum * left_y_sum)) 
            / ((left_wt * left_yy_sum - left_y_sum * left_y_sum)*(left_wt * left_kk_sum - left_k_sum * left_k_sum)- (left_wt*left_ky_sum-left_yy_sum*left_kk_sum)); 
        
    beta_2 =  ((left_wt * left_kz_sum *left_wt* left_kk_sum- left_wt * left_kz_sum * left_y_sum * left_y_sum- left_z_sum * left_k_sum *left_wt *left_yy_sum + left_z_sum * left_k_sum * left_y_sum * left_y_sum)
              -(left_wt * left_yz_sum * left_wt * left_ky_sum-left_wt * left_yz_sum * left_y_sum *left_k_sum - left_z_sum * left_y_sum * left_wt * left_ky_sum + left_z_sum * left_y_sum * left_y_sum * left_k_sum)) 
            / ((left_wt * left_yy_sum - left_y_sum * left_y_sum)*(left_wt * left_kk_sum - left_k_sum * left_k_sum)- (left_wt*left_ky_sum-left_yy_sum*left_kk_sum));
        
    beta_0 = (left_z_sum - beta_1 * left_y_sum -beta_2 * left_k_sum) / left_wt;
        
    left_temp = beta_1;
       
    left_temps = beta_2; 
                    
    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;
    var_beta = eta*(beta1_sqr_sum / left_wt - beta_1 * beta_1 / (left_wt * left_wt) ) + (1-eta)*(beta2_sqr_sum / left_wt - beta_2 * beta_2 / (left_wt * left_wt));

     
    left_effect =  alpha *(eta*left_temp*left_temp+(1-eta)*left_temps*left_temps) * left_wt - (1 - alpha) * (1 + train_to_est_ratio) 
                    * left_wt * (var_beta);

                   
      //left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
                /*left_tr_var = left_tr_sqr_sum / left_tr - 
                    left_tr_sum  * left_tr_sum / (left_tr * left_tr);
                left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
                    - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
                    / ((left_wt - left_tr) * (left_wt - left_tr));        
                left_effect = alpha * left_temp * left_temp * left_wt
                        - (1 - alpha) * (1 + train_to_est_ratio) * left_wt 
                    * (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));*/         

    
     beta_1 = ((right_wt * right_yz_sum *right_wt* right_yy_sum- right_wt * right_yz_sum * right_y_sum * right_y_sum-right_y_sum * right_z_sum *right_wt *right_kk_sum + right_y_sum * right_z_sum * right_k_sum * right_k_sum)
              -(right_wt * right_kz_sum * right_wt * right_ky_sum-right_wt * right_kz_sum * right_y_sum *right_k_sum - right_z_sum * right_k_sum * right_wt * right_ky_sum + right_z_sum * right_k_sum * right_k_sum * right_y_sum)) 
            / ((right_wt * right_yy_sum - right_y_sum * right_y_sum)*(right_wt * right_kk_sum - right_k_sum * right_k_sum)- (right_wt*right_ky_sum-right_yy_sum*right_kk_sum)); 
        
    beta_2 =  ((right_wt * right_kz_sum *right_wt* right_kk_sum- right_wt * right_kz_sum * right_y_sum * right_y_sum- right_z_sum * right_k_sum *right_wt *right_yy_sum + right_z_sum * right_k_sum * right_y_sum * right_y_sum)
              -(right_wt * right_yz_sum * right_wt * right_ky_sum-right_wt * right_yz_sum * right_y_sum *right_k_sum - right_z_sum * right_y_sum * right_wt * right_ky_sum + right_z_sum * right_y_sum * right_y_sum * right_k_sum)) 
            / ((right_wt * right_yy_sum - right_y_sum * right_y_sum)*(right_wt * right_kk_sum - right_k_sum * right_k_sum)- (right_wt*right_ky_sum-right_yy_sum*right_kk_sum));
        
    beta_0 = (right_z_sum - beta_1 * right_y_sum -beta_2 * right_k_sum) / right_wt;
        
    right_temp = beta_1 ;
    right_temps = beta_2 ;
                    
    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;
    var_beta = eta * (beta1_sqr_sum / right_wt - beta_1 * beta_1 / (right_wt * right_wt)) + (1-eta) * (beta2_sqr_sum / right_wt - beta_2 * beta_2 / (right_wt * right_wt));

    right_effect = alpha * (eta*right_temp*right_temp+(1-eta)*right_temps*right_temps)  * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
                    * right_wt * (var_beta);
                    
/* right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
                /*right_tr_var = right_tr_sqr_sum / right_tr -
                    right_tr_sum * right_tr_sum / (right_tr * right_tr);
                right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                    - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
                    / ((right_wt - right_tr) * (right_wt - right_tr));
                right_effect = alpha * right_temp * right_temp * right_wt
                        - (1 - alpha) * (1 + train_to_est_ratio) * right_wt * 
                            (right_tr_var / right_tr + right_con_var / (right_wt - right_tr));*/
                
                temp = left_effect + right_effect - node_effect;
                //Rprintf("The temp in function CT in CT.c is %d\n",temp);
                     Rprintf("temp in cont in CT.c %d.\n", temp); 
		    Rprintf("best in cont in CT.c %d.\n", best); 
                if (temp > best) {Rprintf("cont: compare temp and best\n");
                    best = temp;
                    where = i;     
				  Rprintf("best after in cont %d.\n", best);
                    if (left_temp < right_temp){
                        direction = LEFT;
                    }
                    else{
                        direction = RIGHT;
                    }
                }             
            }
        }
        
        *improve = best;
	    Rprintf("improve in cont in CT.c %d.\n", *improve); 
	    
        if (best > 0) {         /* found something */
        csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2; 
        }
    }
    
    /*
    * Categorical predictor
    */
    else {
	    Rprintf("nclass in CT in CT.c %d.\n", nclass); 
        for (i = 0; i < nclass; i++) {
            countn[i] = 0;
            wts[i] = 0;
            trs[i] = 0;
            sums[i] = 0;
            wtsums[i] = 0;
            trsums[i] = 0;
            wtsqrsums[i] = 0;
            trsqrsums[i] = 0;
          
        y_[i] =  0;
        z_ [i]=  0;
        yz_[i] =  0;
        yy_[i] =  0;
        zz_ [i]=  0;
        k_ [i]=  0;
        kz_[i] =  0;
        ky_ [i]=  0;
        kk_ [i]=  0;
                
        }
        
        /* rank the classes by treatment effect */
        for (i = 0; i < n; i++) {
            j = (int) x[i] - 1;
            countn[j]++;
            wts[j] += wt[i];
            trs[j] += wt[i] * treatment[i];
            sums[j] += *y[i];
            wtsums[j] += *y[i] * wt[i];
            trsums[j] += *y[i] * wt[i] * treatment[i];
            wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
            trsqrsums[j] +=  (*y[i]) * (*y[i]) * wt[i] * treatment[i];

        y_[j] += treatment[i];
        z_[j] += *y[i];   
        yz_[j] += *y[i] * treatment[i];
        yy_[j] += treatment[i] * treatment[i];
        zz_[j] += *y[i] * *y[i];
        k_[j]+= treatments[i];
        kk_[j] += treatments[i] * treatments[i];
        ky_[j]+= treatments[i] * treatment[i];
        kz_[j]+= *y[i] * treatments[i];

                
        }
        
        for (i = 0; i < nclass; i++) {Rprintf("nclass in function CT in CT.c is %d\n", nclass);
				      
            if (countn[i] > 0) { 
				
				
                tsplit[i] = RIGHT;
               //treatment_effect[i] = trsums[j] / trs[j] - (wtsums[j] - trsums[j]) / (wts[j] - trs[j]);
				
            Rprintf("two treatment_effect start in CT.c \n");    
				
     treatment_effect[i] =  ((wts[i]* yz_[i]*wts[i]* kk_[i] - wts[i]* yz_[i] * k_[i] * k_[i] -
      y_[i] * z_[i] *wts[i]* kk_[i] + y_[i] * z_[i] * k_[i] * k_[i])-(wts[i]* kz_[i] *wts[i]* ky_[i]-
      wts[i]* kz_[i] * y_[i] * k_[i] - z_[i] * k_[i] *wts[i]* ky_[i] + z_[i] * k_[i] * k_[i] * y_[i])) 
      /( (wts[i] * yy_[i] - y_[i] * y_[i]) * (wts[i]* kk_[i]- k_[i] * k_[i]) - (wts[i] * ky_[i] - yy_[i] * kk_[i])); 

	Rprintf("treatment_effect[i] in function CT in CT.c is %d\n", treatment_effect[i]); 
		  //  treatments_effect[i]=0;
		    
treatments_effect[i] =  ((wts[i]* kz_[i] *wts[i]* yy_[i]-wts[i]* kz_[i] * y_[i] * y_[i]- 
z_[i] * k_[i] *wts[i]*yy_[i] + z_[i] * k_[i] * y_[i] * y_[i]) -(wts[i]* yz_[i] *wts[i]* ky_[i] -
wts[i]* yz_[i] * y_[i] *k_[i] - z_[i] * y_[i] *wts[i]* ky_[i] + z_[i] * y_[i] * y_[i] * k_[i])) 
 /( (wts[i]* yy_[i] - y_[i] * y_[i]) * (wts[i]* kk_[i] - k_[i] * k_[i]) - ( wts[i] * ky_[i]-yy_[i] * kk_[i]) );
  
Rprintf("treatments_effect[i] in function CT in CT.c is %d\n", treatments_effect[i]);
                 
            } else
                tsplit[i] = 0;
        }
        graycode_init2(nclass, countn, treatment_effect, treatments_effect);
        
        /*
         * Now find the split that we want
         */
        
        left_wt = 0;
        left_tr = 0;
        left_n = 0;
        left_sum = 0;
        left_tr_sum = 0;
        left_sqr_sum = 0.;
        left_tr_sqr_sum = 0.;
        
        best = 0;
        where = 0;
        while ((j = graycode()) < nclass) {
            tsplit[j] = LEFT;
            left_n += countn[j];
            right_n -= countn[j];
            
            left_wt += wts[j];
            right_wt -= wts[j];
            
            left_tr += trs[j];
            right_tr -= trs[j];
            
            left_sum += wtsums[j];
            right_sum -= wtsums[j];
            
            left_tr_sum += trsums[j];
            right_tr_sum -= trsums[j];
            
            left_sqr_sum += wtsqrsums[j];
            right_sqr_sum -= wtsqrsums[j];
            
            left_tr_sqr_sum += trsqrsums[j];
            right_tr_sqr_sum -= trsqrsums[j];

           left_y_sum += y_[j];
            right_y_sum -= y_[j];
            left_z_sum += z_[j];
            right_z_sum -= z_[j];
            left_yz_sum += yz_[j];
            right_yz_sum -= yz_[j];
           
            left_yy_sum += yy_[j];
            right_yy_sum -= yy_[j];
            left_zz_sum += zz_[j];
            right_zz_sum -= zz_[j];
              /* add treatments */  
             left_k_sum += k_[j];
             right_k_sum -= k_[j];
             left_ky_sum += ky_[j];
             right_ky_sum -= ky_[j];
           
            left_kk_sum += kk_[j];
            right_kk_sum -= kk_[j];
            left_kz_sum += kz_[j];
            right_kz_sum -= kz_[j];


		
		
            
            if (left_n >= edge && right_n >= edge &&
              
                (int) left_wt  >= min_node_size &&
                
                (int) right_wt  >= min_node_size) {
	
                /*left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) 
                    / (left_wt - left_tr);
                
                left_tr_var = left_tr_sqr_sum / left_tr 
                    - left_tr_sum  * left_tr_sum / (left_tr * left_tr);
                left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
                    - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
                    / ((left_wt - left_tr) * (left_wt - left_tr));       
                left_effect = alpha * left_temp * left_temp * left_wt
                    - (1 - alpha) * (1 + train_to_est_ratio) * left_wt * 
                        (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));
                
                right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) 
                    / (right_wt - right_tr);
                right_tr_var = right_tr_sqr_sum / right_tr 
                    - right_tr_sum * right_tr_sum / (right_tr * right_tr);
                right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                    - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
                    / ((right_wt - right_tr) * (right_wt - right_tr));
                right_effect = alpha * right_temp * right_temp * right_wt
                        - (1 - alpha) * (1 + train_to_est_ratio) * right_wt *
                            (right_tr_var / right_tr + right_con_var / (right_wt - right_tr));*/
                    
                    /*treatment split change*/
    beta_1 = ((left_n * left_yz_sum *left_n* left_yy_sum- left_n * left_yz_sum * left_y_sum * left_y_sum-left_y_sum * left_z_sum *left_n *left_kk_sum + left_y_sum * left_z_sum * left_k_sum * left_k_sum)
              -(left_n * left_kz_sum * left_n * left_ky_sum-left_n * left_kz_sum * left_y_sum *left_k_sum - left_z_sum * left_k_sum * left_n * left_ky_sum + left_z_sum * left_k_sum * left_k_sum * left_y_sum)) 
            / ((left_n * left_yy_sum - left_y_sum * left_y_sum)*(left_n * left_kk_sum - left_k_sum * left_k_sum)); 
        
    beta_2 =  ((left_n * left_kz_sum *left_n* left_kk_sum- left_n * left_kz_sum * left_y_sum * left_y_sum- left_z_sum * left_k_sum *left_n *left_yy_sum + left_z_sum * left_k_sum * left_y_sum * left_y_sum)
              -(left_n * left_yz_sum * left_n * left_ky_sum-left_n * left_yz_sum * left_y_sum *left_k_sum - left_z_sum * left_y_sum * left_n * left_ky_sum + left_z_sum * left_y_sum * left_y_sum * left_k_sum)) 
            / ((left_n * left_yy_sum - left_y_sum * left_y_sum)*(left_n * left_kk_sum - left_k_sum * left_k_sum));
        
    beta_0 = (left_z_sum - beta_1 * left_y_sum - beta_2 * left_k_sum) / left_n;

    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;
    left_temp = eta * beta1_sqr_sum+ (1-eta)* beta2_sqr_sum;
    
    var_beta = eta * (beta1_sqr_sum / n - beta_1 * beta_1 / (n * n)) + (1-eta)*(beta2_sqr_sum / n - beta_2 * beta_2 / (n * n));
    left_effect = alpha * (eta * beta1_sqr_sum+ (1-eta)* beta2_sqr_sum) * left_wt - (1 - alpha) * (1 + train_to_est_ratio) 
                    * left_wt * (var_beta);
                    
                    
    beta_1 = ((right_n * right_yz_sum *right_n* right_yy_sum- right_n * right_yz_sum * right_y_sum * right_y_sum-right_y_sum * right_z_sum *right_n *right_kk_sum + right_y_sum * right_z_sum * right_k_sum * right_k_sum)
              -(right_n * right_kz_sum * right_n * right_ky_sum-right_n * right_kz_sum * right_y_sum *right_k_sum - right_z_sum * right_k_sum * right_n * right_ky_sum + right_z_sum * right_k_sum * right_k_sum * right_y_sum)) 
            / ((right_n * right_yy_sum - right_y_sum * right_y_sum)*(right_n * right_kk_sum - right_k_sum * right_k_sum)); 
        
    beta_2 =  ((right_n * right_kz_sum *right_n* right_kk_sum- right_n * right_kz_sum * right_y_sum * right_y_sum- right_z_sum * right_k_sum *right_n *right_yy_sum + right_z_sum * right_k_sum * right_y_sum * right_y_sum)
              -(right_n * right_yz_sum * right_n * right_ky_sum-right_n * right_yz_sum * right_y_sum *right_k_sum - right_z_sum * right_y_sum * right_n * right_ky_sum + right_z_sum * right_y_sum * right_y_sum * right_k_sum)) 
            / ((right_n * right_yy_sum - right_y_sum * right_y_sum)*(right_n * right_kk_sum - right_k_sum * right_k_sum));
        
    beta_0 = (right_z_sum - beta_1 * right_y_sum - beta_2 * right_k_sum) / right_n;

    beta1_sqr_sum = beta_1 * beta_1;
    beta2_sqr_sum = beta_2 * beta_2;    
    right_temp = eta * beta_1 * beta_1 + (1-eta) * beta_2 * beta_2;
    
    var_beta = eta*(beta1_sqr_sum / n - beta_1 * beta_1 / (n * n) ) + (1-eta)*(beta2_sqr_sum / n - beta_2 * beta_2 / (n * n));
    
   
    right_effect = alpha * (eta * beta_1 * beta_1 + (1-eta) * beta_2 * beta_2) * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
                    * right_wt * (var_beta);
                    
             
    temp = left_effect + right_effect - node_effect;
         
    Rprintf("temp in cat in CT.c %d.\n", temp); 
    Rprintf("best in cat in CT.c %d.\n", best); 
    
    if (temp > best) {
		    Rprintf("YES!cat: compare temp and best\n");
                    best = temp;
				  
                    Rprintf("best after in cat is %d\n", best);
				  
                    if (left_temp > right_temp)
                        for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
                    else
                        for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
		}
	    }
	}

        *improve = best;
	    Rprintf("improve cat is %d\n", *improve);
    }
        Rprintf("End function CT in CT.c \n");
} /*CT FUNCTION*/


double
    CTpred(double *y, double wt, double treatment, double *yhat, double propensity)
    {
        double ystar;
        double temp;
        
        ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
        temp = ystar - *yhat;
        return temp * temp * wt;
    } 
