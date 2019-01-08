/*
 * This rundown function for CTH.
 *
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

#ifdef NAN
/* NAN is supported */
#endif

void
CTH_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, double alpha, double eta, 
            double xtrain_to_est_ratio, double propensity)
{   
	Rprintf("Enter CTH_rundown.c");
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    int my_leaf_id;
    pNode otree =  tree;
    pNode otree_tmp = tree;
    pNode tree_tmp = tree;
    
    int opnumber = 0;
    int j, s;
    int tmp_obs, tmp_id;
    double tr_mean, con_mean;
    double tr_sqr_sum, con_sqr_sum;
    double consums, trsums, cons, trs;
    double tr_var, con_var;
   /* double  y_sum, z_sum;
    double yz_sum,  yy_sum, zz_sum;
    
    double k_sum  ; 
    double kz_sum ,  ky_sum , kk_sum ;
    int n; */
    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
	
	 
        cons = 0.;
        trs = 0.;
        consums = 0.;
        trsums = 0.;
        tr_sqr_sum = 0.;
        con_sqr_sum = 0.;
        /*n = 0;
	k_sum = 0.;
	kz_sum = 0.;
	ky_sum = 0.;
	kk_sum = 0.;
	y_sum = 0.;
        yz_sum = 0.;
	yy_sum = 0.;
	z_sum = 0.;
	zz_sum = 0.;*/
	    
	    
        while (cp[i] < tree->complexity) { 
		 Rprintf("cp in CTH_rundown.c is %d.\n", cp);
	        tree = branch(tree, obs);
	        if (tree == 0)
		        goto oops;
	        otree = tree;
	    }
	    xpred[i] = tree->response_est[0];
            my_leaf_id = tree->id;
        
        for (s = k; s < ct.n; s++) {

	    tree_tmp = otree_tmp;
            j = ct.sorts[0][s];
            tmp_obs = (j < 0) ? -(1 + j) : j;
		
		Rprintf(" tree_tmp in CTH_rundown.c is %d.\n",  tree_tmp);
		
		
            while (cp[i] < tree_tmp->complexity) {
                tree_tmp = branch(tree_tmp, tmp_obs);
            }
            tmp_id = tree_tmp->id;
            if (tmp_id == my_leaf_id) {
                if (ct.treatment[tmp_obs] == 0) {
                    cons += ct.wt[tmp_obs];
                    consums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    con_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                } else {
                    trs += ct.wt[tmp_obs];
                    trsums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    tr_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                }
		
		/*n++;
		y_sum += ct.treatment[tmp_obs];
                z_sum += *ct.ydata[tmp_obs];
		yy_sum += ct.treatment[tmp_obs]*ct.treatment[tmp_obs];
                yz_sum += ct.treatment[tmp_obs] * *ct.ydata[tmp_obs];
                zz_sum += *ct.ydata[tmp_obs] * *ct.ydata[tmp_obs];
                k_sum += ct.treatments[tmp_obs];
                ky_sum += ct.treatments[tmp_obs] * ct.treatment[tmp_obs];
                kz_sum += ct.treatments[tmp_obs] * *ct.ydata[tmp_obs];
                kk_sum += ct.treatments[tmp_obs] * ct.treatments[tmp_obs];
		Rprintf(" finish define CTH_rundown.c \n");*/
		
            }
		Rprintf("trs in CTH_rundown.c %d.\n", trs);  
        }

       //if (trs == 0) {
	    Rprintf("if trs in CTH_rundown.c %d.\n", trs);   
          //  tr_mean = tree->parent->xtreatMean[0];
	    tr_mean = 0;
            tr_var = 0;
        /*} else {Rprintf("else trs in CTH_rundown.c %d.\n", trs);
            tr_mean = trsums / trs;
            tree->xtreatMean[0] = tr_mean;
            tr_var = tr_sqr_sum / trs - tr_mean * tr_mean;
		
		Rprintf(" tree->parent->xtreatMean[0] in CTH_rundown.c %d.\n", tree->parent->xtreatMean[0]); 
        }*/
            
	    
	    
       // if (cons == 0) {
	    Rprintf("if cons in CTH_rundown.c %d.\n", cons); 
           // con_mean = tree->parent->xcontrolMean[0];
	    con_mean = 0;
            con_var = 0;
        /*} else {Rprintf("else cons in CTH_rundown.c %d.\n", cons); 
            con_mean = consums / cons;
            tree->xcontrolMean[0] = con_mean;
            con_var = con_sqr_sum / cons - con_mean * con_mean;
        }*/
        
        xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], ct.treatments[obs2], tr_mean, 
                    con_mean, trs, cons, alpha, eta, xtrain_to_est_ratio, propensity);
       /*
       double  beta_1 , beta_0 , beta_2; 
       double var_beta , beta1_sqr_sum ; 
       double beta2_sqr_sum ; 
  
       double tmp;
	    
       beta_1 = ((n* yz_sum *n* yy_sum-n* yz_sum * y_sum * y_sum-y_sum * z_sum *n*kk_sum + y_sum * z_sum * k_sum * k_sum)
	              -(n* kz_sum *n* ky_sum-n* kz_sum * y_sum *k_sum - z_sum * k_sum *n* ky_sum + z_sum * k_sum * k_sum * y_sum)) 
	            / ((n* yy_sum - y_sum * y_sum)*(n* kk_sum - k_sum * k_sum)); 
	        
       beta_2 = ((n* kz_sum *n* kk_sum-n* kz_sum * y_sum * y_sum- z_sum * k_sum *n*yy_sum + z_sum * k_sum * y_sum * y_sum)
	              -(n* yz_sum *n* ky_sum-n* yz_sum * y_sum *k_sum - z_sum * y_sum *n* ky_sum + z_sum * y_sum * y_sum * k_sum)) 
	            / ((n* yy_sum - y_sum * y_sum)*(n* kk_sum - k_sum * k_sum)); 
	        
       beta_0 = (z_sum - beta_1 * y_sum -beta_2 * k_sum) / n;
	        
       double effect = beta_1;
       double effects = beta_2;

       beta1_sqr_sum = beta_1 * beta_1;
       beta2_sqr_sum = beta_2 * beta_2;
       var_beta = beta1_sqr_sum /n- beta1_sqr_sum / (n* n) + beta2_sqr_sum /n- beta2_sqr_sum / (n* n);
       tmp=var_beta;

    xtemp[i] = 4 * ct.max_y * ct.max_y - alpha * (effect*effect+effects*effects)  + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) 
       * (1 - alpha) * tmp; */
	    
	    
	Rprintf("xtemp in CTH_rundown.c %d.\n", xtemp);
    }
    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	    Rprintf("Entered CTH_rundown.c. Double check.\n");
	for (i = 0; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];

	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tr_mean, con_mean);
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    warning("Warning message--see rundown.c");
}
