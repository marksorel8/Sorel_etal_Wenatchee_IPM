#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_sl );                // number of streams x life histories
  DATA_INTEGER( n_i );                 // number of years x streams x life histories
  DATA_INTEGER( n_t );                 // number of years 
  DATA_INTEGER( n_l );                 // number of life histories 
  DATA_INTEGER( n_s );                 // number of streams 
  DATA_INTEGER( n_f );                 // number of latent variable factors
  
  DATA_IVECTOR( mod );                 // process model vector (e.g. Beverton Holt) for each stream x LH
  DATA_INTEGER( no_rand_FF  );         // flag for whether to use hierarchical structure for functional form parameters
  
  
  DATA_VECTOR(J_obs);                  // Observed log juveniles (posterior mean from juvenile model)
  DATA_VECTOR(J_obs_sd);               // Observed log juveniles standard error (posterior mean from juvenile modle)
  DATA_VECTOR(log_S_obs);              // log mean of observed spawners (log(redd count)-.5 CV^2)
  DATA_SCALAR(S_obs_CV);               // Observed spawners (from redd counts)
  
  DATA_FACTOR(st_i);                   // stream by year index for observation i
  DATA_FACTOR(sl_i);                   // stream and life history index for observation i
  DATA_FACTOR(t_i);                    // year index for observation i
  DATA_SPARSE_MATRIX(X);               // Design matrix of covariates for process error model
  //---------------------------------------------------------------------------------------------------
  
  //Parameters
  PARAMETER_VECTOR(beta_alpha);        // log alpha hyper means (LHP-specific)
  REPORT(beta_alpha); 
  PARAMETER_VECTOR(beta_Jmax);         // log Jmax hyper means (LHP-specific)
  REPORT(beta_Jmax); 
  PARAMETER_VECTOR(beta_gamma);        // log gamma hyper means (LHP-specific)
  REPORT(beta_gamma); 
  
  PARAMETER_VECTOR(log_sigma_alpha);   // log standard deviations of stream-specific log alpha random effects (LHP-specific)
  REPORT(log_sigma_alpha);
  PARAMETER_VECTOR(log_sigma_Jmax);    // log standard deviations of stream-specific log Jmax random effects (LHP-specific)
  REPORT(log_sigma_Jmax);
  PARAMETER_VECTOR(log_sigma_gamma);   // log standard deviations of stream-specific log gamma random effects (LHP-specific)
  REPORT(log_sigma_gamma);
  
  PARAMETER_VECTOR(beta_e);            // environmental variable coefficients           
  PARAMETER_VECTOR(Loadings_vec);      // latent variable factor loadings
  PARAMETER_VECTOR(log_sigma_eta);     // log idiosyncratic process errors standard deviations
  vector<Type> sigma_eta = exp(log_sigma_eta);
  REPORT(sigma_eta);
  PARAMETER_VECTOR(log_sigma_beta_e);
  vector<Type>sigma_beta_e=exp(log_sigma_beta_e);
  
  PARAMETER_VECTOR(log_sigma_Bgamma);
  PARAMETER_VECTOR(log_sigma_BJmax);
  PARAMETER(log_sigma_loading);
  
  PARAMETER_VECTOR(rate);
  
  // Random effects
  PARAMETER_VECTOR(eps_alpha);         // random effect for stream- and LHP-specific alpa
  REPORT(eps_alpha);
  PARAMETER_VECTOR(eps_Jmax);          // random effect for stream- and LHP-specific Jmax
  REPORT(eps_Jmax);
  PARAMETER_VECTOR(eps_gamma);         // random effect for stream- and LHP-specific gamma
  REPORT(eps_gamma);
  PARAMETER_VECTOR(log_S_hat);         // log latent female spawner abundance
  vector<Type> S_hat =exp(log_S_hat);
  REPORT(S_hat);
  PARAMETER_MATRIX(Omega_xf);          // latent factor variables
  REPORT(Omega_xf);
  PARAMETER_VECTOR(eta);               // latent stream by year by life history (idiosyncratic) process errors in juvenile recruitment
  //---------------------------------------------------------------------------------------------------  
  
  //Variables
  
  //// Initialize likelihood
  vector<Type> jnll_comp(3);   
  jnll_comp.setZero();   
  
  DATA_VECTOR(rate_prior);
  jnll_comp(0)  -=dnorm(exp(rate),rate_prior(0),rate_prior(1),true).sum()+rate.sum(); //same penalties as for rates in MS mod
  // jnll-= beta_alpha.sum();
  jnll_comp(0) -=dnorm(beta_alpha,Type((4)),Type(3),true).sum();//beta_alpha.sum();
  
  
  
  
  //// Unpack factor loadings matrix  (taken from J. Thorson spatial class example)
  matrix<Type> Loadings_pf((n_sl), n_f);
  int Count = 0;
  for(int fac=0; fac<n_f; fac++){
    for(int p=0; p<(n_sl); p++){
      if(p==fac){
        Loadings_pf(p,fac) = exp(Loadings_vec(Count));
        Count++;
      }
      else{
        if(p>=fac){
          Loadings_pf(p,fac) = Loadings_vec(Count);//
          // if(p==fac){Loadings_pf(p,fac) = exp(Loadings_pf(p,fac));}
          Count++;
        }else{
          Loadings_pf(p,fac) = 0.0;
        }
      }
    }
  }
  REPORT(Loadings_pf);
  
  // PARAMETER_VECTOR(log_loadings_sd);
  // jnll_comp(0) -= dnorm(Loadings_vec,Type(0),vector<Type>(exp(log_loadings_sd)),true).sum()+
  //   dexp(exp(log_loadings_sd),Type(1),true).sum()+log_loadings_sd.sum();
  
  
  vector<Type> log_alpha(n_sl); // empty vector to hold log alphas
  vector<Type> log_Jmax(n_sl);  // empty vector to hold log Jmaxes
  vector<Type> log_gamma(n_sl); // empty vector to hold log gammas
  
  
  
  Type alpha_bias_correction = 0;
  Type gamma_bias_correction = 0;
  Type Jmax_bias_correction = 0;
  
  
  
  //// calculate alphas, gammas, and Jmaxes
  //////linear predictors on log scale
  for(int j=0; j<n_l; j++){          // loop over life histories
    alpha_bias_correction = pow(exp(log_sigma_alpha(j)),2)/2.0;
    Jmax_bias_correction = pow(exp(log_sigma_Jmax(j)),2)/2.0;
    gamma_bias_correction = pow(exp(log_sigma_gamma(j)),2)/2.0;
    
    for( int i=0; i<n_s; i++){       // loop over streams
      
      jnll_comp(0) -= dnorm( Type(eps_alpha(i*n_l+j)),
                Type(beta_alpha(j) -alpha_bias_correction),
                Type(exp(log_sigma_alpha(j))),true );
      
      jnll_comp(0) -= dnorm( Type(eps_Jmax(i*n_l+j)),
                Type(beta_Jmax(j) -Jmax_bias_correction),
                Type(exp(log_sigma_Jmax(j))),true );
      
      // log_Jmax(i*n_l+j) =            // log maximum recruitment
      //   beta_Jmax(j)* exp(log_sigma_BJmax(j))+log_Jmax_prior_mean +               // life-history intercept
      //   eps_Jmax(i*n_l+j)*exp(log_sigma_Jmax(j))-
      //   Jmax_bias_correction;            // 
      
      
      jnll_comp(0) -= dnorm( Type(eps_gamma(i*n_l+j)),
                Type(beta_gamma(j) -gamma_bias_correction ),
                Type(exp(log_sigma_gamma(j))),true );
      
      // log_gamma(i*n_l+j) =            // log gamma 
      //   beta_gamma(j)* exp(log_sigma_Bgamma(j)) +               // life-history intercept
      //   eps_gamma(i*n_l+j)*exp(log_sigma_gamma(j))-
      //   gamma_bias_correction;
      
    }
  }
  
  
  
  ////// transform to positive real space
  vector<Type> alpha = exp(eps_alpha);
  vector<Type> Jmax = exp(eps_Jmax);
  vector<Type> gamma =  exp(eps_gamma);
  
  
  // calculate latent juvenile emigrants
  vector<Type> J_hat(n_i);      // vector to hold latent juveniles expectations (without process error)
  vector<Type> J_pred(n_i);     // vector to hold latent juveniles (with process error)
  
  //// covariate effects on process errors
  vector<Type> cov_e = X * vector<Type>(beta_e); // design matrix * coefficients
  
  // latent variable factor effects on process errors
  matrix<Type> LV_effects = Loadings_pf* Omega_xf.transpose();
  REPORT(LV_effects);
  
  for( int i=0; i<n_i; i++){      // loop over all observations of streams, life-histories, and years
    
    ///// calculate expected juveniles based on spawners and functional form (without process error)
    J_hat(i)=  (alpha(sl_i(i))* pow(S_hat(st_i(i)),gamma(sl_i(i)))) /
      (Type(1.0)+   alpha(sl_i(i))* pow(S_hat(st_i(i)),gamma(sl_i(i))) / Jmax(sl_i(i)));
    
    jnll_comp(0)-=dnorm(Type(eta(i)),
              Type(log(J_hat(i))+cov_e(i)+ LV_effects(sl_i(i), t_i(i))),
              Type(sigma_eta(sl_i(i))),true);
    
    //// multiplicitive process error
    
    J_pred(i) = (exp(eta(i)));
    
    
  }
  REPORT(J_pred);
  //---------------------------------------------------------------------------------------------------  
  
  // Probabilities of random effects
  // 
  
  //// stream-specific random effects on log alphas, log gammas, and log Jmaxes
  // if(!no_rand_FF){
  // for ( int i = 0; i<n_sl; i++){ // loop over stream x LHP combinations
  //   jnll_comp(0) -= dnorm(eps_alpha(i),Type(0),Type(1),true);
  //   
  //   
  //   // if((mod(i)==1 ) | (mod(i)==3)| (mod(i)==5)){   // only include if gamma is in the model for stream X LHP i
  //   jnll_comp(0) -= dnorm(eps_gamma(i),Type(0),Type(1),true);
  //   
  //   // }      // jnll_comp(0) -= (dexp(exp(log_sigma_gamma(i)),Type(rate(1)),true) + log_sigma_gamma(i));
  //   //   
  //   // }
  //   // if((mod(i)==1) | (mod(i)==2)| (mod(i)==5)){   // only include if Jmax is in the model for stream X LHP i
  //   jnll_comp(0) -= dnorm(eps_Jmax(i),Type(0),Type(1),true);
  //   // }
  // }
  
  
  
  //regulatizing priors on hyper-standard deviations of random effects on log alphas, log gammas, and log Jmaxes
  jnll_comp(0) -= (dexp(exp(log_sigma_alpha),Type(exp(rate(0))),true).sum() + 
    log_sigma_alpha.sum());
  
  jnll_comp(0) -= (dexp(exp(log_sigma_gamma),Type(exp(rate(0))),true).sum() + 
    log_sigma_gamma.sum());
  
  jnll_comp(0) -= (dexp(exp(log_sigma_Jmax),Type(exp(rate(0))),true).sum() +
    log_sigma_Jmax.sum());
  
  
  vector<Type> sigma_Bgamma= exp(log_sigma_Bgamma);
  vector<Type> sigma_BJmax= exp(log_sigma_BJmax);
  vector<Type> bias_correct_Bgamma =(sigma_Bgamma*sigma_Bgamma)/2.0;
  vector<Type> bias_correct_BJmax =(sigma_BJmax*sigma_BJmax)/2.0;
  
  //reluarizing prior on hyper-means for log gammas
  jnll_comp(0) -= dnorm(beta_gamma,0.0-bias_correct_Bgamma, sigma_Bgamma,true).sum();
  jnll_comp(0) -= (dexp(exp(log_sigma_Bgamma),Type(exp(rate(1))),true).sum() +
    log_sigma_Bgamma.sum());
  
  //regularizing prior on hyper-mean of log-Jmaxes
  jnll_comp(0) -= dnorm(beta_Jmax,Type(log(15000))-
    bias_correct_BJmax,
    sigma_BJmax,true).sum();
  jnll_comp(0) -= (dexp(exp(log_sigma_BJmax),Type(exp(rate(2))),true).sum() +
    log_sigma_BJmax.sum());
  
  //regularizing prior on idiosyncratic process error standard deviations
  jnll_comp(0) -= (dexp(exp(log_sigma_eta),Type(exp(rate(3))),true).sum() +
    log_sigma_eta.sum());
  
  
  // if(n_f>0){
  jnll_comp(0) -= dnorm(Type(exp(Loadings_vec(0))),Type(0),Type(exp(log_sigma_loading)),true)+Loadings_vec(0);
  for(int i =1; i<Loadings_vec.size(); i++){
    jnll_comp(0) -= dnorm(Type(Loadings_vec(i)),Type(0),Type(exp(log_sigma_loading)),true);
  }
  jnll_comp(0) -=dexp(exp(log_sigma_loading),Type(exp(rate(3))))+log_sigma_loading;
  // }
  // }
  // PARAMETER_VECTOR(mu_beta_e);
  // PARAMETER_VECTOR(log_sd_beta_e);
  // DATA_IVECTOR(beta_e_ind);
  // for(int i =0; i<beta_e.size();i++){
  //   jnll_comp(0) -=  (dnorm(Type(beta_e(i)), Type(mu_beta_e(beta_e_ind(i))),
  //                     Type(exp(log_sd_beta_e(beta_e_ind(i)))),true));
  //   }
  
  // jnll_comp(0) -= (dexp(vector<Type>(exp(log_sd_beta_e)),Type(1),true).sum()+log_sd_beta_e.sum());
  
  
  
  REPORT(beta_e);  
  
  
  //// latent factor variables
  for(int t = 0; t<(n_t); t++){ // loop over years
    jnll_comp(1)-= dnorm(vector<Type>(Omega_xf.row(t)),Type(0),Type(1),true).sum();
  }
  
  // idiosyncratic process error
  // jnll_comp(1)-=dnorm(eta ,
  //           Type(0.0) ,
  //           //            sigma_eta,
  //           Type(1.0),
  //           true ).sum();
  
  ////  Latent spawners
  jnll_comp(2) -= dnorm(log_S_obs ,log_S_hat  ,S_obs_CV , true ).sum();
  
  //Likelihood
  
  ////  Latent juveniles
  vector<Type> log_J_pred=log(J_pred);
  vector<Type> juv_like= dnorm(J_obs ,log_J_pred ,J_obs_sd , true );
  REPORT(juv_like);
  
  // DATA_IVECTOR(fold_exclude);
  // Type ll_exclude=0;
  // for(int i = 0; i<fold_exclude.size(); i++){
  // ll_exclude +=juv_like(fold_exclude(i));
  // }
  jnll_comp(2) -= juv_like.sum();//+ll_exclude);
  // REPORT(ll_exclude);
  
  // Ridge prior on environmental covariates
  jnll_comp(2) -=(dnorm(beta_e,0,sigma_beta_e,true).sum()+
    dexp(sigma_beta_e,exp(rate(4)),true).sum()+log_sigma_beta_e.sum());
  
  //Return objective function 
  // ADREPORT(beta_e);
  
  // ADREPORT(log_alpha);
  // ADREPORT(log_Jmax);
  // ADREPORT(log_gamma); 
  
  REPORT(gamma);
  REPORT(Jmax);
  REPORT(alpha);
  
  REPORT(jnll_comp);
  Type obj_fun = jnll_comp.sum();
  ADREPORT(log(J_pred));   //get standard errors for log unobserved juveniles
  ADREPORT(log_S_hat);     //get standard errors for log unobserved spaweners
  return(obj_fun); 
}
