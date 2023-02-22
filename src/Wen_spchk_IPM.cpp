#include <TMB.hpp>

//Integrated Population Model for Wenatchee River spring Chinook salmon
// 
// Copyright (C) 2021  Mark Sorel
// 
// This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU Affero General Public License as
//   published by the Free Software Foundation, either version 3 of the
//   License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Affero General Public License for more details.
//   
//   You should have received a copy of the GNU Affero General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//     
//     I can be contacted at marks6@uw.edu or at:
//       Mark Sorel
//       1122 NE Boat Street,
//       Seattle, WA 98105
// 


// Parts of this program were inspired by packages marked by Jeff Laake, Devin Johnson, and Paul Conn 
// <https://cran.r-project.org/web/packages/marked/index.html>  and
// package glmmTMB by Mollie Brooks et al. <https://cran.r-project.org/web/packages/glmmTMB/index.html>
// Code chunks were copied directly from the glmmTMB source code with permission form
// some glmmTMB coauthors. 



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Functions copied from glmmTMB for calculating random effect probabilities
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//valid covariance structures
enum valid_covStruct {
  diag_covstruct = 0, //diagonal
  us_covstruct   = 1, //unstructured
};

//defines elements of list data structure
template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};


//translates r list data structure to C/TMB list data structure. 
//Returns a vector of list, where each  list is for one random effect "component" i.e. inside one parenthesid e.g.(LH|Year)
template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!isNull(t)){
        RObjectTestExpectedType(t, &isNumeric, "times");
        (*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!isNull(d)){
        RObjectTestExpectedType(d, &isMatrix, "dist");
        (*this)(i).dist = asMatrix<Type>(d);
      }
    }
  }
};


//function that calculates the probability of random effects for many different random effects structures.
//Returns negative log prob of random effects for a given random effect compnenet e.g. (LH|year)
template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false,
                  Type pen = 1) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    ans -= (dexp(sd,pen,true).sum() +theta.sum()); //penalize complexity
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      U.col(i)=U.col(i);
      if (do_simulate) {
        U.col(i) = rnorm(Type(0), sd);
      }
    }
    term.sd = sd; // For report
  }
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    ans -= (dexp(sd,pen,true).sum() +logsd.sum()); //penalize complexity
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      // ans += nldens(U.col(i));
      
      // U.col(i)=U.col(i)*sd;
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  
  else error("covStruct not implemented!");
  return ans;
};


//function that creats the structures and call termwise_nll for all random effects. 
//Returns negative log prob of random effects.
template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
                  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false, 
                  Type pen = 0) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate, Type((pen)));
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
};

// data structure that holds the parameter index matrices
template<class Type>
struct pim: vector<matrix<int> > {
  
  pim(SEXP x){ // Constructor
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP m = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<int>(m);
    }
    
  }
};


//----------------------------------------------------
//----------------------------------------------------
//----------------------------------------------------
//Objective funtion

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  //~~~~~~~~~~~~~~~~~~~
  ////Spawner recruit model
  //~~~~~~~~~~~~~~~~~~~
  DATA_INTEGER( n_sl );                // number of streams x life histories
  DATA_INTEGER( n_t );                 // number of years
  DATA_INTEGER( n_l );                 // number of life histories
  DATA_INTEGER( n_s );                 // number of streams
  DATA_INTEGER(n_slt);                  // number of combinations of stream and life history and spawning year modeled
  DATA_INTEGER( n_f );                 // number of latent variable factors
  
  DATA_IVECTOR( mod );                 // process model vector (e.g. Beverton Holt) for each stream x LH

  DATA_VECTOR(J_obs);                  // Observed log juveniles (posterior mean from juvenile model)
  DATA_VECTOR(J_obs_sd);               // Observed log juveniles standard error (posterior mean from juvenile modle)
  DATA_MATRIX(log_S_obs);              // Observed spawners (from redd counts)
  DATA_SPARSE_MATRIX(X);               // Design matrix of covariates for process error model
  DATA_SPARSE_MATRIX(X_proj);               // Design matrix of covariates for process error model
  
  DATA_VECTOR(stream_length);           // length (km) of spawning habitat in each natal stream
  //~~~~~~~~~~~~~~~~~~~
  //// MSCJS model
  //~~~~~~~~~~~~~~~~~~~  
  DATA_INTEGER(n_OCC);         //number of total survival/ recapture occasions;
  DATA_INTEGER(nDS_OCC);       //number of downstream survival/ recapture occasions;
  DATA_INTEGER(nUS_OCC);       //number of upstream survival/ recapture occasions;
  DATA_INTEGER(n_states);      //number of possible adult return ages 
  DATA_INTEGER(n_unique_CH);   //number of unique capture occasions
  DATA_IVECTOR(f);             //release occasion
  
  //CH data
  DATA_IMATRIX(CH);           //capture histories (excluding occasion at marking (which we are conditioning on))
  DATA_IVECTOR(freq);         //frequency of capture histories
  //design matrices fixed effects
  DATA_SPARSE_MATRIX(X_phi);        //fixed effect design matrix for phi
  DATA_SPARSE_MATRIX(X_p);          // fixed effect design matrix for p
  DATA_SPARSE_MATRIX(X_psi);        // fixed effect design matrix for psi
  DATA_SPARSE_MATRIX(X_phi_proj);        //fixed effect design matrix for phi
  DATA_SPARSE_MATRIX(X_psi_proj);        // fixed effect design matrix for psi
  //design matrices random effects
  DATA_SPARSE_MATRIX(Z_phi); //random effect design matrix for phi
  DATA_SPARSE_MATRIX(Z_p);   //random effect design matrix for p
  DATA_SPARSE_MATRIX(Z_psi); //random effect design matrix for psi
  DATA_SPARSE_MATRIX(Z_phi_proj); //random effect design matrix for phi
  DATA_SPARSE_MATRIX(Z_psi_proj); //random effect design matrix for psi
  //PIMS
  DATA_STRUCT(Phi_pim, pim); //index vector of matrices for phi parameter vector for a given Ch x occasion
  DATA_STRUCT(p_pim, pim);   //index vector of matrices for p parameter vector for a given Ch x occasion
  DATA_IVECTOR(Psi_pim);
  // Covariance structures 
  DATA_STRUCT(phi_terms, terms_t);//  Covariance structure for the Phi model
  DATA_STRUCT(p_terms, terms_t);  //  Covariance structure for the p model
  DATA_STRUCT(psi_terms, terms_t);//  Covariance structure for the Psi model
  // for simulation
  DATA_IVECTOR(n_released);   // number of fish released in each cohort (LH x stream x year)
  DATA_IMATRIX(phi_pim_sim);  // index of phi parameters for the simulation 
  DATA_IMATRIX(p_pim_sim);    // index of p parameters for the simulation 
  DATA_IVECTOR(psi_pim_sim);  // index of psi parameters for the simulation 
  DATA_IVECTOR(f_rel);        // occasion of release for each cohort
  //Covariance structures for projection
  DATA_STRUCT(phi_terms_proj, terms_t);//  Covariance structure for the Phi model
  DATA_STRUCT(psi_terms_proj, terms_t);//  Covariance structure for the Psi model
  
  //~~~~~~~~~~~~~~~~~~~
  //// Other
  //~~~~~~~~~~~~~~~~~~~ 
  
  //for calculating fry survival
  DATA_INTEGER(n_sum_0_1);    // number of rows of design matrix corresponding with summer subyearling time 1 survival
  DATA_INTEGER(ind_sum_0_1);  // index of first row of design matrix corresponding with summer subyearling time 1 survival
  DATA_INTEGER(n_sum_0_1_proj);    // number of rows of design matrix corresponding with summer subyearling time 1 survival
  DATA_INTEGER(ind_sum_0_1_proj);  // index of first row of design matrix corresponding with summer subyearling time 1 survival
  DATA_VECTOR(beta_DOY_prior);      // effect of scaled emigration DOY on time one survival for summer subyearling emigrants
  DATA_SCALAR(spr0_DOY);      // difference in mean emigration DOY between spring and summer subyearling emigrants
  // carcass data that infoms pHOS
  DATA_VECTOR(WH_carc);      // number of wild and hatchery carcasses recovered
  DATA_VECTOR(H_carc);       // number of wild carcasses recovered

  
  //~~~~~~~~~~~~~~~~~~~
  //// "Other"
  //~~~~~~~~~~~~~~~~~~~
  DATA_IMATRIX(phi_ind);         //survival index
  DATA_IVECTOR(psi_ind);         // return age index
  DATA_MATRIX(broodstock);       // natural origin broodstock removals
  DATA_IVECTOR(first_t);           // first year modeled by stream
  DATA_INTEGER(last_t);            // last year modeled
  DATA_INTEGER(proj_years);       //number of projection years
  DATA_INTEGER(n_slt_proj);                  // number of combinations of stream and life history and spawning year modeled
  DATA_SCALAR(p_fem_obs_CV);             //observation error on proportion female
  DATA_VECTOR(logit_p_fem_obs);         // observations of proportion female
  DATA_VECTOR(f_5_prior);               // prior on the productivity of 5 year old females relative to four years olds, based on fecundity
  DATA_VECTOR(H_RRS_prior);             // prior on hatchery relative reporductive succss
  DATA_INTEGER(n_ages); //number of adult return ages
  
  DATA_VECTOR(Hmax);       // maximum number of hatchery origin spawners on spawning ground
  DATA_VECTOR(Hslope);     // rate at which HO spawners decreases with increasing NO spawners
  DATA_MATRIX(stray_home); // proprtions of fish from each hatchery program going to each trib
  DATA_VECTOR(bs_prop); //proportion of NO return taken for broodstock
  DATA_VECTOR(BS_Max);  //total number of fish needed for broodstock in each stream 
  
  //penalty on logit average prespawn survival
  DATA_SCALAR(mu_pss_prior_mean);
  DATA_SCALAR(mu_pss_prior_sd);
  // prior on redd observation error standard deviation
  DATA_SCALAR(log_S_obs_CV_sd);
  DATA_SCALAR(log_S_obs_CV_mu);
  
  DATA_SCALAR(rate);  // exponential rate parameters in PC priors
  //---------------------------------------------------------------------------------------------------
  
  //Parameters
  //~~~~~~~~~~~~~~~~~~~
  ////Spawner recruit model
  //~~~~~~~~~~~~~~~~~~
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
  
  PARAMETER_VECTOR(log_sigma_eta);     // log idiosyncratic process errors standard deviations
  vector<Type> sigma_eta = exp(log_sigma_eta);
  REPORT(sigma_eta);
  
  PARAMETER(log_S_obs_CV);
  Type S_obs_CV = exp(log_S_obs_CV);
  
  
  
  //~~~~~~~~~~~~~~~~~~~
  ////MSCJS model
  //~~~~~~~~~~~~~~~~~~~
  //fixed effects
  PARAMETER_VECTOR(beta_phi_ints);  //Phi fixed effect coefficients
  PARAMETER_VECTOR(beta_p_ints);    //p fixed effect coefficients
  PARAMETER_VECTOR(beta_psi_ints);   //psi fixed effect coefficients
  // log standard deviations of PC priors
  PARAMETER_VECTOR(log_pen_sds_phi);   //penalty log SDs
  PARAMETER_VECTOR(log_pen_sds_p);   //penalty log SDs
  PARAMETER_VECTOR(log_pen_sds_psi);   //penalty log SDs
  // Joint vector of covariance parameters
  PARAMETER_VECTOR(theta_phi);  //phi
  PARAMETER_VECTOR(theta_p);    //p
  PARAMETER_VECTOR(theta_psi);  //psi
  PARAMETER(theta_psi_cor);

  
  //~~~~~~~~~~~~~~~~~~~
  ////"Other"  
  //~~~~~~~~~~~~~~~~~~~
  //paramaters
  PARAMETER(beta_DOY);            // effect of scaled DOY on time one survival for spring and summer subyearling emigrants
  PARAMETER_VECTOR(logit_pHOS);   // proportion of hatchery origin spawners
  
  PARAMETER(mu_pss);           // mean prespawn survival on logit scale
  Type PSS =invlogit(mu_pss);
  
  PARAMETER(mu_fem)           // mean logit proportion of females
    PARAMETER(log_p_fem_sd);  // sd logit proportion of females
  Type p_fem_sd = exp(log_p_fem_sd);
  PARAMETER(mu_pHOS);          // mean logit proportion of hatchery origin spawners
  PARAMETER(log_sd_pHOS);      // log sd logit proportion of hatchery origin spawners
 
  PARAMETER_VECTOR(log_sigma_Bgamma);  // hyper-distribution SD for gamma parameters in juvenile production function
  PARAMETER_VECTOR(log_sigma_BJmax);   // hyper-distribution SD for Jmax parameters in juvenile production function
  
  PARAMETER_VECTOR(log_sigma_beta_e);  // log SD of PC prior on environmental covariate effects on juvenile production
  vector<Type>sigma_beta_e=exp(log_sigma_beta_e);
  PARAMETER(log_sigma_loading);        // log SD of PC prior on loadings in juvenile production function
  
  
  PARAMETER(log_f_5_fec);             // productivity of 5 year old spawners relative to 4 year old spawners
  Type f_5_fec=exp(log_f_5_fec);
  
  PARAMETER(logit_RRS);               // relative reproductive success
  Type RRS = invlogit(logit_RRS);
  //---------------------------------------------------------------------------------------------------
  
  //Random effects
  //~~~~~~~~~~~~~~~~~~~
  ////Juvenile production function
  //~~~~~~~~~~~~~~~~~~~
  PARAMETER_VECTOR(eps_alpha);         // random effect for stream- and LHP-specific alpa
  REPORT(eps_alpha);
  PARAMETER_VECTOR(eps_Jmax);          // random effect for stream- and LHP-specific Jmax
  REPORT(eps_Jmax);
  PARAMETER_VECTOR(eps_gamma);         // random effect for stream- and LHP-specific gamma
  REPORT(eps_gamma);
  PARAMETER_MATRIX(Omega_xf);          // latent factor variables
  REPORT(Omega_xf);
  PARAMETER_VECTOR(eta);               // latent stream by year by life history (idiosyncratic) process errors in juvenile recruitment
  
  PARAMETER_VECTOR(beta_e);            // environmental variable coefficients           
  PARAMETER_VECTOR(Loadings_vec);      // latent variable factor loadings
  
  
  
  //~~~~~~~~~~~~~~~~~~~
  ////MSCJS model
  //~~~~~~~~~~~~~~~~~~~
  //penalized effects
  PARAMETER_VECTOR(beta_phi_pen);  //Phi effect coefficients
  PARAMETER_VECTOR(beta_p_pen);    //p effect coefficients
  PARAMETER_VECTOR(beta_psi_pen);   //psi effect coefficients
  //random effects
  PARAMETER_VECTOR(b_phi);      //phi random effects
  PARAMETER_VECTOR(b_p);        //p random effects
  PARAMETER_VECTOR(b_psi);      //psi random effects
  PARAMETER_VECTOR(b_phi_proj);      //phi random effects for projection (not optimized...mapped)
  PARAMETER_VECTOR(b_psi_proj);      //psi random effects for projection (not optimized...mapped)
  
  
  //~~~~~~~~~~~~~~~~~~~
  ////other
  //~~~~~~~~~~~~~~~~~~~
  PARAMETER_VECTOR(eps_p_fem);  // year and stream specific proportion female deviations
  PARAMETER_VECTOR(log_S_init); //  spawner abundance in first 5 years for each stream
  vector<Type>S_init=exp(log_S_init);
  //---------------------------------------------------------------------------------------------------  
  

  
  //// Initialize likelihood
  // Joint negative log-likelihood
  parallel_accumulator<Type> jnll(this);
  
  //concatenate return age marginal standard deviations and correltion 
  vector<Type> theta_psi2(3);
  theta_psi2 << theta_psi,theta_psi_cor;
  
  // multistate model random effect probabilities (allterms_nll returns the nll and also simulates new values of the random effects)
  jnll += allterms_nll(b_phi, theta_phi, phi_terms, false,rate);//);//phi
  jnll += allterms_nll(b_p, theta_p, p_terms, false,rate);//);//p
  jnll += allterms_nll(b_psi, theta_psi2, psi_terms, false,rate);//);//psi

  
  //  Prior on hyper-means of log productivity (alpha) in juvenle production function; ~normal(4,3)
  // jnll-=dnorm(beta_alpha,Type((4)),Type(3),true).sum();
  
  //---------------------------------------------------------------------------------------------------  

  
  //~~~~~~~~~~~~~~~~~~~
  ////Spawner recruit model variables
  //~~~~~~~~~~~~~~~~~~~
  
  //// Unpack recruitment deviation factor loadings matrix  (taken from J. Thorson spatial class example)
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
      
      jnll -= dnorm( Type(eps_alpha(i*n_l+j)),
                Type(beta_alpha(j) -alpha_bias_correction),
                Type(exp(log_sigma_alpha(j))),true );
      
      jnll -= dnorm( Type(eps_Jmax(i*n_l+j)),
                Type(beta_Jmax(j) -Jmax_bias_correction),
                Type(exp(log_sigma_Jmax(j))),true );
       
      
      jnll-= dnorm( Type(eps_gamma(i*n_l+j)),
                Type(beta_gamma(j) -gamma_bias_correction ),
                Type(exp(log_sigma_gamma(j))),true );
      
    }
  }
  
  
  
  ////// transform to positive real space
  vector<Type> alpha = exp(eps_alpha);
  vector<Type> Jmax = exp(eps_Jmax);
  vector<Type> gamma =  exp(eps_gamma);
  
  REPORT(alpha);
  REPORT(Jmax);
  REPORT(gamma);
  
  //// covariate effects on process errors
  vector<Type> cov_e(n_slt+n_slt_proj);
  
  cov_e.head(n_slt) = X * vector<Type>(beta_e); // design matrix * coefficients
  
  SIMULATE { //for population projection (everying for population projection wrapped in SIMULATE so not evaluated during model fitting)
    cov_e.tail(n_slt_proj) = X_proj * vector<Type>(beta_e); // design matrix * coefficients
  }
  REPORT(cov_e);
  
  // latent variable factor effects on process errors
  matrix<Type> LV_effects(12,(Omega_xf.rows()+proj_years));
  ////retrospecitive
  LV_effects.block(0,0,12,Omega_xf.rows()) = matrix<Type>(Loadings_pf * Omega_xf.transpose());
  
  
  //prospective
  SIMULATE {
    matrix<Type> Omega_xf_proj(proj_years,1);
    for(int i=0; i<proj_years; i++) {Omega_xf_proj(i,0)=rnorm(0.0,1.0);}
    LV_effects.block(0,Omega_xf.rows(),12,proj_years) = Loadings_pf * Omega_xf_proj.transpose();
  }
  REPORT(LV_effects);
  
  //independent process error compoenent  
  vector<Type> eta_long(n_slt+n_slt_proj);
  ////retrospective
  eta_long.head(n_slt)=eta;
  ////prospective
  // SIMULATE {
  //   for(int i=n_slt; i<(n_slt_proj+n_slt); i++) {eta_long(i)=rnorm(0.0,1.0);}
  // }
  //~~~~~~~~~~~~~~~~~~~
  ////MSCJS model
  //~~~~~~~~~~~~~~~~~~~
  //allocate vector to concatenate intercepts and penalized coefficient
  vector<Type> beta_phi(beta_phi_ints.size()+beta_phi_pen.size()); 
  /////concatenate intercepts and penalized effects
  beta_phi << beta_phi_ints,beta_phi_pen;
  
  //same process for detection probability
  vector<Type> beta_p(X_p.cols());
  beta_p << beta_p_ints,beta_p_pen;
  
  //and for return are probability
  vector<Type> beta_psi(X_psi.cols());
  beta_psi << beta_psi_ints,beta_psi_pen;
  
  //some reporting
  ADREPORT(beta_phi);
  ADREPORT(beta_p);
  ADREPORT(beta_psi);
  REPORT(beta_phi);
  REPORT(beta_p);
  REPORT(beta_psi);
  
  // Linear predictors for survival (phi) and return age proportions (psi)
  ///allocate vectors to hold
  vector<Type> eta_phi(X_phi.rows()+X_phi_proj.rows()+n_sum_0_1+n_sum_0_1_proj);
  vector<Type> eta_psi(X_psi.rows());
  vector<Type> eta_psi_proj(X_psi_proj.rows());
  //// Fixed component
  eta_phi.head(X_phi.rows()) = X_phi*
    vector<Type>(beta_phi.head(X_phi.cols()));
  vector<Type> eta_p = X_p*beta_p;
  eta_psi = X_psi*beta_psi;
  
  
  SIMULATE {
    ////projection
    eta_phi.segment((X_phi.rows()+n_sum_0_1),X_phi_proj.rows()) = vector<Type>(X_phi_proj*
      vector<Type>(beta_phi.head(X_phi.cols())));
    eta_psi_proj = vector<Type>(X_psi_proj*beta_psi);
  }
  REPORT(eta_psi_proj);
  //// Random component
  eta_phi.head(X_phi.rows()) += Z_phi*b_phi;
  eta_p += Z_p*b_p;
  eta_psi += Z_psi*b_psi;
  
  
  SIMULATE {
    //////projection
    Type nothing=0;
    ///simulate random effects for projection
    nothing += allterms_nll(b_phi_proj, theta_phi, phi_terms_proj, true,rate);//phi
    nothing += allterms_nll(b_psi_proj, theta_psi2, psi_terms_proj, true,rate);//psi
    eta_phi.segment((X_phi.rows()+n_sum_0_1),X_phi_proj.rows()) += Z_phi_proj*b_phi_proj;
    eta_psi_proj += Z_psi_proj*b_psi_proj;
  }
  
  REPORT(b_phi_proj);
  REPORT(b_psi_proj);
  
  
  // calculate spring subyearling emigrant time 1 survival
  eta_phi.segment(X_phi.rows(),n_sum_0_1)=
    vector<Type>(vector<Type>(eta_phi.segment(ind_sum_0_1,n_sum_0_1))+
    Type(beta_DOY*spr0_DOY));
  
  
  SIMULATE {
    eta_phi.tail(n_sum_0_1_proj)=
      vector<Type>(vector<Type>(eta_phi.segment(ind_sum_0_1_proj,n_sum_0_1_proj))+
      Type(beta_DOY*spr0_DOY));
  }
  
  
  
  // Apply link
  vector<Type> phi=invlogit(eta_phi);
  vector<Type> p(eta_p.size()+1);
  p.head(eta_p.size())=invlogit(eta_p);
  p.tail(1)=Type(0);
  REPORT(phi);
  REPORT(p);
  ////psi inverse multinomial logit
  int n_groups = eta_psi.size()/2;
  int n_groups_proj = eta_psi_proj.size()/2;
  matrix<Type> psi(n_groups+n_groups_proj,n_states);
  eta_psi= exp(eta_psi);
  vector<Type> denom = eta_psi.segment(0,n_groups)+eta_psi.segment(n_groups,n_groups)+Type(1);
  psi.block(0,0,n_groups,1)= eta_psi.segment(0,n_groups)/denom;        //return after 1 year
  psi.block(0,1,n_groups,1)= Type(1)/denom;                            //return after 2 year
  psi.block(0,2,n_groups,1)= eta_psi.segment(n_groups,n_groups)/denom; //return after 3 year
  
  
  SIMULATE {
    eta_psi_proj= exp(eta_psi_proj);
    vector<Type> denom_proj = eta_psi_proj.segment(0,n_groups_proj)+
      eta_psi_proj.segment(n_groups_proj,n_groups_proj)+Type(1);
    psi.block(n_groups,0,n_groups_proj,1)= eta_psi_proj.segment(0,n_groups_proj)/denom_proj;        //return after 1 year
    psi.block(n_groups,1,n_groups_proj,1)= Type(1)/denom_proj;                            //return after 2 year
    psi.block(n_groups,2,n_groups_proj,1)= eta_psi_proj.segment(n_groups_proj,n_groups_proj)/denom_proj; //return after 3 year
  }
  
  REPORT(psi); 
  
  //calculate the average return age proportion to use in first 5 years
  vector<Type> temp = exp(beta_psi_ints);
  Type temp2 =  temp.sum()+1;
  temp(0) = 1/temp2;  //proportion of age 4
  temp(1) =  temp(1)/temp2; // proportion of age 5
  Type prop_5_init = temp(1)/(temp(0)+temp(1));
  //~~~~~~~~~~~~~~~~~~~
  ////Other 
  //~~~~~~~~~~~~~~~~~~~
  
  matrix<Type> S_hat(last_t+proj_years,n_s);           //female spawners  stream x year
  matrix<Type> S_eff(last_t+proj_years,n_s);           //effective female spawners stream x year
   matrix<Type> S_A_hat(last_t+proj_years, n_ages*n_s ); //Wild spawners (male + female) at Tumwater Dam stream x year x age
  S_A_hat.setZero();                  //start at 0 because we will be adding to it from different LHs
  
  Type HOS=0; //hatchery origin spawners
  
  matrix<Type> broodstock_proj(proj_years,n_s-1); // projected natural origin broodstock collections
  
  
  matrix<Type> tum_H_ad(last_t+proj_years,n_s);  // number of hatchery adults released upstream of Tumwater dam headed for each stream
  tum_H_ad.setZero();
  matrix<Type> tum_W_ad(last_t+proj_years,n_s);  // number wild adult released upstream of tumwater headed for each strema
  vector<Type> forecast(3);  // forecast (with observation error) of adult returns                  
  
  array<Type> A_tum_y(last_t+proj_years,n_s,n_l,n_ages); //adult returns to TUmwater by year, strea, LH, and age
  

// prior on prespawn survival (assumed constant across years and streams)
  jnll -= dnorm (mu_pss,mu_pss_prior_mean,mu_pss_prior_sd,true);
  
  

  // proportion Female random effects
  vector<Type> eps_p_fem2(eps_p_fem.size()+proj_years*n_s);
  eps_p_fem2.head(eps_p_fem.size())=eps_p_fem;
  SIMULATE{
    for(int i=eps_p_fem.size(); i<((proj_years*n_s)+eps_p_fem.size()); i++) {eps_p_fem2(i)=rnorm(mu_fem,p_fem_sd);}
  }
  
  vector<Type> logit_p_female = eps_p_fem2;
  vector<Type> p_female =invlogit(logit_p_female); // prespawn survival (Tumwater to spawn)
    REPORT(p_female);
  jnll -= dnorm (vector<Type>(logit_p_female.head(eps_p_fem.size())),logit_p_fem_obs,p_fem_obs_CV,true).sum();
  
  // random year and stream effect probabilities
  jnll -= dnorm(eps_p_fem,mu_fem,p_fem_sd,true).sum();

  //// random pHOS probabilities
  jnll -= dnorm(logit_pHOS,mu_pHOS,Type(exp(log_sd_pHOS)),true).sum();
  
  vector<Type> pHOS(logit_pHOS.size()+proj_years*n_s);
  pHOS.head(logit_pHOS.size())=invlogit(logit_pHOS);
  

  
  // prior on spawner observation error CV
  jnll -= dnorm(log_S_obs_CV,log_S_obs_CV_mu,log_S_obs_CV_sd,true);
  
  
  
  vector<Type> J_hat(n_slt+n_slt_proj);          // expected number of juvenile emigrants based on Beverton hold
  vector<Type> J_pred(n_slt+n_slt_proj);         // juvenile emigrant abundance. stream x LH x brood year
  
  vector<Type> JBon(n_slt+n_slt_proj);           // smolts at Bonneville dam across years. stream x LH x brood year
  vector<Type> J_LW(n_slt+n_slt_proj);           // smolts at the lower wenatchee across years. stream x LH x brood year
  matrix<Type> A_tum(n_slt+n_slt_proj, n_ages);  // fish returning as adult to Tumwater stream X LH x brood year, age.
  
  
  Type S_hat_stan = 0;               // spawners per kilometer
  int st =0;                         // counter for index of BY and streams
  int st2 =0;                         // counter for index of BY and streams
  int init_ind =0;                    // counter for initial spawner abundance vector
  int sl = 0;                        // counter for index of stream, LH
  int slt = 0;                       // counter for index of stream, LH, BY
  Type prop_nas = 0;                 // proportion of returning fish from the nason
  Type prop_5 = 0;                   // proportion of fish that are age 5
  Type new_bs_prop=0;
  Type for_chiw_nas = 0;
 
  // ~~~ Begin population process model ~~~~~
  
  for(int t = 0; t<(last_t+proj_years); t++){   //brood year loop
    
    //Tumwater to spawning females
    for(int s =0; s<n_s; s++){    //stream loop
      if (t < first_t(s)) continue; //don't do adult stuff if previous to first year modeled for that stream
      // natural origin female spawners
      if(t<(first_t(s)+5)){
        tum_W_ad(t,s) = S_init(init_ind); //wild returns to tumwater dam in first 5 years
        init_ind++;
      }
      
      //subtract broodstock removals
      if(((t>=(first_t(s)+5))&(t<last_t))){ // HISTORICAL
        tum_W_ad(t,s)= Type(S_A_hat(t,s*n_s)+S_A_hat(t,s*n_s+1)+S_A_hat(t,s*n_s+2)); //add three four and five year olds
        
        for(int a =0; a<n_ages; a++){  //adult ages loop
          
          if((s==0)|((s==1)&(t<20))){ //nason before 2015 (i.e. t ==20) used broodstock exclusively form Nason
            
            tum_W_ad(t,s) -= Type(broodstock(st,a));
            
          }else{
            if(s==1){ // nason in 2015 on takes broodstock in proportion to abundance in nason and Chiwawa
              prop_nas=tum_W_ad(t,s)/(tum_W_ad(t,0)+tum_W_ad(t,1));
              
              tum_W_ad(t,1) -= Type(broodstock(st,a)*prop_nas);//take away from nason
              tum_W_ad(t,0) -= Type(broodstock(st,a)*(1-prop_nas));// take away from Chiwawa
              
            }
          }
        }

        
      }
      st++;
    }
    
    if(t>=last_t){       //// PROJECTION
      for(int s =0; s<n_s; s++){ // stream loop for projection broodstock and hatchery releases
        tum_W_ad(t,s)= Type(S_A_hat(t,s*n_s)+S_A_hat(t,s*n_s+1)+S_A_hat(t,s*n_s+2));
        ////broodstock collection
        if(s==0){
          //Chiwawa

          forecast(s) = tum_W_ad(t,s)*exp(rnorm(0.0,0.10));
          broodstock_proj((t-last_t),s)=bs_prop(s)*forecast(s);
          if(broodstock_proj((t-last_t),s)>BS_Max(s)){
            broodstock_proj((t-last_t),s)=BS_Max(s);
          }
          
          tum_W_ad(t,s) -= broodstock_proj((t-last_t),s);
        }else{
          //Nason
          if(s==1){
            forecast(s) = tum_W_ad(t,s)*exp(rnorm(0.0,0.10));
            new_bs_prop =  (broodstock_proj((t-last_t),0)/forecast(0));
            new_bs_prop = bs_prop(s)-new_bs_prop;
             
            for_chiw_nas = (forecast(s)+forecast(0));
            broodstock_proj((t-last_t),s)=new_bs_prop * for_chiw_nas;
            if(broodstock_proj((t-last_t),s)>BS_Max(s)){
              broodstock_proj((t-last_t),s)=BS_Max(s);
            }
            

            prop_nas=tum_W_ad(t,1)/(tum_W_ad(t,0)+tum_W_ad(t,1));
            tum_W_ad(t,1) -= Type( broodstock_proj((t-last_t),s)*prop_nas);//take away from nason
            tum_W_ad(t,0) -= Type( broodstock_proj((t-last_t),s)*(1-prop_nas));// take away from Chiwawa
            
          }}
        

        
        /////calculate hatchery origin spawners (from Chiwawa and nason programs)
        if(s<2){
          HOS = Hmax(s)-forecast(s)*Hslope(s);
          if(HOS<0.0){
            HOS=0.0;
          }
          
          for(int s_to =0; s_to<n_s; s_to++){ 
            tum_H_ad(t,s_to)+= HOS*stray_home(s,s_to);
          }
        }
        // TODO: add in a step where HOS can't be more than the number of adult HOS which is a funciton of broodstock size in previous years?
        
      }//end of stream loop
    }
    
    
    for(int s =0; s<n_s; s++){ // stream loop adding hatchery and wild and multiplying times PSM and p female
      if (t < (first_t(s))) continue;
      
      if(t<last_t){ // historical 
        
        tum_H_ad(t,s) =  (tum_W_ad(t,s)*pHOS(st2))/Type(1-pHOS(st2));
        S_hat(t,s) =  tum_W_ad(t,s)+tum_H_ad(t,s);      // natural plus hatchery origin spawners
        S_hat(t,s) *= (PSS*p_female(st2));
        // prespawn mortality and percent female
        
        ////  Latent spawners likelihood
        jnll -= dnorm(Type(log_S_obs(t,s)) , 
                      Type(log(S_hat(t,s)))  ,S_obs_CV , true );
        
        if(t<(first_t(s)+5)){
          prop_5 = prop_5_init;
        }else{
          prop_5 = (S_A_hat(t,s*n_s+2)/Type(S_A_hat(t,s*n_s+1)+S_A_hat(t,s*n_s+2)));
        }
        S_eff(t,s) = (tum_W_ad(t,s)*(1.0+(prop_5*f_5_fec)));
        S_eff(t,s) += tum_H_ad(t,s)*RRS;
        S_eff(t,s) *= (PSS*p_female(st2));
        
        
        
        
      }else{// PROJECTION
        
        
        //sum across ages and multiply by prespawn mortality rate and percent females and proportion hatchery
        S_hat(t,s) = tum_W_ad(t,s)+ tum_H_ad(t,s);          // plus hatchery origin spawners
        pHOS(st2) = tum_H_ad(t,s)/S_hat(t,s);
        S_hat(t,s) *= (PSS*p_female(st2));                  // prespawn mortality and percent female
        
        prop_5 = (S_A_hat(t,s*n_s+2)/Type(S_A_hat(t,s*n_s+1)+S_A_hat(t,s*n_s+2)));
        S_eff(t,s) = (tum_W_ad(t,s)*(1.0+(prop_5*f_5_fec)));
        S_eff(t,s) += tum_H_ad(t,s)*RRS;
        S_eff(t,s) *= (PSS*p_female(st2));
        
        
      }
      st2++;
      
    }// end of stream loop
    
    
    
    
    //***************************************************************************************************
    //spawning females to juvenile transition
    
    for(int s =0; s<n_s; s++){    //stream loop
      //don't calculate juveniles if previous to first year modeled for that stream
      // model brood year productions through 2016
      if ((t < first_t(s))|(t>=((last_t+proj_years)-3))) continue;
      S_hat_stan= S_eff(t,s)/stream_length(s);
      for(int l =0; l<n_l; l++){    //juvenile life history strategies
        sl=  s*n_l+l;
        ////calculate number of juvenile emigrants from natal stream
        ///// calculate expected juveniles based on spawners and functional form (without process error)
        /////depensatory Beverton holt II (Myers et al 1995)
        J_hat(slt)=  (alpha(sl)* pow(S_hat_stan,gamma(sl))) /
          (Type(1.0)+   alpha(sl)* pow(S_hat_stan,gamma(sl)) /
            Jmax(sl));
        
        // //// multiplicitive process error
        
        if(slt<(n_slt)){ //idiosyncratic process error
          jnll-=dnorm(Type(eta(slt)),
                      Type(log(J_hat(slt))+(cov_e(slt)+ LV_effects(sl, t))),
                      Type(sigma_eta(sl)),true);
          
        }else{
          SIMULATE {
            eta_long(slt)=rnorm( Type(log(J_hat(slt))+(cov_e(slt)+ LV_effects(sl, t))),
                     Type(sigma_eta(sl)));
          }
        }
  
        
        J_pred(slt) = exp(eta_long(slt));

        
        //***************************************************************************************************
        //juvenile survival
        
        if (t<((last_t+proj_years)-3)){//survival through BY 2015
          ////downstream
          J_LW(slt)  =  J_pred(slt)*stream_length(s);
          J_LW(slt) *= Type(invlogit(eta_phi(phi_ind(slt,0))));
          JBon(slt) = J_LW(slt);
          for(int o=1; o<(nDS_OCC); o++){ //downstream occasions
            ////spawner density effect on survival
              JBon(slt) *= Type(invlogit(eta_phi(phi_ind(slt,o))));
          }
          
          ////Smolt to adult
          ////return age proportion
          A_tum.row(slt) = vector<Type>(vector<Type>(psi.row(psi_ind(slt))) * JBon(slt)*
            Type(invlogit(eta_phi(phi_ind(slt,nDS_OCC)))));
          
          ////upstream survival
          for(int a =0; a<n_ages; a++){  //adult ages a=0 is age 3
            for(int o=0; o<nUS_OCC; o++){ //upstream occasions
              A_tum(slt,a) *= Type(phi(phi_ind(slt,nDS_OCC+1+a*nUS_OCC+o)));
            }// end loop over upstream occasions
            
            if((((t+a+3)<(last_t+proj_years)))){ //don't include jacks in returning spawners, don't bother with fish returning >3 years after the last modeled brood year.
              SIMULATE{
              A_tum_y(t+3+a,s,l,a)= Type(A_tum(slt,a));
              }
              S_A_hat(t+3+a,s*n_s+a)+= Type(A_tum(slt,a)); //add natural origin returns across life histories and translate to index by return years
            
            }
          } //end loop over adult ages
        }//if statement for survival of juvenile to adulthood
        slt++;
      }  // end of loop over life histories
    } //end of loop over streams
  } // end of loop over brood years
  
  //~~~~~~~ end of population process model ~~~~~~~~~~~~~~~~\\
  
  

  REPORT(S_hat);
  REPORT(S_eff);   

  REPORT(tum_W_ad);
  REPORT(tum_H_ad);
  REPORT(JBon);
  REPORT(S_A_hat);
  REPORT(PSS);
  REPORT(pHOS);
  REPORT(eta_phi);
  REPORT(A_tum);
  REPORT(J_pred);
  REPORT(prop_5);
  REPORT(st2)
  REPORT(J_LW);
  REPORT(broodstock_proj);
  REPORT(A_tum_y);

  //---------------------------------------------------------------------------------------------------  
  //---------------------------------------------------------------------------------------------------  
  //---------------------------------------------------------------------------------------------------  
  
  
  // Probabilities of random effects
  
  
  //~~~~~~~~~~~~~~~~~~~
  ////juvenile production function
  //~~~~~~~~~~~~~~~~~~~
  
  // exponential (regularizing) priors on hyperdistribution standard deviations for shape parameters 
  jnll -= (dexp(exp(log_sigma_alpha),rate,true).sum() + 
    log_sigma_alpha.sum());
  
  jnll -= (dexp(exp(log_sigma_gamma),rate,true).sum() + 
    log_sigma_gamma.sum());
  
   jnll  -= (dexp(exp(log_sigma_Jmax),rate,true).sum() +
    log_sigma_Jmax.sum());
  
  vector<Type> sigma_Bgamma= exp(log_sigma_Bgamma);
  vector<Type> sigma_BJmax= exp(log_sigma_BJmax);
  vector<Type> bias_correct_Bgamma =(sigma_Bgamma*sigma_Bgamma)/2.0;
  vector<Type> bias_correct_BJmax =(sigma_BJmax*sigma_BJmax)/2.0;
  
  
  // penaltized complexity prior on gammas (favoring less density dependence)
  jnll -= dnorm(beta_gamma,0.0-bias_correct_Bgamma, sigma_Bgamma,true).sum();
  jnll -= (dexp(exp(log_sigma_Bgamma),rate,true).sum() +
    log_sigma_Bgamma.sum());
  
  // penaltized complexity prior on Jmax (favoring less density dependence)
  jnll -= dnorm(beta_Jmax,Type(log(15000))-
    bias_correct_BJmax,
    sigma_BJmax,true).sum();
  jnll -= (dexp(exp(log_sigma_BJmax),rate,true).sum() +
    log_sigma_BJmax.sum());
  
  // exponential (regularizing) priors on hyperdistribution standard deviations for asynchronous random year effects
  jnll -= (dexp(exp(log_sigma_eta),rate,true).sum() +
    log_sigma_eta.sum());
  
  
  
  
  //// latent factor variables
  for(int t = 0; t<(n_t); t++){ // loop over years
    jnll-= dnorm(vector<Type>(Omega_xf.row(t)),Type(0),Type(1),true).sum();
  }
  


// penalized complexity prior on loadings  
  jnll -= dnorm(Type(exp(Loadings_vec(0))),Type(0),Type(exp(log_sigma_loading)),true)+Loadings_vec(0);
  for(int i =1; i<Loadings_vec.size(); i++){
    jnll -= dnorm(Type(Loadings_vec(i)),Type(0),Type(exp(log_sigma_loading)),true);
  }
  jnll -=dexp(exp(log_sigma_loading),exp(rate))+log_sigma_loading;
  
  
  //Ridge prior on environmental covariates
  jnll-=(dnorm(beta_e,0,sigma_beta_e,true).sum()+
    dexp(sigma_beta_e,rate,true).sum()+log_sigma_beta_e.sum());
  
  
  //~~~~~~~~~~~~~~~~~~~
  ////MSCJS model
  //~~~~~~~~~~~~~~~~~~~
  
  
  // PC priors (Simpson et al 2018)
  //// concatenate all covariates to apply PC priors too
  // Phi
  for (int i =0; i <beta_phi_pen.size(); i++){
    jnll -= (dnorm(beta_phi_pen(i),Type(0),Type(exp(log_pen_sds_phi(i))),true)+
      dexp(exp(log_pen_sds_phi(i)),rate,true)+ 
      log_pen_sds_phi(i)); //for change of variables (log_pen_sds is parameter but penalizing pen_sd)
  }
  
//p
  for (int i =0; i <beta_p_pen.size(); i++){
    jnll -= (dnorm(beta_p_pen(i),Type(0),Type(exp(log_pen_sds_p(i))),true)+
      dexp(exp(log_pen_sds_p(i)),rate,true)+ 
      log_pen_sds_p(i)); //for change of variables (log_pen_sds is parameter but penalizing pen_sd)
  }
  
  //psi
  jnll -= (dnorm(beta_psi_pen,Type(0),vector<Type>(exp(log_pen_sds_psi)),true).sum()+
    dexp(vector<Type>(exp(log_pen_sds_psi)),rate,true).sum()+ 
    log_pen_sds_psi.sum()); //for change of variables (log_pen_sds is parameter but penalizing pen_sd)
  
  //more priors
  jnll -= dnorm(beta_DOY, beta_DOY_prior(0), beta_DOY_prior(1),true);  // prior on effect of emigration doy on spring/summer subyearling survival in interval 1
  jnll -= dnorm(log_f_5_fec, f_5_prior(0), f_5_prior(1),true);     // prior on productivity of 5 year olds relative to 4
  jnll -= dnorm(logit_RRS, H_RRS_prior(0), H_RRS_prior(1),true);  // prior on hatchery relative reporductive success
  
   // regularizing prior for return age probability correlation
  jnll -= dnorm(theta_psi_cor,Type(0.0),Type(0.625),true);
  //---------------------------------------------------------------------------------------------------
  
  //Likelihood
  
  //~~~~~~~~~~~~~~~~~~~
  ////Juvenile emigrant abundance
  //~~~~~~~~~~~~~~~~~~~
  // vector<Type> log_J_pred=log(J_pred.head(n_slt));
  // REPORT(log_J_pred);
  jnll -= dnorm(J_obs , eta  ,J_obs_sd , true ).sum();
  
  //~~~~~~~~~~~~~~~~~~~
  ////hatchery carcasses
  //~~~~~~~~~~~~~~~~~~~
  
  vector<Type> biased_phos =pHOS.head(logit_pHOS.size())*WH_carc;
  jnll -= dpois(H_carc,biased_phos,true).sum();
  //~~~~~~~~~~~~~~~~~~~
  ////MSCJS model
  //~~~~~~~~~~~~~~~~~~~
  ////Variables
  vector<Type> pS(4); //state probs: dead, 1, 2, 3
  Type u = 0;         // holds the sum of probs after each occasion
  Type NLL_it=0;      // holds the NLL for each CH
  Type tmp = 0;       // holds the prob of a given state during observation process in upstream migration
  
  // SIMULATE{
  // vector<Type> NLL_it_vec(n_unique_CH); // holds likelihood of each unique CH
  // }
  
  for(int n=0; n<n_unique_CH; n++){ // loop over individual unique capture histories
    pS.setZero(); //initialize at 0,1,0,0 (conditioning at capture)
    pS(1)=Type(1);
    NLL_it=Type(0); //initialize capture history NLL at 0
    
    //downstream migration
    for(int t=f(n); t<nDS_OCC; t++){       //loop over downstream occasions (excluding capture occasion)
      //survival process
      pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1)); //prob die or stay dead
      pS(1) *= Type(phi(Phi_pim(0)(n,t))); //prob stay alive
      
      //observation process
      pS(1) *= Type(p(p_pim(0)(n,t))*CH(n,t)+ (Type(1)-p(p_pim(0)(n,t)))*(Type(1)-CH(n,t))); //prob observation given alive
      pS(0) *= Type(Type(1)-CH(n,t)); //prob observation given dead
      //acculate NLL
      u = pS.sum();  //sum of probs
      pS = pS/u; //normalize probs
      NLL_it  +=log(u);    //accumulate nll
    }
    
    //ocean occasion
    int t = nDS_OCC;  //set occasion to be ocean occasion
    ////survival process
    pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1)); //prob die or stay dead in ocean
    pS(1) *= Type(phi(Phi_pim(0)(n,t))); //prob survive ocean
    //maturation age process
    pS(2) = pS(1) * psi(Psi_pim(n),1);
    pS(3) = pS(1) * psi(Psi_pim(n),2);
    pS(1) *= psi(Psi_pim(n),0);
    
    
    for(int t=(nDS_OCC+1); t<n_OCC; t++){       //loop over upstream occasions
      
      ////observation process at t-1 (Obs_t below), because I'm going to fix the detection prob at 1 for the last occasion after this loop
      int Obs_t=t-1;
      if(!CH(n,t-1)){
        pS(1) *= Type(Type(1)-p(p_pim(0)(n,Obs_t)));
        pS(2) *= Type(Type(1)-p(p_pim(1)(n,Obs_t)));
        pS(3) *= Type(Type(1)-p(p_pim(2)(n,Obs_t)));
      } else{
        tmp=pS(CH(n,Obs_t))*p(p_pim((CH(n,Obs_t)-1))(n,Obs_t));
        pS.setZero();
        pS(CH(n,Obs_t))=tmp;
      }
      //accumlate NLL
      u = pS.sum();  //sum of probs
      pS = pS/u; //normalize probs
      NLL_it  +=log(u);    //accumulate nll
      //end ocean occasion
      
      //upstream migration
      ////survival process at time t
      pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1))+
        Type((Type(1)-phi(Phi_pim(1)(n,t)))*pS(2))+
        Type((Type(1)-phi(Phi_pim(2)(n,t)))*pS(3));  // sum(prob vec * 1, 1-phi_1, 1-phi_2, 1-phi_3)
      pS(1) *= Type(phi(Phi_pim(0)(n,t)));                 // sum(prob vec * 0,   phi_1,       0,       0)
      pS(2) *=  Type(phi(Phi_pim(1)(n,t)));                 // sum(prob vec * 0,       0,   phi_2,       0)
      pS(3) *=  Type(phi(Phi_pim(2)(n,t)));                 // sum(prob vec * 0,       0,       0,   phi_3)
      
    }
    
    ////observation process at final time assuming detection probability is 1
    if(!CH(n,(n_OCC-1))){
      pS(1) =  Type(0);
      pS(2) =  Type(0);
      pS(3) =  Type(0);
    }else{
      tmp=pS(CH(n,(n_OCC-1)));
      pS.setZero();
      pS(CH(n,(n_OCC-1)))=tmp;
    }
    //accumulate NLL
    u = pS.sum();  //sum of probs
    pS = pS/u; //normalize probs
    NLL_it  +=log(u);    //accumulate nll
    //end observation process at final time
    
    //multiply the NLL of an individual CH by the frequency of that CH and subtract from total jnll
    jnll-=(NLL_it*freq(n));
  }

  //end of likelihood
  
  
  
  //---------------------------------------------------------------------------------------------------
  
  // // Reporting correlation in return ages
  vector<matrix<Type> > corr_psi(psi_terms.size());
  vector<vector<Type> > sd_psi(psi_terms.size());
  for(int i=0; i<psi_terms.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(psi_terms(i).blockNumTheta > 0){
      corr_psi(i) = psi_terms(i).corr;
      sd_psi(i) = psi_terms(i).sd;
    }
  }

  REPORT(corr_psi);
  REPORT(eta_long);

  //return log likelihood
  return(jnll); 
  
  //THe end. Easy ;)
}
