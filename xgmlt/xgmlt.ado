/************************************************************************/
/* xgmlt.ado                                                            */
/*				ESTIMATION OF A LOGIT MODEL FOR X-SECTIONAL DATA			   */
/*								 WITH LOCATION-BASED DEPENDENCE						*/
/*for STATA 6.0                                                         */
/*			by	Robert Vigfusson (email r-vigfusson@nwu.edu)               	*/
/*				Northwestern University                 		               */
/*				October 11, 1999                                            */
/*   based on code 	for x_gmm.ado written by Jean-Pierre Dube	        	*/		
/*                                                                      */
/*			reference:                                     		            */
/*											                                    	*/
/*	Conley, Timothy G.[1999].  "GMM Estimation with Cross		      */
/*	Sectional Dependence." J. of Econometrics 92:1, 1-45.      	*/
/*											                                    	*/
/************************************************************************/
/************************************************************************/


* WARNING: This program is offered without any guarantee. 
* If you find errors, please contact r-vigfusson@nwu.edu.



/*  To invoke this command type:                                        */
/*	>>xgmlt coordlist cutofflist depvar regressorlist,  xreg()  coord()  */
/*																								*/
/*  NOTE: (1) If you want a constant in the regression, specify one of  */
/*	your input variables as a 1. (ie. include it in list of					*/
/*	regressors).																			*/
/*																								*/
/*	(2) MUST specify positive value for xreg() and	coord() options.		*/									*/
/*	(3)	xreg() denotes # regressors												*/
/*	                       															   */
/*		coord()	denotes dimension of coordinates                			*/
/*																								*/
/*																								*/
/*	(5) Your cutofflist must correspond to coordlist (same order)			*/
/* cutofflist is the values of LM. Zero weight is put on terms          */
/* more than LM units apart in that dimmension in the calculation       */
/* of the covariance  term																*/
/*																								*/
/*  OUTPUT:	                                                            */
/*		betalgt= Maximum Likelihood Logit estimator                    	*/
/*		cov_nd= variance-covariance matrix 	of betalgt 							*/
/*							without correction for spatial correlation			*/
/*		cov_dep= variance-covariance matrix 	of betalgt 							*/
/*							without correction for spatial correlation			*/
																						*/
/************************************************************************/

program define xgmlt
version 6.0
syntax varlist(min=1) [, xreg(int -1) COord(int -1)]
/* everything is required. */
#delimit ;				/*sets `;' as end of line*/


/* The following IF statments check to make sure everything is set up correctly. */

if `xreg'<1{;
	if `xreg'==-1{;
		di in red "option xreg() required!!!";
		exit 198};
	di in red "xreg(`xreg') is invalid";
	exit 198};	


if `coord'<1{;
	if `coord'==-1{;
		di in red "option coord() required!!!";
		exit 198};
	di in red "coord(`coord') is invalid";
	exit 198};

/*Separate input variables: coordinates, cutoffs, dependent, regressors */

   
parse "`varlist'", parse(" ");	
local a=1;
while `a'<=`coord'{;
	tempvar coord`a';
	gen `coord`a''=``a'';	/*get coordinates*/
local a=`a'+1};

local aa=1;
while `aa'<=`coord'{;
	tempvar cut`aa';
	gen `cut`aa''=``a'';	/*get cutoffs*/
	local a=`a'+1;
local aa=`aa'+1};

tempvar Y;
gen `Y'=``a'';			/*get dep variable*/
local depend : word `a' of `varlist';

local a=`a'+1;

local b=1;
while `b'<=`xreg'{;
	tempvar X`b';
	local ind`b'="`b'";
	gen `X`b''= ``a'';
	local ind`b' : word `a' of `varlist';
	local a=`a'+1;
local b=`b'+1};			/*get indep variable(s)*/






/*ESTIMATE THE LOGIT Model*/
quietly{

if `xreg'==1{;
		logit `Y' `X1',noconstant};
else{;
		logit `Y' `X1'-`X`xreg'' , noconstant};

/*save the parameter estimates and the variance covariance matrix */
matrix betalgt = e(b);
matrix cov_nd = e(V);

/* Create the residual term */
tempvar phat uhat;
predict `phat'; 			/* Phat is exp(XB)/(1+exp(XB)) */
gen `uhat' = `Y' - `phat';   /* The residual */


};


/*Create the Variance Covariance matrix that corrects for SPATIAL DEPENDENCE*/
tempname XUUX XUUk window XUUX1 XUUX2 XUUXt fix ; 	/* Declare a set of variables */
matrix `XUUX' = J(`xreg',`xreg',0);  /* Intializes the matrix */

quietly{

gen `XUUk'=0;
gen `window'=1;		/*initializes mat.s/var.s to be used*/
local i=1;

while `i'<=_N{;		/*loop through all observations*/
		local j=1;
		replace `window'=1;
		while `j'<=`coord'{;	/*loop through coordinates*/
				if `i'==1{;
						tempvar dis`j';
						gen `dis`j''=0};

				replace `dis`j''=abs(`coord`j''-`coord`j''[`i']);
				replace `window'=`window'*(1-`dis`j''/`cut`j'');
				replace `window'=0 if `dis`j''>=`cut`j'';
				local j=`j'+1};			/*create window*/

		/* End of j loop */

		capture mat drop `XUUX2';
		local k=1;
		while `k'<=`xreg'{;
				replace `XUUk'=`X`k''[`i']*`uhat'*`uhat'[`i']*`window';
				if `xreg'==1{;
						mat vecaccum `XUUX1'=`XUUk' `X1', noconstant};
				else{;
						mat vecaccum `XUUX1'=`XUUk' `X1'-`X`xreg'', noconstant};
				
				mat `XUUX2'=nullmat(`XUUX2') \ `XUUX1';
				local k=`k'+1};
		
		/* End of k loop */

		mat `XUUXt'=`XUUX2'';
		mat `XUUX1'=`XUUX2' +`XUUXt';
		scalar `fix'=.5;	/*to correct for double-counting*/
		mat `XUUX1'=`XUUX1'*`fix';
		mat `XUUX'=`XUUX' +`XUUX1';
		local i=`i'+1};
		
/* End of i loop */

matrix `XUUX'=`XUUX'/_N;

/* We just calculated the Variance Covariance matrix and it is stored in */
/* in the local variable XUUX. */
   
/* The next step is to calculate dg and construct the var-covar matrix for the parameters */
/* for the parameter estimates. */
      

 
/* Compute dg */
/* g = uhat*X so dg = d(uhat)/db * X        */
/* d(uhat)/db(j) = -1*(phat)*(1-phat)*X(j)  */    


/* use phat calculated above using predict  */     
tempvar sqdphat dg ;
gen `sqdphat' = sqrt(`phat'*(1-`phat'));  /* I take the square root to make following calculations easier */
      
local k=1;

while `k'<=`xreg'{;
		tempvar STF`k';
		gen `STF`k''=`X`k''*`sqdphat';
		local k=`k'+1};

/* End of While loop */
		
if `xreg'==1{;
	 	mat accum `dg'=``STF1', noconstant};		
else{;
		mat accum `dg'=`STF1'-`STF`xreg'', noconstant};

 
matrix `dg' = `dg'/_N;           
/* Create the var-covariance matrix. */            
mat cov_dep = inv(`dg'*inv(`XUUX')*`dg'')/_N; 	/*corrected covariance matrix*/


/*end of the quietly section */
};


/*THIS PART CREATES AND PRINTS THE OUTPUT TABLE IN STATA*/
local z=1;
local v=`a';
di _newline(2) _skip(5)
"Results for Logit with Spatial Dependence ";
di _newline	_col(20)	" number of observations=  "  _N;
di _newline "Dependent variable= `depend'";
di _newline
"variable" _col(13) "Coef Est." _col(29) "Standard SE" _col(45) "Spatial SE";
di 
"--------" _col(13) "-------------" _col(29) "---------" _col(45) "------------------";

while `z'<=`xreg'{;
	tempvar beta`z' beta1`z' se`z' see`z' se1`z' se2`z';
	gen `beta`z''=betalgt[1,`z'];
	gen `se`z''=cov_nd[`z',`z'];
	gen `see`z''=sqrt(`se`z'');
	gen `se1`z''=cov_dep[`z',`z'];
	gen `se2`z''=sqrt(`se1`z'');
	di "`ind`z''" _col(13)  `beta`z''  _col(29)   `see`z''  _col(45)  `se2`z'';
local z=`z'+1};
end;

exit;


