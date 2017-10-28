c From aagiunt@sahp4947.sandia.gov Mon Dec  1 11:12:15 2003
c Received: from waffle.cs.vt.edu (waffle.cs.vt.edu [128.173.40.51])
c 	by cayuga.cs.vt.edu (8.12.1/8.12.1) with ESMTP id hB1GCEpn012225
c	for <ltw@cayuga.cs.vt.edu>; Mon, 1 Dec 2003 11:12:14 -0500 (EST)
c Received: from mm02snlnto.son.sandia.gov (mm02snlnto.sandia.gov [132.175.109.21])
c 	by waffle.cs.vt.edu (8.12.10/8.12.10) with ESMTP id hB1GC55g030646
c 	for <ltw@cs.vt.edu>; Mon, 1 Dec 2003 11:12:06 -0500
c Received: from 132.175.109.1 by mm02snlnto.son.sandia.gov with ESMTP (
c  Tumbleweed MMS SMTP Relay 01 (MMS v5.5.3)); Mon, 01 Dec 2003 09:11:55
c  -0600
c Received: from sahp4947.sandia.gov (sahp4947.sandia.gov [134.253.159.14]
c  ) by sass165.sandia.gov (8.12.10/8.12.10) with ESMTP id hB1GBrTC004807;
c  Mon, 1 Dec 2003 09:11:53 -0700 (MST)
c Received: (from aagiunt@localhost) by sahp4947.sandia.gov (8.11.6/8.11.6
c  ) id hB1GBqs05915; Mon, 1 Dec 2003 09:11:52 -0700
c Date: Mon, 1 Dec 2003 09:11:52 -0700
c From: "Anthony A. Giunta" <aagiunt@sandia.gov>
c Message-ID: <200312011611.hB1GBqs05915@sahp4947.sandia.gov>
c To: ltw@cs.vt.edu
c Subject: file: mars36_fort.f
c cc: aagiunt@sandia.gov
c X-PMX-Version: 4.1.0.80455
c MIME-Version: 1.0
c X-WSS-ID: 13D5B7C1236725-01-01
c Content-Type: text/plain
c Content-Transfer-Encoding: 7bit
c Status: R
c 
c 
c file: mars36_fort.f
c ----------------------
c 
c
c Multivariate Adaptive Regression Splines (MARS modeling, version 3.6).
c
c
c Coded and copywrite (c) by Jerome H. Friedman (3/25/93).
c
c
c
c                         A Micro User's Guide
c                                  to
c                               MARS 3.6
c
c                          Jerome H. Friedman
c                          Stanford University
c
c    MARS 3.6 is a collection of subroutines that implement the multivariate
c adaptive regression spline strategy for data fitting and function
c approximation described in Friedman (1991a, 1991b, 1993). It is a general-
c ization of MARS 1.0 described in Friedman (1988) and MARS 2.5 described in
c Friedman (1991a). It implements as a special case a palindromically
c invariant version of the TURBO fitting technique for smoothing and additive
c modeling described in Friedman and Silverman (1989).
c
c    These subroutines represent a set of tools that can be invoked from a user
c coded program to perform various analyses. The user routine is responsible for
c reading the data into memory and passing it to the MARS subroutines, along
c with the various parameter settings, as arguments. This set of subroutines 
c can also form the basis for incorporating this methodology into a
c statistical language or package.
c
c    The user interface subroutines are:
c
c MARS: computes the mars model from the data and provides summary
c     information for interpreting it.
c
c MISS: sets up mars to handle missing values.
c
c PLOT: constructs graphical output, useful for interpreting the
c     continuous part of the mars model, in a format that can be plotted
c     with a local graphics package.
c
c CATPRT: prints tables useful for interpreting the purely categorical parts
c     of the mars model.
c
c SLICE: produces user selected lower dimensional representations of the
c    mars model to aid in interpretation.
c
c FMOD: computes response estimates for given predictor variable vectors
c    given the mars model.
c
c   It should be noted that this is a new methodology with which there is,
c as of this time, little collective experience.
c
c
c References:
c
c [1] Friedman, J. H. (1988). Fitting functions to noisy data in high
c     dimensions. Proc., Twentyth Symposium on the Interface, Wegman, Gantz,
c     and Miller, eds. American Statistical Association, Alexandria, VA. 3-43.
c
c [2] Friedman, J. H. (1991a). Multivariate adaptive regression splines
c     (with discussion).  Annals of Statistics, 19, 1-141 (March).
c
c [3] Friedman, J. H. (1991b). Estimating functions of mixed ordinal and
c     categorical variables using adaptive splines. Department of Statistics,
c     Stanford University, Tech. Report LCS108.
c
c [4] Friedman, J. H. (1993). Fast MARS. Department of Statistics,
c     Stanford University, Tech. Report LCS110.
c
c [5] Friedman, J. H. and Silverman, B. W. (1989). Flexible parsimonious
c     smoothing and additive modeling (with discussion). TECHNOMETRICS, 31,
c     3-39 (Feburary).
c
c
c
c User interface subroutines:
c
c       All arguments in the calling sequence of user called MARS
c       subroutines must be of the same type in the calling program
c       as indicated in the subroutine's documentation below. All
c       calling sequence arrays must be dimensioned  in the calling
c       program as indicated in the documentation below. This includes
c       all workspace arrays. The leading dimensions of multidimensional
c       arrays must match, whereas the last dimension and singley 
c       dimensioned arrays can be as large or larger than that indicated.
c
c       More detailed explanations for some of the quanities below can
c       be found in the indicated sections of references [2], [3], and/or
c       [4] above.
c
c
c
c call mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm):
c
c input:
c n = number of observations.
c p = number of predictor variables per observation (integer).
c x(n,p) = predictor variable data matrix.
c y(n) = response value for each observation.
c w(n) = weight (mass) for each observation.
c nk = maximum number of basis functions.(Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c mi = maximum number of variables per basis function (interaction level).
c      mi=1 => additive modeling (main effects only);
c      mi>1 => up to mi-variable interactions allowed.
c lx(p) = predictor variable flags; lx(i) corresponds to the ith variable:
c    lx(i) =  0 : exclude variable from model.
c             1 : ordinal variable - no restriction.
c             2 : ordinal variable that can only enter additively;
c                 no interactions with other variables.
c             3 : ordinal variable that can enter only linearly.
c            -1 : categorical variable - no restriction.
c            -2 : categorical variable that can only enter additively;
c                 no interactions with other variables.
c
c output:
c fm(3+nk*(5*mi+nmcv+6)+2*p+ntcv), im(21+nk*(3*mi+8)) = mars model.
c    (nmcv = maximum number of distinct values for any categorical variable;
c     ntcv = total number of distinct values over all categorical variables.)
c
c    note: upon return im(1) and im(2) contain the lengths of the fm and im
c          arrays (respectively) actually used by the program.
c
c workspace:
c sp(n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk) : real.
c dp(max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)) : double precision.
c mm(n*p+2*max(mi,nmcv)) : integer.
c
c
c defaults:
c    the following quanities are set to default values in mars. their values
c    can be changed by executing any of the following statements before the
c    call to mars, with the new value as the argument.
c
c call speed(is):
c is = speed acceleration factor (1-5).
c    larger values progressively sacrifice optimization thuroughness for
c    computational speed advantage. this usually results in marked decrease
c    in computing time with little or no effect on resulting approximation
c    accuracy (especially useful for exploratory work).
c    is = 1 => no acceleration.
c    is = 5 => maximum speed advantage.
c    (default: is=4) (Ref [4] Secs. 3.0 - 4.0)
c
c call logit(il):
c il = ordinary/logistic regression flag.
c    il=0 => ordinary least-squares regression.
c    il=1 => logistic regression.
c    (default: il=0). If logistic regression is selected (il=1) then
c    the response variable is assumed to take on only the values 0/1 and
c    the mars model is for the log-odds: f(x) = log (Pr(Y=1:x)/Pr(Y=0:x)).
c    (Ref[2] Sec. 4.5)
c
c call setdf(df):
c df = number of degrees-of-freedom charged for (unrestricted)
c      knot optimization. (default: df=3.0)
c      (Ref[2] Sec. 3.6, Ref[3] Sec. 2.3)
c
c call xvalid(ix):
c ix = control parameter for sample reuse technique used to automatically
c    estimate smoothing parameter df (see above) from the data.
c    ix = 0 => no effect (default). value used for df is set by user if
c            setdf(df) is called (see above), otherwise default value 
c            (df=3.0) is used.
c    ix > 0 => ix - fold cross-validation.
c    ix < 0 => single validation pass using every (-ix)th (randomly selected)
c            observation as an independent test set.
c    if ix.ne.0 then call setdf(df) (see above) has no effect. if ix > 0,
c    computation increases roughly by a factor of ix over that for ix = 0.
c    for ix < 0 computation increases approximately by a factor of two.
c    (Ref[3] Sec. 2.3)
c 
c call stseed(is):
c is = seed for internal random number generator used to group observation
c      subsets for validation (ix.ne.0). (default: is=987654321).
c
c call print(it):
c it = fortran file number for printed output. (it.le.0 => no printed output.)
c      note that this controls printed output for all user called mars
c      routines. (default: it=6).
c
c call setfv(fv):
c fv = (fractional) incremental penalty for increasing the number of variables
c      in the mars model. sometimes useful with highly collinear designs as it
c      may produce nearly equivalent models with fewer predictor variables,
c      aiding in interpretation. (fv .ge. 0)
c    fv=0.0  => no penalty (default).
c    fv=0.05 => moderate penalty.
c    fv=0.1  => heavy penality.
c    the best value depends on the specific situation and some user
c    experimentation using different values is usually required. this option
c    should be used with some care. (Ref[2] Sec. 5.3)
c
c call setic(ic):
c ic = flag restricting categorical - ordinal interactions.
c    ic=0 => no effect (default).
c    ic=1 => interactions between categorical and ordinal variables prohibited.
c    ic=2 => maximum number of ordinal variables participating in any
c            interaction is restricted to two. categorical interactions are
c            unrestricted.
c    the restrictions associated with a value of ic are imposed in addition
c    to those that are controlled by the mi and lx flags (see above).
c
c call setint(i,j,k):
c i,j = predictor variable numbers.
c   k = interaction flag:
c      k.eq.0 => interactions between variables i and j are prohibited.
c      k.ne.0 => interactions between variables i and j are permitted
c                if allowed by the mi, lx, and ic parameter values (see above).
c                (default)
c   i.eq.0 .or. j.eq.0 => reset to defaults for all predictor variables.
c
c   this call can be executed repeatedly to invoke (or remove) multiple
c   constraints.
c
c call nest (n,i,j,nv,vals);
c nests predictor variable i within categorical predictor variable j.
c links the existance of a value for var(i) to a subset of values for var(j).
c variable j must be categorical (lx(j) < 0, see above). (Ref[3] Sec. 3.3)
c
c n = same as in mars (see above).
c i = variable to be nested (ordinal or categorical).
c j = variable to which i is nested (categorical).
c nv = size of corresponding subset of values for variable j.
c vals(nv) = specific values of variable j for which values of variable i
c            are defined to exist.
c
c setting nv = 0 removes the nesting associated with i and j previously set.
c setting i = 0 and/or j = 0 removes all nesting (default).
c
c this call can be executed repeatedly to invoke (or remove) nestings.
c (recursive nesting is not implemented. j may not itself be nested to
c another variable.)
c
c note: variable nesting does NOT override the interaction constraints set by
c       calling setic or setint (see above). the mi and lx flags (see above)
c       do not limit interactions on variables to which others are nested.
c       they control nested and other nonnested variables as described above
c       except that interactions with variables to which others are nested
c       are ignored in applying the constraints.
c
c
c call setms(ms): 
c ms = minimum span (minimum number of observations between each knot).
c      ms .le. 0 => default value (depending on n and p) is used.
c      (default: ms=0). (Ref[2] Sec. 3.8)
c
c
c the following three routines (miss, mkmiss, xmiss) are used to enable
c mars to deal with various aspects of missing predictor values in the
c training data and/or future data to be predicted. for problems with
c no such missing values, these three routines are of no concern.
c
c
c call miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp);
c
c called  (optionally) before mars (see above) to indicate the presence of
c missing values in the predictor variable data matrix.
c
c sets up mars to accomodate missing values.
c produces as output transformations of some of the original mars input
c quanities defining the problem. the new transformed quanities define
c the same problem - but with missing values - and replace the corresponding
c original quanities in the call to mars (see above) and plot and slice
c (see below). in particular, a new predictor data matrix is created
c with extra (dummy) variables - one for each original variable with
c missing values - to indicate observations for which the corresponding
c original variable is not missing. each such original variable is
c automatically nested (see above) to this (corresponding) dummy variable.
c also produces as output a slicing vector to be used as input to
c slice (see below) to produce the mars model corresponding to no missing
c values on any of the original predictor variables, for interpretation.
c (Ref[3] Sec. 3.4)
c
c input:
c n,p,x,lx = original intended input to mars (see above).
c xm(p) = vector giving missing value flag for each original variable.
c         x(i,j) = xm(j) => x(i,j) is missing.
c flg = input to slice (see below).
c
c output:
c pn,xn,lxn = corresponding transformed quanities to be used as input to mars.
c pn = number of predictor variables in transformed data matrix (integer,
c      .le. 2*p)
c xn(n,pn) = transformed data matrix.
c lxn(pn) = predictor variable flags for transformed data matrix.
c xs(pn) = input for slice (see below). slicing vector that produces a slice
c          of the augmented predictor variable space to produce the nonmissing
c          mars model.
c xp(2*p+1) = input for xmiss (see below).
c
c notes:
c    (1) the value of the output quanity pn is less than or equal to 2*p.
c        the other output quanities (arrays) should be dimensioned
c        xn(n,2*p), lxn(2*p), xs(2*p) in the calling program to be safe.
c    (2) if there are no missing values then the output quanities 
c        pn,xn,lxn will be identical to p,x,lx respectively, and
c        xs(j)=flg, j=1,p.
c    (3) the corresponding quanities for input and output can be the same in
c        the calling program. however, they should have the corresponding
c        dimension equal to 2*p, and they will be altered.
c    (4) dimensions of all relevant workspace arrays in mars and other
c        user callable routines must be large enough to accomodate the
c        the increased number of variables pn (.le. 2*p) produced.
c
c
c call mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc);
c
c called (optionally) before miss (see above) to generate additional data
c with missing values by resampling from original (input) data.
c
c used to train mars for missing predictor values when they are under
c represented in the training data. takes as input original data and
c produces as output a new (larger) data set containing both the original
c data and additional data resampled from it with predictor variable values
c flagged as missing. (Ref[3] Sec. 3.4)
c
c input:
c n,p,x,y,w = original data set (see mars above).
c xm(p) = same as in miss (see above).
c pm(p) = vector of fractions of missing values in output sample:
c    pm(j) = fraction of the total output data with jth predictor variable
c            value missing.
c nnx = maximum sample size of new (output) data.
c
c output:
c nn = sample size of generated output data set (input to miss, mars, and
c      other user called routines, in place of n).
c xnn(nn,p) = new generated predictor data matrix to be used as input to
c             miss (see above) in place of x.
c yn(nn),wn(nn) = new generated data to be used as input to mars and other
c                 user called routines in place of y and w.
c
c workspace:
c sc(p*nnx): real.
c
c notes:
c    (1) the output arrays should be dimensioned xnn(nnx,p),yn(nnx),wn(nnx)
c        in the calling program to be safe.
c    (2) if much larger fractions of missing values are requested than
c        exist in the orginal training data then the size of the output
c        (resampled) data set can become very large (nn = big).
c    (3) if the size of the output (resampled) data set reaches the
c        specified maximum (nnx) then the actual fractions of missing
c        values for each variable will be less than those specified in the
c        input vector pm(p).
c    (4) dimensions of all relevant workspace arrays in mars and other
c        user callable routines must be large enough to accomodate the
c        the increased number of observations nn (.le. nnx) produced.
c
c
c call xmiss (n,x,xm,xp,xn);
c
c must be called before fmod (see below) if miss (see above) was called
c before mars (see above) to handle missing values. produces a new set
c of covariate vectors, from the original set intended for fmod, to replace
c the original set in the call to fmod. (Ref[3] Sec. 3.4)
c
c input:
c n = number of covariate vectors (see fmod below).
c x(n,p) = original covariate vectors intended for fmod.
c xm, xp = same as in miss (see above).
c
c output:
c xn(n,pn) = new set of covariate vectors to be used as input to fmod.
c
c notes:
c    (1) the value of the quanity pn is less than or equal to 2*p (p = the
c        number of original predictor variables). the output quanity (array)
c        should be dimensioned xn(n,2*p) in the calling program to be safe.
c    (2) the corresponding quanities x and xn can be the same in
c        the calling program. however, it should have the second
c        dimension equal to 2*p, and it will be altered.
c
c
c
c
c The following subroutines can be called only after mars.
c
c
c
c call plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm):
c
c computes plots of all purely additive and all purely bivariate ordinal
c contributions to the mars model,in a form suitable for displaying
c with a computer graphics package. If there are no interactions between
c categorical and ordinal variables present in the mars model, this
c subroutine returns plots of the (compound) anova functions (Ref[2] Sec. 3.5,
c eqns (25) and (27)) for those ordinal variables involved in at most two
c (ordinal) variable interactions. If categorical-ordinal interactions are
c present, it returns the corresponding plots for each ordinal contribution
c in the categorical-ordinal decomposition (Ref[3] Sec. 3.2).  since the
c locations of the plotted functions are arbitrary, they are all translated
c to have zero minimum value.
c
c input:
c   m = model flag:
c     = 1 => plot piecewise-linear mars model.
c     = 2 => plot piecewise-cubic mars model. (Ref[2] Sec. 3.7)
c   x,fm,im = same as in mars (see above).
c   ngc = number of raster points for computing curve estimates.
c   ngs = number of raster points on each axis for computing surface estimates.
c   icx = convex hull flag:
c       = 0 => plot surface estimates over entire range of argument limits.
c       > 0 => plot surface estimates only inside the convex hull of the
c             (bivariate) point set.
c
c output:
c    nc = number of curves (purely additive ordinal anova functions).
c    crv(ngc,2,nc) = additive ordinal anova functions in anova
c                    decomposition order.
c       crv(.,1,m) = ordered abscissa values for mth anova function.
c       crv(.,2,m) = corresponding ordinate values.
c    ns = number of surfaces (purely bivariate ordinal anova functions).
c    srf(ngs,ngs,ns) = bivariate (plus associated univariate) ordinal
c                      anova functions in anova decomposition order.
c       srf(.,.,m) = total contribution (bivariate + univariate) of ordinal
c                  variables associated with mth bivariate anova function of the
c                  mars model, in a square raster format (ngs x ngs) over the
c                  ranges of the two variables. the first and second indicies
c                  correspond to the first and second variables respectively.
c
c workspace:
c    sp(max(4*ngs*ngs,ngc,2*n)) : real.
c    mm(max(2*(mi+1),nmcv)) : integer.
c
c
c
c call catprt (m,fm,im,sp,mm):
c
c prints all univariate and bivariate contributions to the purely categorical
c part of the mars model that are not involved in higher order (categorical) 
c interactions. These are the (compound) anova functions (Ref[2] Sec. 3.5,
c eqns (25) and (27)) for the purely categorical part of the categorical-
c ordinal decomposition (Ref[3] Sec. 3.2, eqns (32b) and (38b)). since
c the locations of the printed functions are arbitrary they are all
c translated to have zero minimum value. the functions are printed as tables
c of integers in the range [0,99] or [0,9]. the function value for each
c entry is the product of the corresponding integer and the scale factor
c printed above the table.
c
c input:
c m = model flag:
c   = 1 => piecewise-linear model.
c   = 2 => piecewise-cubic  model. (Ref[2] Sec. 3.7)
c fm,im = same as in mars (see above).
c
c workspace:
c sp(nmcv*nmcv) : real.
c mm(2*nk+nmcv) : integer.
c
c
c
c call slice (flg,xs,x,fm,im,fmn,imn,sp,mm):
c
c computes the mars model within a lower dimensional subspace defined by an
c axis oriented slice of the predictor variable space. the slice is selected
c by assigning specific values to a subset of the predictor variables. the
c returned model is a function of the variables complement to the selected
c subset, and represents the mars model conditioned on the specified values
c for the selected variables. this new lower dimensional sliced model can
c be input to plot and/or catprt for interpretation in the same manner as the
c original mars model. (Ref[2] Sec. 4.7, Ref[3] Sec. 2.4)
c
c input:
c flg = flag for indicating that a predictor variable is not in the subset
c    defining the slice. its value should be outside the range of values for
c    all predictor variables.
c xs(p) = vector defining the slice:
c    xs(i).eq.flg => do not condition on ith variable.
c    xs(i).ne.flg => condition on ith variable at value stored in xs(i).
c x = same as in mars (see above).
c fm,im = arrays defining mars model (output from mars - see above).
c
c output:
c fmn,imn = corresponding arrays defining the sliced model
c             (dimensioned same as in mars - see above).
c
c workspace:
c sp(2*nk+2*p+max(nk,3*p)) : real.
c mm(2*p) : integer.
c
c
c
c call fmod (m,n,x,fm,im,f,sp):
c
c calculates mars model response estimates for sets of covariate vectors.
c
c input:
c m = model flag:
c   = 1 => piecewise-linear mars model.
c   = 2 => piecewise-cubic mars model. (Ref[2] Sec. 3.7)
c n = number of covariate vectors.
c x(n,p) = covariate vectors.
c fm,im = same as in mars (see above).
c
c output:
c f(n) = value of the mars model estimate for each covariate vector.
c
c workspace:
c sp(n,2) : real.
c
c
c
c call cvinfo (dfs,pse,nbf):
c 
c returns results of sample reuse procedure for estimating optimal smoothing
c parameter df (see above). can only be called if xvalid(ix) was called
c before mars with ix.ne.0 (see above). (Ref[2] Sec 3.6, Ref[3] Sec. 2.3)
c
c output:
c dfs = optimal smoothing parameter estimate.
c pse = estimate of corresponding predictive-squared-error.
c nbf = estimate of associated number of (nonconstant) basis functions.
c
c 
c
      subroutine mars (n,p,x,y,w,nk,mi,lx,fm,im,sp,dp,mm)                   1
      integer p,lx(*),im(*),mm(*)                                           2
      real x(*),y(*),w(*),fm(*),sp(*)                                       3
      double precision dp(*)                                                4
      im(3)=n                                                               5
      im(4)=p                                                               6
      im(5)=nk                                                              7
      im(6)=mi                                                              8
      im(7)=16                                                              9
      im(8)=im(7)+5*nk                                                     10
      im(9)=im(8)+2*nk*mi                                                  11
      im(10)=im(9)+3*(nk+2)                                                12
      im(2)=im(10)+nk*mi-1                                                 13
      im(11)=1                                                             14
      im(12)=2                                                             15
      im(13)=im(12)+5*nk                                                   16
      im(14)=im(13)+1                                                      17
      im(15)=im(14)+nk*(5*mi+1)                                            18
      call mars1(n,p,x,y,w,nk,mi,lx,fm(im(11)),fm(im(12)),fm(im(15)),im(   19
     1im(7)),  im(im(8)),im(im(9)),im(im(10)),fm(im(13)),fm(im(14)),sp,d   20
     1p,mm)                                                                21
      im(1)=im(15)+lcm(p,nk,fm(im(12)),fm(im(15)))-1                       22
      return                                                               23
      end                                                                  24
      subroutine plot (m,x,fm,im,ngc,ngs,icx,nc,crv,ns,srf,sp,mm)          25
      integer im(*),mm(*)                                                  26
      real x(*),fm(*),crv(*),srf(*),sp(*)                                  27
      if(m .ne. 1) go to 1                                                 28
      call plotl(im(3),im(4),x,im(5),im(im(7)),im(im(8)),im(im(9)),im(im   29
     1(10)),  fm(im(12)),fm(im(15)),ngc,ngs,icx,nc,crv,ns,srf,sp,mm)       30
      return                                                               31
    1 call plotc(im(3),im(4),x,im(5),im(im(7)),im(im(8)),im(im(9)),im(im   32
     1(10)),  fm(im(14)),fm(im(15)),ngc,ngs,icx,nc,crv,ns,srf,sp,mm)       33
      return                                                               34
      end                                                                  35
      subroutine catprt (m,fm,im,sp,mm)                                    36
      integer im(*),mm(*)                                                  37
      real fm(*),sp(*)                                                     38
      call ctprt1(m,im(5),im(im(7)),im(im(8)),fm(im(12)),fm(im(15)),fm(i   39
     1m(14)),sp,mm)                                                        40
      return                                                               41
      end                                                                  42
      subroutine slice (flg,xs,x,fm,im,fmn,imn,sp,mm)                      43
      integer im(*),imn(*),mm(*)                                           44
      real xs(*),x(*),fm(*),fmn(*),sp(*)                                   45
      do 1 i=1,15                                                          46
      imn(i)=im(i)                                                         47
    1 continue                                                             48
      i=im(15)                                                             49
      go to 3                                                              50
    2 i=i+1                                                                51
    3 if((i).gt.(im(1))) go to 4                                           52
      fmn(i)=fm(i)                                                         53
      go to 2                                                              54
    4 call slice1(flg,xs,im(3),im(4),x,im(5),fm(im(11)),fm(im(12)),fm(im   55
     1(15)),  im(im(7)),im(im(8)),im(im(9)),im(im(10)),fm(im(13)),fm(im(   56
     114)),  fmn(im(11)),fmn(im(12)),imn(im(7)),imn(im(8)),imn(im(9)),im   57
     1n(im(10)),  fmn(im(13)),fmn(im(14)),sp,mm)                           58
      return                                                               59
      end                                                                  60
      subroutine fmod (m,n,x,fm,im,f,sp)                                   61
      integer im(*)                                                        62
      real x(*),fm(*),f(*),sp(*)                                           63
      if(m .ne. 1) go to 1                                                 64
      call fmrs(n,x,im(5),fm(im(11)),fm(im(12)),fm(im(15)),f)              65
      return                                                               66
    1 call cmrs(n,x,fm(im(15)),im(im(7)),im(im(8)),im(im(9)),im(im(10)),   67
     1  fm(im(13)),fm(im(14)),f,sp)                                        68
      return                                                               69
      end                                                                  70
      subroutine print(it)                                                 71
      call printm(it)                                                      72
      call printg(it)                                                      73
      call printc(it)                                                      74
      call prtslc(it)                                                      75
      return                                                               76
      end                                                                  77
      subroutine setint(i,j,k)                                             78
      parameter(mlist=1000)                                                79
      integer m(2,mlist)                                                   80
      save m                                                               81
      data il /0/                                                          82
      if((i .ne. 0) .and. (j .ne. 0)) go to 1                              83
      il=0                                                                 84
      return                                                               85
    1 if(i.eq.j) return                                                    86
      m1=min0(i,j)                                                         87
      m2=max0(i,j)                                                         88
      if(k .ne. 0) go to 6                                                 89
      l=1                                                                  90
      go to 3                                                              91
    2 l=l+1                                                                92
    3 if((l).gt.(il)) go to 4                                              93
      if(m1.eq.m(1,l).and.m2.eq.m(2,l)) return                             94
      go to 2                                                              95
    4 il=il+1                                                              96
      if(il .le. mlist) go to 5                                            97
c     write(6,  '('' increase parameter mlist in subroutine setint to gr   98
c    1eater than'',           i5,/,'' and recompile.'')') il               99
      stop                                                                100
    5 m(1,il)=m1                                                          101
      m(2,il)=m2                                                          102
      return                                                              103
    6 ig=0                                                                104
      l=1                                                                 105
      go to 8                                                             106
    7 l=l+1                                                               107
    8 if((l).gt.(il)) go to 10                                            108
      if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 7                      109
      ig=1                                                                110
   10 if(ig.eq.0) return                                                  111
      il=il-1                                                             112
      ll=l                                                                113
      go to 12                                                            114
   11 ll=ll+1                                                             115
   12 if((ll).gt.(il)) go to 13                                           116
      m(1,ll)=m(1,ll+1)                                                   117
      m(2,ll)=m(2,ll+1)                                                   118
      go to 11                                                            119
   13 return                                                              120
      entry intlst(it)                                                    121
      if(it.le.0) return                                                  122
      if(il.eq.0) return                                                  123
      write(it,'(/,'' interactions prohibited between:'')')               124
c     do 14 l=1,il                                                        125
c     write(it,'(''    var('',i3,'')  and  var('',i3,'')'')') m(1,l),m(2  126
c    1,l)                                                                 127
c  14 continue                                                            128
      return                                                              129
      entry intalw(i,j,k)                                                 130
      k=1                                                                 131
      m1=min0(i,j)                                                        132
      m2=max0(i,j)                                                        133
      l=1                                                                 134
      go to 16                                                            135
   15 l=l+1                                                               136
   16 if((l).gt.(il)) go to 18                                            137
      if(m1 .ne. m(1,l) .or. m2 .ne. m(2,l)) go to 15                     138
      k=0                                                                 139
   18 return                                                              140
      end                                                                 141
      subroutine mars1 (n,p,x,y,w,nk,mi,lx,az,tb,cm,kp,kv,lp,lv,bz,tc,sp  142
     1,dp,mm)                                                             143
      integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),mm(n,*),lx(p)               144
      real x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),sp(*)                    145
      double precision dp(*)                                              146
      data ms,df,il,fv,it,ic,ix /0,3.0,0,0.0,6,0,0/                       147
c     if(it.gt.0) write(it,11)                                            148
c     if(it.gt.0) write(it,10) n,p,nk,ms,mi,df,il,fv,ic                   149
c     if(it.gt.0) write(it,12)                                            150
c     if(it.gt.0) write(it,'('' var: '',5('' '',20i3,/))') (i,i=1,p)      151
c     if(it.gt.0) write(it,'('' flag:'',5('' '',20i3,/))') (lx(i),i=1,p)  152
c     print *, ' '
c     do 321 i = 1, n
c        print *,'M1 ',x(i,1),' ',x(i,2),' ',x(i,3),' ',x(i,4),
c    1         ' ',y(i)
c321  continue
      call intlst(it)                                                     153
      call nstlst(it)                                                     154
      i1=max0(n*(nk+1),2*n)+1                                             155
      im=i1+n+max0(3*n+5*nk,2*p,4*n,2*n+5*nk+p)                           156
      is=im+p                                                             157
      i2=max0(n*nk,(nk+1)*(nk+1))+1                                       158
      call rspnpr(it,il,n,y,w,mm)                                         159
      do 2 j=1,p                                                          160
      do 1 i=1,n                                                          161
      mm(i,j)=i                                                           162
    1 continue                                                            163
      call psort(x(1,j),mm(1,j),1,n)                                      164
    2 continue                                                            165
      call ordpr(it,n,p,x,lx,mm)                                          166
      call atoscl (n,p,w,x,lx,mm,sp(im),sp(is),cm,x)                      167
      call catpr(it,n,p,x,cm,mm(1,p+1))                                   168
      call oknest(it,p,lx,cm)                                             169
      if(ix.ne.0) call cvmars (ix,n,p,x,y,w,nk,ms,df,fv,mi,lx,it,sp(im),  170
     1sp(is),tb,cm,sp,dp,dp(i2),mm, sp(is+p),sp(is+p+2*n))                171
      call marsgo  (n,p,x,y,w,nk,ms,df,fv,mi,lx,it,sp(im),sp(is),az,tb,c  172
     1m,sp,dp,dp(i2),mm)                                                  173
      if(il .le. 0) go to 6                                               174
      call logitl(n,x,y,w,nk,il,az,tb,cm,sp,dp)                           175
      if(it .le. 0) go to 6                                               176
      sw=0.0                                                              177
      wn=sw                                                               178
      do 3 i=1,n                                                          179
      sw=sw+w(i)                                                          180
      wn=wn+w(i)**2                                                       181
    3 continue                                                            182
      wn=sw**2/wn                                                         183
      ef=1.0                                                              184
      do 4 k=1,nk                                                         185
      if(tb(1,k).ne.0.0) ef=ef+tb(5,k)                                    186
    4 continue                                                            187
      ef=1.0/(1.0-ef/wn)**2                                               188
      s=0.0                                                               189
      t=s                                                                 190
      call fmrs(n,x,nk,az,tb,cm,sp)                                       191
      do 5 i=1,n                                                          192
      yh=1.0/(1.0+exp(-sp(i)))                                            193
      gcv=ef*(y(i)-yh)**2                                                 194
      s=s+w(i)*gcv                                                        195
      t=t+w(i)*yh*(1.0-yh)                                                196
    5 continue                                                            197
      s=s/sw                                                              198
      t=t/sw                                                              199
c     write(it,13) s,t                                                    200
    6 if(it .le. 0) go to 7                                               201
      if(il.eq.0) call anova (n,x,y,w,nk,it,tb,cm,lp,lv,sp,dp)            202
      if(il.gt.0) call anoval(n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sp,dp)      203
    7 call ccoll (nk,tb,cm,kp,kv,lp,lv,mm)                                204
      call cubic (n,p,x,y,w,nk,it,tb,cm,kp,kv,lp,lv,bz,tc,sp,sp(i1),sp(i  205
     11+2*p),mm,dp)                                                       206
      if(il .le. 0) go to 9                                               207
      call logitc(n,x,y,w,nk,il,cm,kp,kv,lp,lv,bz,tc,sp,sp(i1+4*n),dp)    208
      if(it .le. 0) go to 9                                               209
      call cmrs(n,x,cm,kp,kv,lp,lv,bz,tc,sp,sp(n+1))                      210
      s=0.0                                                               211
      t=s                                                                 212
      do 8 i=1,n                                                          213
      yh=1.0/(1.0+exp(-sp(i)))                                            214
      gcv=ef*(y(i)-yh)**2                                                 215
      s=s+w(i)*gcv                                                        216
      t=t+w(i)*yh*(1.0-yh)                                                217
    8 continue                                                            218
      s=s/sw                                                              219
      t=t/sw                                                              220
c     write(it,14) s,t                                                    221
    9 if(it.gt.0) call varimp (n,p,x,y,w,nk,il,it,az,tb,cm,sp,sp(p+1),dp  222
     1)                                                                   223
      call orgpl(sp(im),sp(is),nk,tb,cm)                                  224
      call orgpc(sp(im),sp(is),lp,lv,tc)                                  225
      call sclato(n,p,x,sp(im),sp(is),cm,x)                               226
      return                                                              227
      entry setms(mal)                                                    228
      ms=mal                                                              229
      return                                                              230
      entry setdf(val)                                                    231
      df=val                                                              232
      return                                                              233
      entry printm(mal)                                                   234
      it=mal                                                              235
      return                                                              236
      entry logit(mal)                                                    237
      il=mal                                                              238
      return                                                              239
      entry setfv(val)                                                    240
      fv=val                                                              241
      return                                                              242
      entry setic(mal)                                                    243
      ic=mal                                                              244
      z00001=stelg(ic)                                                    245
      return                                                              246
      entry xvalid(mal)                                                   247
      ix=mal                                                              248
      call xvmrgo(ix)                                                     249
      return                                                              250
   10 format(/' input parameters (see doc.):',/,  '    n     p    nk      251
     1ms    mi     df    il    fv     ic',/,  ' ',i5,i5,i6,i6,i6,f8.3,i5  252
     1,f7.3,i6)                                                           253
   11 format(//,' MARS modeling, version 3.6 (3/25/93)',/)                254
   12 format(/' predictor variable flags:')                               255
   13 format(/' piecewise-linear logistic gcv =',g12.4,'   ave var =',g1  256
     12.4)                                                                257
   14 format(/' piecewise-cubic logistic gcv =',g12.4,'   ave var =',g12  258
     1.4)                                                                 259
      end                                                                 260
      subroutine plotc (n,p,x,nk,kp,kv,lp,lv,tc,cm,ngc,ngs,icx,nc,crv,ns  261
     1,srf,sp,mm)                                                         262
      integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),mm(*)                       263
      real x(n,p),tb(5,nk),tc(*),cm(*),crv(ngc,2,*),srf(ngs,ngs,*),sp(*)  264
     1,zl(2),zu(2)                                                        265
      data big,it /1.e30,6/                                               266
c     if(it.gt.0) write(it,'(/'' mars graphics (piecewise-cubic):'',/)')  267
      jnt=2                                                               268
      go to 1                                                             269
      entry plotl (n,p,x,nk,kp,kv,lp,lv,tb,cm,ngc,ngs,icx,nc,crv,ns,srf,  270
     1sp,mm)                                                              271
c     if(it.gt.0) write(it,'(/'' mars graphics (piecewise-linear):'',/)'  272
c    1)                                                                   273
      jnt=1                                                               274
    1 ngsq=ngs**2                                                         275
      iz=2*ngsq                                                           276
      d=1.0/(ngs-1)                                                       277
      dc=1.0/(ngc-1)                                                      278
      ll=1                                                                279
      nc=0                                                                280
      ns=nc                                                               281
    2 if(kp(1,ll).lt.0) go to 36                                          282
      if(kp(3,ll) .gt. 0) go to 3                                         283
      ll=ll+1                                                             284
      go to 2                                                             285
    3 nf=kp(3,ll)                                                         286
      k4=kp(4,ll)-1                                                       287
      k1=kp(1,ll)                                                         288
      k2=kp(2,ll)                                                         289
      if(it .le. 0) go to 7                                               290
      if(k1 .ne. 0) go to 4                                               291
c     write(it,'('' pure ordinal contribution:'')')                       292
      go to 7                                                             293
    4 continue
c     write(it,'('' categorical - ordinal interaction:'')')               294
      do 6 i=1,k1                                                         295
      jj=kv(1,k2+i-1)                                                     296
      j=iabs(jj)                                                          297
      k=kv(2,k2+i-1)                                                      298
      ncx=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                             299
      do 5 l=1,ncx                                                        300
      mm(l)=cm(k+l)+.1                                                    301
      if(jj.lt.0) mm(l)=mod(mm(l)+1,2)                                    302
    5 continue                                                            303
c     write(it,'('' x('',i3,'') ='',70i1/80i1)') j,(mm(l),l=1,ncx)        304
    6 continue                                                            305
    7 do 35 k=1,nf                                                        306
      l=lp(1,k+k4)                                                        307
      if(l.gt.2) go to 35                                                 308
      ko=lp(2,k+k4)                                                       309
      if(l .ne. 1) go to 17                                               310
      j=0                                                                 311
      jv=lv(ko)                                                           312
      do 9 m=k,nf                                                         313
      l1=lp(1,m+k4)                                                       314
      if(l1.eq.1) go to 9                                                 315
      l2=lp(2,m+k4)-1                                                     316
      do 8 i=1,l1                                                         317
      if(jv.eq.lv(l2+i)) j=1                                              318
    8 continue                                                            319
      if(j.eq.1) go to 10                                                 320
    9 continue                                                            321
   10 if(j.eq.1) go to 35                                                 322
      nc=nc+1                                                             323
      zl(1)=big                                                           324
      zu(1)=-big                                                          325
      do 11 i=1,n                                                         326
      r=x(i,jv)                                                           327
      zl(1)=amin1(zl(1),r)                                                328
      zu(1)=amax1(zu(1),r)                                                329
   11 continue                                                            330
      dl=(zu(1)-zl(1))*dc                                                 331
      do 12 i=1,ngc                                                       332
      crv(i,1,nc)=zl(1)+dl*(i-1)                                          333
   12 continue                                                            334
      if(jnt .ne. 1) go to 13                                             335
      call fun(l,jv,ngc,crv(1,1,nc),nk,tb,cm,k1,kv(1,k2),crv(1,2,nc),mm)  336
      go to 14                                                            337
   13 call cfun (l,jv,ngc,crv(1,1,nc),nf,lp(1,k4+1),lv,tc(kp(5,ll)),  cr  338
     1v(1,2,nc),sp,mm)                                                    339
   14 dl=big                                                              340
      do 15 i=1,ngc                                                       341
      dl=amin1(dl,crv(i,2,nc))                                            342
   15 continue                                                            343
      fx=0.0                                                              344
      do 16 i=1,ngc                                                       345
      crv(i,2,nc)=crv(i,2,nc)-dl                                          346
      fx=amax1(fx,crv(i,2,nc))                                            347
   16 continue                                                            348
c     if(it.gt.0) write(it,39) nc,jv,fx                                   349
      go to 35                                                            350
   17 j=0                                                                 351
      mm(1)=lv(ko)                                                        352
      mm(2)=lv(ko+1)                                                      353
      do 19 m=k,nf                                                        354
      l1=lp(1,m+k4)                                                       355
      if(l1.le.2) go to 19                                                356
      l2=lp(2,m+k4)-1                                                     357
      do 18 i=1,l1                                                        358
      if(mm(1).eq.lv(l2+i).or.mm(2).eq.lv(l2+i)) j=1                      359
   18 continue                                                            360
      if(j.eq.1) go to 20                                                 361
   19 continue                                                            362
   20 if(j.eq.1) go to 35                                                 363
      ns=ns+1                                                             364
      zl(1)=big                                                           365
      zl(2)=zl(1)                                                         366
      zu(1)=-big                                                          367
      zu(2)=zu(1)                                                         368
      do 22 j=1,2                                                         369
      do 21 i=1,n                                                         370
      r=x(i,mm(j))                                                        371
      zl(j)=amin1(zl(j),r)                                                372
      zu(j)=amax1(zu(j),r)                                                373
   21 continue                                                            374
   22 continue                                                            375
      do 23 j=1,2                                                         376
      dl=(zu(j)-zl(j))/(ngs-3)                                            377
      zu(j)=zu(j)+dl                                                      378
      zl(j)=zl(j)-dl                                                      379
   23 continue                                                            380
      ne=0                                                                381
      d1=d*(zu(1)-zl(1))                                                  382
      d2=d*(zu(2)-zl(2))                                                  383
      do 25 j=1,ngs                                                       384
      do 24 i=1,ngs                                                       385
      ne=ne+1                                                             386
      sp(iz+ne)=zl(1)+d1*(i-1)                                            387
      sp(iz+ngsq+ne)=zl(2)+d2*(j-1)                                       388
   24 continue                                                            389
   25 continue                                                            390
      dl=big                                                              391
      if(jnt .ne. 1) go to 26                                             392
      call pair(mm,ngsq,sp(iz+1),nk,tb,cm,k1,kv(1,k2),  srf(1,1,ns),sp,m  393
     1m(3))                                                               394
      go to 27                                                            395
   26 call cpair(mm,ngsq,sp(iz+1),nf,lp(1,k4+1),lv,  tc(kp(5,ll)),srf(1,  396
     11,ns),sp)                                                           397
   27 if(icx .le. 0) go to 29                                             398
      call cvxhul(n,x(1,mm(1)),x(1,mm(2)),big,nh,sp)                      399
      if(it .le. 0 .or. 3*nh .lt. iz) go to 28                            400
      nxs=sqrt(float(3*nh)*0.5)+1.1                                       401
c     write(it,38) nxs                                                    402
   28 call hulset(ngsq,sp(iz+1),big,nh,sp,srf(1,1,ns))                    403
   29 do 31 j=1,ngs                                                       404
      do 30 i=1,ngs                                                       405
      if(i.eq.1.or.j.eq.1.or.i.eq.ngs.or.j.eq.ngs.or.srf(i,j,ns).ge.big)  406
     1 go to 30                                                           407
      dl=amin1(dl,srf(i,j,ns))                                            408
   30 continue                                                            409
   31 continue                                                            410
      fx=0.0                                                              411
      do 34 j=1,ngs                                                       412
      do 33 i=1,ngs                                                       413
      if((i .ne. 1) .and. ((j .ne. 1) .and. ((i .ne. ngs) .and. ((j .ne.  414
     1 ngs) .and. (srf(i,j,ns) .lt. big))))) go to 32                     415
      srf(i,j,ns)=0.0                                                     416
      go to 33                                                            417
   32 srf(i,j,ns)=srf(i,j,ns)-dl                                          418
      fx=amax1(fx,srf(i,j,ns))                                            419
   33 continue                                                            420
   34 continue                                                            421
c     if(it.gt.0) write(it,40) ns,mm(1),mm(2),fx                          422
   35 continue                                                            423
      ll=ll+1                                                             424
      go to 2                                                             425
   36 continue
c     if(it.gt.0) write(it,37) nc,ns                                      426
      return                                                              427
      entry printg(nal)                                                   428
      it=nal                                                              429
      return                                                              430
   37 format(/,' ',i3,' curves and',i3,' surfaces.'/)                     431
   38 format(' plot: convex hull too large. increase ngs to',i6)          432
   39 format('   crv',i3,':  x(',i2,').  max =',g12.4)                    433
   40 format('   srf',i3,':  x(',i2,'), x(',i2,').  max =',g12.4)         434
      end                                                                 435
      subroutine ctprt1 (m,nk,kp,kv,tb,cm,tc,sc,js)                       436
      integer kp(5,*),kv(2,*),js(*)                                       437
      real cm(*),tc(*),sc(*)                                              438
      data big,it /9.9e30,6/                                              439
      if(it.le.0) return                                                  440
      nc=ncat(kp)                                                         441
      if(nc.eq.0) return                                                  442
      write(it,'(/,'' there are'',i3,'' purely categorical basis functio  443
     1ns.'')') nc                                                         444
c     write(it,'('' purely additive and bivariate contributions follow''  445
c    1)')                                                                 446
      if(m .ne. 1) go to 1                                                447
c     write(it,'('' (piecewise-linear fit):'')')                          448
      go to 2                                                             449
    1 continue
c     write(it,'('' (piecewise-cubic fit):'')')                           450
    2 call catv(1,kp,kv,nv,js)                                            451
      do 8 jj=1,nv                                                        452
      j=js(jj)                                                            453
      xm=big                                                              454
      xx=-big                                                             455
      nl=int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                              456
      do 3 i=1,nl                                                         457
      sc(i)=cvlv(m,1,j,i,nk,kp,kv,tb,cm,tc)                               458
      xm=amin1(xm,sc(i))                                                  459
      xx=amax1(xx,sc(i))                                                  460
    3 continue                                                            461
      px=99.0                                                             462
      if(nl.gt.26) px=9.0                                                 463
      rx=xx-xm                                                            464
      rxp=rx/px                                                           465
c     write(it,'(/,'' f( x('',i3,'') ) : scale ='',g12.4)') j,rxp         466
      if(rxp.le.0.0) go to 8                                              467
      do 4 i=1,nl                                                         468
      js(i+nv)=(sc(i)-xm)/rxp+.5                                          469
    4 continue                                                            470
      if(nl .gt. 26) go to 5                                              471
c     write(it,28) (i,i=1,nl)                                             472
c     write(it,28) (js(i+nv),i=1,nl)                                      473
      go to 8                                                             474
    5 if(nl .gt. 38) go to 6                                              475
c     write(it,29) (i,i=1,nl)                                             476
c     write(it,29) (js(i+nv),i=1,nl)                                      477
      go to 8                                                             478
    6 if(nl .gt. 78) go to 7                                              479
c     write(it,30) (mod(i,10),i=1,nl)                                     480
c     write(it,30) (js(i+nv),i=1,nl)                                      481
      go to 8                                                             482
    7 continue
c     write(it,37) 78                                                     483
    8 continue                                                            484
      call catv(2,kp,kv,nv,js)                                            485
      do 27 jj=1,nv                                                       486
      j1=js(2*jj-1)                                                       487
      j2=js(2*jj)                                                         488
      xm=big                                                              489
      xx=-big                                                             490
      n1=int(cm(2*j1+1)+.1)-int(cm(2*j1)+.1)+1                            491
      n2=int(cm(2*j2+1)+.1)-int(cm(2*j2)+.1)+1                            492
      k=0                                                                 493
      do 10 i=1,n1                                                        494
      s1=cvlv(m,1,j1,i,nk,kp,kv,tb,cm,tc)                                 495
      js(2*nv+1)=i                                                        496
      do 9 j=1,n2                                                         497
      js(2*nv+2)=j                                                        498
      k=k+1                                                               499
      sc(k)=s1+cvlv(m,2,js(2*jj-1),js(2*nv+1),nk,kp,kv,tb,cm,tc)          500
    9 continue                                                            501
   10 continue                                                            502
      do 12 j=1,n2                                                        503
      s1=cvlv(m,1,j2,j,nk,kp,kv,tb,cm,tc)                                 504
      do 11 i=1,n1                                                        505
      k=j+n2*(i-1)                                                        506
      sc(k)=sc(k)+s1                                                      507
   11 continue                                                            508
   12 continue                                                            509
      k=0                                                                 510
      do 14 i=1,n1                                                        511
      do 13 j=1,n2                                                        512
      k=k+1                                                               513
      x=sc(k)                                                             514
      xx=amax1(xx,x)                                                      515
      xm=amin1(xm,x)                                                      516
   13 continue                                                            517
   14 continue                                                            518
      na=min0(n1,n2)                                                      519
      nb=max0(n1,n2)                                                      520
      if(na .ne. n1) go to 15                                             521
      ja=j1                                                               522
      jb=j2                                                               523
      go to 16                                                            524
   15 ja=j2                                                               525
      jb=j1                                                               526
   16 px=99.0                                                             527
      if(na.gt.25) px=9.0                                                 528
      rx=xx-xm                                                            529
      rxp=rx/px                                                           530
c     write(it,'(/,'' f( x('',i3,''), x('',i3,'') ) : scale ='',g12.4)')  531
c    1ja,jb,rxp                                                           532
      if(rxp.le.0.0) go to 27                                             533
      if(na .le. 75) go to 17                                             534
c     write(it,37) 75                                                     535
      go to 27                                                            536
   17 if(na .gt. 25) go to 18                                             537
c     write(it,34) (i,i=1,na)                                             538
      go to 20                                                            539
   18 if(na .gt. 37) go to 19                                             540
c     write(it,35) (i,i=1,na)                                             541
      go to 20                                                            542
   19 continue
c     write(it,36) (mod(i,10),i=1,na)                                     543
   20 do 26 j=1,nb                                                        544
      do 23 i=1,na                                                        545
      if(na .ne. n1) go to 21                                             546
      k=j+n2*(i-1)                                                        547
      go to 22                                                            548
   21 k=i+n2*(j-1)                                                        549
   22 js(i+2*nv)=(sc(k)-xm)/rxp+.5                                        550
   23 continue                                                            551
      if(na .gt. 25) go to 24                                             552
c     write(it,31) j,(js(i+2*nv),i=1,na)                                  553
      go to 26                                                            554
   24 if(na .gt. 37) go to 25                                             555
c     write(it,32) j,(js(i+2*nv),i=1,na)                                  556
      go to 26                                                            557
   25 continue
c     write(it,33) j,(js(i+2*nv),i=1,na)                                  558
   26 continue                                                            559
   27 continue                                                            560
      return                                                              561
      entry printc(nal)                                                   562
      it=nal                                                              563
      return                                                              564
   28 format(' ',26i3)                                                    565
   29 format(' ',38i2)                                                    566
   30 format(' ',78i1)                                                    567
   31 format(' ',i3,' ',25i3)                                             568
   32 format(' ',i3,' ',37i2)                                             569
   33 format(' ',i3,' ',75i1)                                             570
   34 format('     ',25i3)                                                571
   35 format('     ',37i2)                                                572
   36 format('     ',75i1)                                                573
   37 format(' function not printed (more than',i3,' categorical values)  574
     1.')                                                                 575
      end                                                                 576
      subroutine slice1 (flg,xs,n,p,x,nk,az,tb,cm,kp,kv,lp,lv,bz,tc,azn,  577
     1tbn,kpn,kvn,  lpn,lvn,bzn,tcn,sp,mm)                                578
      integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),kpn(5,*),kvn(2,*),lpn(3,*)  579
     1,lvn(*),mm(*)                                                       580
      real xs(p),x(n,p),tb(5,nk),cm(*),tc(*),tbn(5,nk),tcn(*),sp(*)       581
      data it,big /6,9.9e30/                                              582
      ni=0                                                                583
      do 1 m=1,nk                                                         584
      if(tb(1,m).ne.0.0) ni=ni+1                                          585
    1 continue                                                            586
      if(ni .ne. 0) go to 3                                               587
      kpn(1,1)=-1                                                         588
      lpn(1,1)=0                                                          589
      do 2 m=1,nk                                                         590
      tbn(1,m)=0.0                                                        591
    2 continue                                                            592
      azn=0.0                                                             593
      bzn=azn                                                             594
c     if(it.gt.0) write(it,'('' slice: original mars model = constant.''  595
c    1)')                                                                 596
      return                                                              597
    3 if(it .le. 0) go to 5                                               598
c     write(it,'(/,'' sliced mars model: flag ='',g12.4)') flg            599
c     write(it,'(/,'' slice:'')')                                         600
      do 4 j=1,p                                                          601
      if(xs(j).eq.flg) go to 4                                            602
c     write(it,'('' x('',i3,'') ='',g12.4)') j,xs(j)                      603
    4 continue                                                            604
    5 i1=2*nk+1                                                           605
      i2=i1+2*p                                                           606
      i3=max0(i2+p,i1+nk)                                                 607
      do 7 j=1,p                                                          608
      xl=big                                                              609
      xr=-xl                                                              610
      do 6 i=1,n                                                          611
      xl=amin1(xl,x(i,j))                                                 612
      xr=amax1(xr,x(i,j))                                                 613
    6 continue                                                            614
      sp(j+i3-1)=xr-xl                                                    615
      sp(j+i3-1+p)=xl                                                     616
    7 continue                                                            617
      call reducq(flg,xs,nk,tb,cm,tc,kp,kv,lp,lv,sp(i3),sp,sp(i1),sp(i2)  618
     1)                                                                   619
      call reducl(flg,xs,nk,az,tb,cm,bz,sp,sp(i3),azn,tbn,bzn,sp(i1))     620
      ni=0                                                                621
      do 8 m=1,nk                                                         622
      if(tbn(1,m).ne.0.0) ni=ni+1                                         623
    8 continue                                                            624
      if(ni .ne. 0) go to 10                                              625
      kpn(1,1)=-1                                                         626
      lpn(1,1)=0                                                          627
      do 9 m=1,nk                                                         628
      tbn(1,m)=0.0                                                        629
    9 continue                                                            630
      azn=0.0                                                             631
      bzn=azn                                                             632
c     if(it.gt.0) write(it,'('' sliced mars model = constant.'')')        633
      return                                                              634
   10 if(it.gt.0) call slova(nk,it,tbn,ni,lpn,lvn)                        635
      call ccoll(nk,tbn,cm,kpn,kvn,lpn,lvn,mm)                            636
      call qslice(p,nk,tbn,cm,sp,kpn,kvn,lpn,lvn,tcn,sp(i3),sp(i1),mm)    637
      return                                                              638
      entry prtslc(ig)                                                    639
      it=ig                                                               640
      return                                                              641
      end                                                                 642
      subroutine cmrs (n,x,cm,kp,kv,lp,lv,bz,tc,y,sc)                     643
      integer kp(5,*),kv(2,*),lp(3,*),lv(*)                               644
      real x(n,*),tc(*),cm(*),y(n),sc(n,2)                                645
      data ifg /0/                                                        646
      do 1 i=1,n                                                          647
      y(i)=bz                                                             648
    1 continue                                                            649
      ll=1                                                                650
      la=ll                                                               651
      l1=la                                                               652
    2 if(kp(1,ll).lt.0) go to 19                                          653
      do 3 i=1,n                                                          654
      sc(i,1)=1.0                                                         655
    3 continue                                                            656
      if(kp(1,ll) .le. 0) go to 11                                        657
      jl=kp(1,ll)                                                         658
      do 10 il=1,jl                                                       659
      k=kp(2,ll)+il-1                                                     660
      jj=kv(1,k)                                                          661
      j=iabs(jj)                                                          662
      kk=kv(2,k)                                                          663
      do 9 i=1,n                                                          664
      if(sc(i,1).eq.0.0) go to 9                                          665
      if(ifg .ne. 0) go to 4                                              666
      ic=icat(x(i,j),j,cm)                                                667
      go to 5                                                             668
    4 ic=x(i,j)+.1                                                        669
    5 if(ic .ne. 0) go to 6                                               670
      sc(i,1)=0.0                                                         671
      go to 7                                                             672
    6 sc(i,1)=cm(ic+kk)                                                   673
    7 if(jj .ge. 0) go to 9                                               674
      if(sc(i,1) .ne. 0.0) go to 8                                        675
      sc(i,1)=1.0                                                         676
      go to 9                                                             677
    8 sc(i,1)=0.0                                                         678
    9 continue                                                            679
   10 continue                                                            680
      go to 12                                                            681
   11 if(kp(3,ll) .gt. 0) go to 12                                        682
      ll=ll+1                                                             683
      go to 2                                                             684
   12 if(kp(3,ll) .ge. 0) go to 14                                        685
      k=-kp(3,ll)                                                         686
      ll=ll+1                                                             687
      do 13 i=1,n                                                         688
      if(sc(i,1).eq.0.0) go to 13                                         689
      y(i)=y(i)+tc(k)                                                     690
   13 continue                                                            691
      go to 2                                                             692
   14 kp3=kp(3,ll)                                                        693
      do 18 m=1,kp3                                                       694
      l=lp(1,l1)                                                          695
      nt=lp(3,l1)                                                         696
      lb=la+5*l*nt-1                                                      697
      do 17 j=1,nt                                                        698
      do 15 i=1,n                                                         699
      sc(i,2)=sc(i,1)                                                     700
   15 continue                                                            701
      call que(j,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,2))                    702
      do 16 i=1,n                                                         703
      y(i)=y(i)+tc(lb+j)*sc(i,2)                                          704
   16 continue                                                            705
   17 continue                                                            706
      la=lb+nt+1                                                          707
      l1=l1+1                                                             708
   18 continue                                                            709
      ll=ll+1                                                             710
      go to 2                                                             711
   19 return                                                              712
      entry stcmrs(i1)                                                    713
      ifg=i1                                                              714
      return                                                              715
      end                                                                 716

      subroutine fmrs (n,x,nk,az,tb,cm,y)                                 717
      real x(n,*),tb(5,nk),cm(*),y(n)                                     718
      double precision s                                                  719
      data ifg /0/                                                        720
      do 13 i=1,n                                                         721
      s=az                                                                722
      do 12 m=1,nk                                                        723
      if(tb(1,m).eq.0.0) go to 12                                         724
      phi=1.0                                                             725
      ip=m                                                                726
    1 if(ip.le.0) go to 11                                                727
      t=tb(2,ip)                                                          728
      j=abs(t)+.1                                                         729
      if(cm(2*j) .le. 0.0) go to 8                                        730
      if(ifg .ne. 0) go to 2                                              731
      k=icat(x(i,j),j,cm)                                                 732
      go to 3                                                             733
    2 k=x(i,j)+.1                                                         734
    3 if(k .ne. 0) go to 4                                                735
      u=0.0                                                               736
      go to 5                                                             737
    4 u=cm(k+int(tb(3,ip)+.1))                                            738
    5 if(t .ge. 0.0) go to 9                                              739
      if(u .ne. 0.0) go to 6                                              740
      u=1.0                                                               741
      go to 9                                                             742
    6 u=0.0                                                               743
      go to 9                                                             744
    8 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,ip)))                          745
    9 if(u .ne. 0.0) go to 10                                             746
      phi=0.0                                                             747
      go to 11                                                            748
   10 phi=phi*u                                                           749
      ip=tb(4,ip)+.1                                                      750
      go to 1                                                             751
   11 s=s+tb(1,m)*phi                                                     752
   12 continue                                                            753
      y(i)=s                                                              754
   13 continue                                                            755
      return                                                              756
      entry stfmrs(i1)                                                    757
      ifg=i1                                                              758
      return                                                              759
      end                                                                 760

      subroutine marsgo (n,p,x,y,w,nk,ms,df,fv,mi,lx,it,xm,xs,az,tb,cm,s  761
     1c,db,d,mm)                                                          762
      integer p,mm(n,p),lx(p)                                             763
      logical elg,newbf                                                   764
      real x(n,p),y(n),w(n),xm(p),xs(p),tb(5,nk),cm(*),sc(n,*),vcst(3),t  765
     1x(5)                                                                766
      double precision db(n,*),d(nk,*)                                    767
      double precision yb,yv,sw,s,t,u,v,we,sy,a,b,xb,xx,xd,ssq,alr        768
      double precision dx,wn,se,tt,txt,xt,st,su,yc,eps,rsq,dy,dv,asq0     769
      character*28 hol                                                    770
      data ix,alr,eps,big,fln,nmin,alf,vcst  /0,1.d-7,1.d-4,9.9e30,-1.0,  771
     15,.05,1.0,.666667,.333333/                                          772
c     if(it.gt.0) write(it,97)                                            773
      mk=nk                                                               774
      df1=0.0                                                             775
      nep=0                                                               776
      do 1 j=1,p                                                          777
      if(lx(j).eq.0) go to 1                                              778
      if(x(mm(1,j),j).ge.x(mm(n,j),j)) go to 1                            779
      nep=nep+1                                                           780
      cst=vcst(iabs(lx(j)))                                               781
      if(mi.eq.1) cst=amin1(cst,vcst(2))                                  782
      df1=df1+cst                                                         783
    1 continue                                                            784
      if(nep .ne. 0) go to 2                                              785
c     if(it.gt.0) write(it,'('' no predictor variables.'')')              786
      stop                                                                787
    2 if(nep.eq.1) df1=vcst(3)                                            788
      cfac=df1/nep                                                        789
      df1=df*cfac                                                         790
      mkp1=mk+1                                                           791
      mkp2=mk+2                                                           792
      sw=0.d0                                                             793
      wn=sw                                                               794
      yb=wn                                                               795
      s=yb                                                                796
      do 3 i=1,n                                                          797
      sw=sw+w(i)                                                          798
      wn=wn+w(i)**2                                                       799
      yb=yb+w(i)*y(i)                                                     800
    3 continue                                                            801
      yb=yb/sw                                                            802
      wn=sw**2/wn                                                         803
      do 4 i=1,n                                                          804
      s=s+w(i)*(y(i)-yb)**2                                               805
    4 continue                                                            806
      yv=s/sw                                                             807
      tcst=1.0                                                            808
      tcmx=wn-df1*vcst(1)-2.0                                             809
      if(cm(1) .le. 0.0) go to 7                                          810
      i=2                                                                 811
      go to 6                                                             812
    5 i=i+(2)                                                             813
    6 if((2)*((i)-(2*p)).gt.0) go to 7                                    814
      if(cm(i).gt.0.0) kcp0=cm(i+1)+.1                                    815
      go to 5                                                             816
    7 m=0                                                                 817
      mtot=m                                                              818
      txm=yv/(1.d0-1.d0/wn)**2                                            819
      rsq=yv*sw                                                           820
      kr=0                                                                821
      nopt=0                                                              822
c     if(it.gt.0) write(it,98) m,txm,0.0,1.0                              823
      if(fln.lt.0.0) fln=1.0+4.0/wn                                       824
      call addpar(0)                                                      825
    8 if(m.ge.mk.or.tcst.ge.tcmx) go to 69                                826
      nopt=nopt+1                                                         827
      call itrpar(nopt)                                                   828
      mm1=m                                                               829
      m=m+1                                                               830
      txi=big                                                             831
      kcp=kcp0                                                            832
      asq0=rsq/sw                                                         833
    9 call nxtpar(l,jq)                                                   834
      if(l.lt.0) go to 53                                                 835
      txl=big                                                             836
      if(nnord(l,tb) .lt. mi) go to 10                                    837
      call updpar(0,-1.d0)                                                838
      go to 9                                                             839
   10 call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))                           840
      lbf=0                                                               841
      if(nnt .gt. nmin) go to 11                                          842
      call updpar(0,-1.d0)                                                843
      go to 9                                                             844
   11 nep=0                                                               845
      do 12 jp=1,p                                                        846
      if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 12                       847
      if(jf(l,jp,tb).ne.0) go to 12                                       848
      call isfac(l,jp,mm1,tb,cm,ja)                                       849
      if(ja.lt.0) go to 12                                                850
      if(.not.elg(jp,l,lx,tb,cm)) go to 12                                851
      nep=nep+1                                                           852
   12 continue                                                            853
      if(nep .ne. 0) go to 13                                             854
      call updpar(0,-1.d0)                                                855
      go to 9                                                             856
   13 call mnspan(ms,alf,nep,nnt,mn,me,mel)                               857
      if(nnt .gt. max0(me,mel)) go to 14                                  858
      call updpar(0,-1.d0)                                                859
      go to 9                                                             860
   14 if(jq .ne. 0) go to 15                                              861
      jd1=1                                                               862
      jd2=p                                                               863
      go to 16                                                            864
   15 jd1=jq                                                              865
      jd2=jd1                                                             866
   16 do 52 jp=jd1,jd2                                                    867
      if(x(mm(1,jp),jp).ge.x(mm(n,jp),jp)) go to 52                       868
      if(jf(l,jp,tb).ne.0) go to 52                                       869
      call isfac(l,jp,mm1,tb,cm,ja)                                       870
      if(ja.lt.0) go to 52                                                871
      if(.not.elg(jp,l,lx,tb,cm)) go to 52                                872
      if(ja .ne. 0) go to 18                                              873
      if(lbf .eq. 0) go to 19                                             874
      call blf0(l,0,n,x,w,cm,sc,nnt,sc(1,mkp1))                           875
      lbf=0                                                               876
      call mnspan(ms,alf,nep,nnt,mn,me,mel)                               877
      go to 19                                                            878
   18 call blf0(l,ja,n,x,w,cm,sc,nnt,sc(1,mkp1))                          879
      lbf=1                                                               880
      if(nnt.le.nmin) go to 52                                            881
      call mnspan(ms,alf,nep,nnt,mn,me,mel)                               882
      if(nnt.le.max0(me,mel)) go to 52                                    883
   19 fvr=1.0                                                             884
      if(jft(mm1,jp,tb).eq.0) fvr=1.0+fv                                  885
      ict=0                                                               886
      if(lx(jp) .ge. 0) go to 20                                          887
      ict=1                                                               888
      nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1                            889
      call csp(jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,nnt,  sw,me,mkp2,nop  890
     1,sc(1,mkp1),db,d(1,3),mm(1,p+1))                                    891
      if(nop.eq.0) go to 52                                               892
      go to 45                                                            893
   20 tb(2,m)=jp                                                          894
      tb(3,m)=x(mm(1,jp),jp)                                              895
      tb(4,m)=l                                                           896
      k1=kr                                                               897
      ssq=rsq                                                             898
      call update(1,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))   899
      if(kr .le. k1) go to 21                                             900
      rsq=rsq-d(kr,1)**2                                                  901
      tb(1,m)=rsq/sw                                                      902
      go to 22                                                            903
   21 tb(1,m)=big                                                         904
   22 if((lx(jp) .ne. 3) .and. ((m .lt. mk) .and. (nnt .gt. me+mel))) go  905
     1 to 26                                                              906
      tb(1,m)=rsq/sw                                                      907
      newbf=newb(m,tb).eq.0                                               908
      if(fvr*tb(1,m) .gt. txl .or. .not.(newbf)) go to 23                 909
      txl=fvr*tb(1,m)                                                     910
      tx1=tb(1,m)                                                         911
      jq=jp                                                               912
   23 if(fvr*tb(1,m) .gt. txi .or. .not.(newbf)) go to 25                 913
      txi=fvr*tb(1,m)                                                     914
      tx(1)=tb(1,m)                                                       915
      do 24 i=2,4                                                         916
      tx(i)=tb(i,m)                                                       917
   24 continue                                                            918
      jas=ja                                                              919
   25 kr=k1                                                               920
      rsq=ssq                                                             921
      go to 52                                                            922
   26 mm1=m                                                               923
      m=m+1                                                               924
      tb(1,m)=big                                                         925
      xa=0.0                                                              926
      j=n                                                                 927
      nnl=nnt                                                             928
      nst=0                                                               929
      nnr=-1                                                              930
   27 j0=j                                                                931
   28 mj=mm(j,jp)                                                         932
      h=sc(mj,mkp1)                                                       933
      if(w(mj) .le. 0.0 .or. h .le. 0.0) go to 29                         934
      nst=nst+1                                                           935
      nnl=nnl-1                                                           936
      nnr=nnr+1                                                           937
   29 if(x(mm(j-1,jp),jp).lt.x(mm(j,jp),jp) .and.nst.ge.mn.and.nnl.ge.me  938
     1l.and.nnr.ge.me) go to 30                                           939
      j=j-1                                                               940
      if(j.le.1) go to 30                                                 941
      go to 28                                                            942
   30 if(j.le.1) go to 45                                                 943
      nst=0                                                               944
      xb=xa                                                               945
      xa=x(mm(j,jp),jp)                                                   946
      if(j0 .ne. n) go to 34                                              947
      v=0.d0                                                              948
      u=v                                                                 949
      t=u                                                                 950
      we=t                                                                951
      se=we                                                               952
      sy=se                                                               953
      dy=sy                                                               954
      i=1                                                                 955
      go to 32                                                            956
   31 i=i+1                                                               957
   32 if((i).gt.(kr)) go to 33                                            958
      d(i,2)=0.d0                                                         959
      d(i,3)=d(i,2)                                                       960
      go to 31                                                            961
   33 txt=x(mm(1,jp),jp)+x(mm(n,jp),jp)                                   962
      xt=0.5*txt                                                          963
      go to 37                                                            964
   34 dx=xb-xa                                                            965
      dy=dy+dx*sy                                                         966
      we=we+dx*se                                                         967
      v=v+dx*(2.d0*u-(xb+xa-txt)*t)                                       968
      i=1                                                                 969
      go to 36                                                            970
   35 i=i+1                                                               971
   36 if((i).gt.(kr)) go to 37                                            972
      d(i,2)=d(i,2)+dx*d(i,3)                                             973
      go to 35                                                            974
   37 do 40 k=j,j0                                                        975
      mj=mm(k,jp)                                                         976
      h=sc(mj,mkp1)                                                       977
      if(w(mj).le.0.0.or.h.le.0.0) go to 40                               978
      xx=x(mj,jp)                                                         979
      xd=xx-xa                                                            980
      su=w(mj)*h                                                          981
      st=su*xd                                                            982
      yc=y(mj)-yb                                                         983
      dy=dy+st*yc                                                         984
      sy=sy+su*yc                                                         985
      we=we+st                                                            986
      se=se+su                                                            987
      sj=w(mj)*h**2                                                       988
      v=v+sj*xd**2                                                        989
      t=t+sj                                                              990
      u=u+sj*(xx-xt)                                                      991
      i=1                                                                 992
      go to 39                                                            993
   38 i=i+1                                                               994
   39 if((i).gt.(kr)) go to 40                                            995
      tt=db(mj,i)                                                         996
      d(i,2)=d(i,2)+st*tt                                                 997
      d(i,3)=d(i,3)+su*tt                                                 998
      go to 38                                                            999
   40 continue                                                           1000
      dv=v-we**2/sw                                                      1001
      if(dv .le. 0.d0) go to 44                                          1002
      a=0.d0                                                             1003
      b=a                                                                1004
      i=1                                                                1005
      go to 42                                                           1006
   41 i=i+1                                                              1007
   42 if((i).gt.(kr)) go to 43                                           1008
      s=d(i,2)                                                           1009
      a=a+s*d(i,1)                                                       1010
      b=b+s**2                                                           1011
      go to 41                                                           1012
   43 b=dv-b                                                             1013
      if(b .le. eps*dv) go to 44                                         1014
      b=-(dy-a)**2/b                                                     1015
      if(b .ge. tb(1,m)) go to 44                                        1016
      tb(1,m)=b                                                          1017
      tb(3,m)=xa                                                         1018
   44 j=j-1                                                              1019
      if(j.le.1) go to 45                                                1020
      go to 27                                                           1021
   45 tb(2,m)=jp                                                         1022
      tb(4,m)=l                                                          1023
      tb(1,m)=(rsq+tb(1,m))/sw                                           1024
      if(ict .ne. 0 .or. tb(1,mm1) .gt. fln*tb(1,m)) go to 46            1025
      mp=mm1                                                             1026
      go to 47                                                           1027
   46 mp=m                                                               1028
   47 newbf=newb(mp,tb).eq.0                                             1029
      if(fvr*tb(1,mp) .ge. txl .or. .not.(newbf)) go to 48               1030
      txl=fvr*tb(1,mp)                                                   1031
      tx1=tb(1,mp)                                                       1032
      jq=jp                                                              1033
   48 if(fvr*tb(1,mp) .ge. txi .or. .not.(newbf)) go to 51               1034
      txi=fvr*tb(1,mp)                                                   1035
      tx(1)=tb(1,mp)                                                     1036
      do 49 i=2,4                                                        1037
      tx(i)=tb(i,mp)                                                     1038
   49 continue                                                           1039
      jas=ja                                                             1040
      if(ict .eq. 0) go to 51                                            1041
      do 50 i=1,nc                                                       1042
      cm(kcp0+i)=cm(kcp+i)                                               1043
   50 continue                                                           1044
      kcp=kcp0+nc                                                        1045
      tx(3)=kcp0                                                         1046
   51 if(ict .ne. 0) go to 52                                            1047
      m=mm1                                                              1048
      mm1=m-1                                                            1049
      kr=k1                                                              1050
      rsq=ssq                                                            1051
   52 continue                                                           1052
      call updpar(jq,asq0-tx1)                                           1053
      go to 9                                                            1054
   53 jp=tx(2)+.1                                                        1055
      call selpar(int(tx(4)+.1))                                         1056
      if(cm(2*jp) .le. 0.) go to 54                                      1057
      nc=int(cm(2*jp+1)+.1)-int(cm(2*jp)+.1)+1                           1058
      kcp0=kcp0+nc                                                       1059
   54 if(jas .le. 0) go to 60                                            1060
      call getnst(jas,cm,jn,kcp,cm(kcp0+1))                              1061
      tb(2,m)=jn                                                         1062
      tb(3,m)=kcp0                                                       1063
      kcp0=kcp0+kcp                                                      1064
      tb(4,m)=tx(4)                                                      1065
      k1=kr                                                              1066
      call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))                            1067
      tx(4)=m                                                            1068
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  1069
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    1070
      call addpar(m)                                                     1071
      if(m .ge. mk) go to 58                                             1072
      m=m+1                                                              1073
      tb(2,m)=-tb(2,m-1)                                                 1074
      do 55 i=3,4                                                        1075
      tb(i,m)=tb(i,m-1)                                                  1076
   55 continue                                                           1077
      if(ibfext(m,tb,cm) .eq. 0) go to 56                                1078
      m=m-1                                                              1079
      go to 58                                                           1080
   56 do 57 i=1,n                                                        1081
      sc(i,m)=phi(m,i,n,x,tb,cm)                                         1082
   57 continue                                                           1083
      call addpar(m)                                                     1084
   58 if(it .le. 0) go to 59                                             1085
      mp=m-1                                                             1086
      tcst=(nopt-1)*df1+kr+1.0                                           1087
      fjn=jn                                                             1088
      fkr=kr                                                             1089
      gcv=(rsq/sw)/(1.d0-tcst/wn)**2                                     1090
      call holl(jn,cm,tb(3,m),hol)                                       1091
c     if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,fjn,hol,tb(4,m)       1092
c     if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,fjn,hol,tb(4,m)     1093
   59 mtot=m                                                             1094
      m=m+1                                                              1095
      if(m.gt.mk) go to 69                                               1096
   60 do 61 i=1,5                                                        1097
      tb(i,m)=tx(i)                                                      1098
   61 continue                                                           1099
      k1=kr                                                              1100
      call blf(int(tx(4)+.1),n,sc,sc(1,mkp1))                            1101
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  1102
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    1103
      call addpar(m)                                                     1104
      if(m .ge. mk .or. (cm(2*jp) .le. 0.0) .and. (tx(3) .le. x(mm(1,jp) 1105
     1,jp))) go to 66                                                    1106
      m=m+1                                                              1107
      do 62 i=1,4                                                        1108
      tb(i,m)=tx(i)                                                      1109
   62 continue                                                           1110
      tb(2,m)=-tb(2,m)                                                   1111
      if(cm(2*jp) .le. 0.0) go to 64                                     1112
      do 63 i=1,n                                                        1113
      sc(i,m)=phi(m,i,n,x,tb,cm)                                         1114
   63 continue                                                           1115
      go to 65                                                           1116
   64 k1=kr                                                              1117
      call update(2,n,m,kr,x,y,w,sw,yb,tb,cm,sc,sc(1,mkp1),db,d,d(1,3))  1118
      if(kr.gt.k1) rsq=rsq-d(kr,1)**2                                    1119
   65 call addpar(m)                                                     1120
   66 tcst=nopt*df1+kr+1.0                                               1121
      if(it .le. 0) go to 68                                             1122
      mp=m-1                                                             1123
      jp=abs(tx(2))+.1                                                   1124
      fkr=kr                                                             1125
      gcv=(rsq/sw)/(1.d0-tcst/wn)**2                                     1126
      if(cm(2*jp) .le. 0.0) go to 67                                     1127
      call holl(jp,cm,tx(3),hol)                                         1128
c     if(m.eq.mtot+1) write(it,100) m,gcv,fkr,tcst,tx(2),hol,tx(4)       1129
c     if(m.eq.mtot+2) write(it,99) m,mp,gcv,fkr,tcst,tx(2),hol,tx(4)     1130
      go to 68                                                           1131
   67 xk=xm(jp)+xs(jp)*tx(3)                                             1132
c     if(m.eq.mtot+1) write(it,93) m,gcv,fkr,tcst,tx(2),xk,tx(4)         1133
c     if(m.eq.mtot+2) write(it,94) m,mp,gcv,fkr,tcst,tx(2),xk,tx(4)      1134
   68 mtot=m                                                             1135
      go to 8                                                            1136
   69 mk=min0(m,mk)                                                      1137
      m=mk+1                                                             1138
      k=m                                                                1139
      go to 71                                                           1140
   70 k=k+1                                                              1141
   71 if((k).gt.(nk)) go to 72                                           1142
      tb(1,k)=0.0                                                        1143
      go to 70                                                           1144
   72 call sscp(n,m,sc,y,w,yb,yv,sw,db,d)                                1145
      call lsf1(db,m,d,yb,alr,b,d(1,2),a,d(1,3))                         1146
      nli=0                                                              1147
      do 73 k=1,mk                                                       1148
      if(d(k,2).ne.0.d0) nli=nli+1                                       1149
   73 continue                                                           1150
      df1=df1*nopt+nli                                                   1151
      tcst=df1+1.0                                                       1152
      df1=df1/nli                                                        1153
      do 74 k=1,nk                                                       1154
      tb(5,k)=df1                                                        1155
   74 continue                                                           1156
      asm=(b/sw)/(1.d0-tcst/wn)**2                                       1157
      tcsts=tcst                                                         1158
      az=a                                                               1159
      do 75 k=1,mk                                                       1160
      tb(1,k)=0.0                                                        1161
      if(d(k,2).ne.0.d0) tb(1,k)=d(k,2)                                  1162
   75 continue                                                           1163
      if(ix .eq. 0) go to 81                                             1164
      sc(1,1)=(cfac*nopt)/nli                                            1165
      sc(2,1)=wn                                                         1166
      sc(3,1)=yv                                                         1167
      sc(4,1)=yb                                                         1168
      do 80 k=nli,nk                                                     1169
      call array(k+4,n,i,j)                                              1170
      sc(i,j)=b/sw                                                       1171
      k1=k*(nk+1)+3                                                      1172
      l=0                                                                1173
      go to 77                                                           1174
   76 l=l+1                                                              1175
   77 if((l).gt.(nk)) go to 80                                           1176
      k1=k1+1                                                            1177
      call array(k1,n,i,j)                                               1178
      if(l .ne. 0) go to 78                                              1179
      sc(i,j)=a                                                          1180
      go to 76                                                           1181
   78 if(l .le. mk) go to 79                                             1182
      sc(i,j)=0.0                                                        1183
      go to 76                                                           1184
   79 sc(i,j)=d(l,2)                                                     1185
      go to 76                                                           1186
   80 continue                                                           1187
      call array((nk+1)**2+4,n,i,j)                                      1188
      sc(i,j)=mk                                                         1189
      kl=nli                                                             1190
   81 do 88 ll=2,nli                                                     1191
      call bkstp(db,m,d,yb,alr,b,d(1,2),a,k,d(1,3))                      1192
      if(k.eq.0) go to 89                                                1193
      if(ix .eq. 0) go to 86                                             1194
      call array(kl+3,n,i,j)                                             1195
      sc(i,j)=b/sw                                                       1196
      kl=kl-1                                                            1197
      k1=kl*(nk+1)+3                                                     1198
      l=0                                                                1199
      go to 83                                                           1200
   82 l=l+1                                                              1201
   83 if((l).gt.(nk)) go to 86                                           1202
      k1=k1+1                                                            1203
      call array(k1,n,i,j)                                               1204
      if(l .ne. 0) go to 84                                              1205
      sc(i,j)=a                                                          1206
      go to 82                                                           1207
   84 if(l .le. mk) go to 85                                             1208
      sc(i,j)=0.0                                                        1209
      go to 82                                                           1210
   85 sc(i,j)=d(l,2)                                                     1211
      go to 82                                                           1212
   86 tcst=tcst-df1                                                      1213
      b=(b/sw)/(1.d0-tcst/wn)**2                                         1214
      if(b.ge.asm) go to 88                                              1215
      asm=b                                                              1216
      tcsts=tcst                                                         1217
      az=a                                                               1218
      do 87 i=1,mk                                                       1219
      tb(1,i)=0.0                                                        1220
      if(d(i,2).ne.0.d0) tb(1,i)=d(i,2)                                  1221
   87 continue                                                           1222
   88 continue                                                           1223
   89 if(txm .gt. asm) go to 91                                          1224
      asm=txm                                                            1225
      tcsts=1.0                                                          1226
      az=yb                                                              1227
      do 90 i=1,nk                                                       1228
      tb(1,i)=0.0                                                        1229
   90 continue                                                           1230
   91 if(it .le. 0) go to 92                                             1231
c     write(it,95)                                                       1232
      call coefpr(it,mk,az,tb,cm,xs)                                     1233
c     write(it,96) asm,tcsts                                             1234
   92 return                                                             1235
      entry setfln(val)                                                  1236
      fln=val                                                            1237
      return                                                             1238
      entry setalf(val)                                                  1239
      alf=val                                                            1240
      return                                                             1241
      entry setmin(nal)                                                  1242
      nmin=nal                                                           1243
      return                                                             1244
      entry setcta(val)                                                  1245
      vcst(2)=val                                                        1246
      return                                                             1247
      entry setctl(val)                                                  1248
      vcst(3)=val                                                        1249
      return                                                             1250
      entry xvmrgo(nal)                                                  1251
      ix=nal                                                             1252
      return                                                             1253
      entry setalr(val)                                                  1254
      alr=val                                                            1255
      return                                                             1256
   93 format('   ',i3,'    ',g12.4,2('   ',f5.1),'       ',f3.0,  '      1257
     1 ',g12.4,'          ',f3.0)                                        1258
   94 format(' ',i3,' ',i3,'  ',g12.4,2('   ',f5.1),'       ',f3.0,  '   1259
     1    ',g12.4,'          ',f3.0)                                     1260
   95 format(/,' final model after backward stepwise elimination:')      1261
   96 format(/,'   (piecewise linear) gcv = ',g12.4,'   #efprms = ',f5.1 1262
     1)                                                                  1263
   97 format(//,' forward stepwise knot placement:',//  '  basfn(s)    g 1264
     1cv      #indbsfns  #efprms',  '   variable      knot            pa 1265
     1rent')                                                             1266
   98 format('   ',i3,'    ',g12.4,2('   ',f5.1))                        1267
   99 format(' ',i3,' ',i3,'  ',g12.4,2('   ',f5.1),'       ',f3.0,a28,f 1268
     13.0)                                                               1269
  100 format('   ',i3,'    ',g12.4,2('   ',f5.1),'       ',f3.0,a28,f3.0 1270
     1)                                                                  1271
      end                                                                1272
      subroutine addpar (ib)                                             1273
      parameter(maxdph=1000)                                             1274
      real que(2,maxdph),sp(maxdph)                                      1275
      integer m(maxdph),n(maxdph),jp(2,maxdph)                           1276
      double precision val                                               1277
      save mpr,mtr,ktr,big,beta,lq,kp,itr,jp,que,n,m                     1278
      data big,mpr,mtr,beta /9.9e30,10,5,1.0/                            1279
      if(ib .ne. 0) go to 1                                              1280
      lq=1                                                               1281
      que(1,1)=big                                                       1282
      que(2,1)=0.0                                                       1283
      m(1)=1                                                             1284
      kp=0                                                               1285
      itr=kp                                                             1286
      n(1)=itr                                                           1287
      jp(1,1)=n(1)                                                       1288
      jp(2,1)=jp(1,1)                                                    1289
      ktr=0.5*(mpr-1)+.1                                                 1290
      return                                                             1291
    1 if(que(1,lq).ge.-0.5) go to 2                                      1292
      lq=lq-1                                                            1293
      go to 1                                                            1294
    2 i=1                                                                1295
      go to 4                                                            1296
    3 i=i+1                                                              1297
    4 if((i).gt.(lq)) go to 7                                            1298
      if(que(1,i).ge.-0.5) go to 3                                       1299
      lq=lq-1                                                            1300
      do 6 j=i,lq                                                        1301
      n(j)=n(j+1)                                                        1302
      do 5 k=1,2                                                         1303
      jp(k,j)=jp(k,j+1)                                                  1304
      que(k,j)=que(k,j+1)                                                1305
    5 continue                                                           1306
    6 continue                                                           1307
      i=i-1                                                              1308
      go to 3                                                            1309
    7 lq=lq+1                                                            1310
      if(lq .le. maxdph) go to 8                                         1311
c     write(6, '('' increase parameter maxdph in subroutine addpar to '' 1312
c    1)')                                                                1313
c     write(6,'('' '',i10,''  or larger, and recompile mars.'')') lq     1314
      stop                                                               1315
    8 que(1,lq)=big                                                      1316
      que(2,lq)=0.0                                                      1317
      n(lq)=ib                                                           1318
      jp(1,lq)=0                                                         1319
      jp(2,lq)=jp(1,lq)                                                  1320
      do 9 i=1,lq                                                        1321
      m(i)=i                                                             1322
      sp(i)=que(1,i)                                                     1323
    9 continue                                                           1324
      call psort(sp,m,1,lq)                                              1325
      do 10 i=1,lq                                                       1326
      j=m(i)                                                             1327
      sp(j)=i+beta*(itr-que(2,j))                                        1328
   10 continue                                                           1329
      call psort(sp,m,1,lq)                                              1330
      kp=max0(0,lq-mpr)                                                  1331
      return                                                             1332
      entry nxtpar (l,jq)                                                1333
      kp=kp+1                                                            1334
      if(kp .le. lq) go to 11                                            1335
      l=-1                                                               1336
      return                                                             1337
   11 l=n(m(kp))                                                         1338
      if(itr-jp(2,m(kp)).gt.mtr.or.itr.le.ktr) jp(1,m(kp))=0             1339
      jq=jp(1,m(kp))                                                     1340
      return                                                             1341
      entry updpar (jq,val)                                              1342
      que(1,m(kp))=val                                                   1343
      que(2,m(kp))=itr                                                   1344
      if(jp(1,m(kp)) .ne. 0) go to 12                                    1345
      jp(1,m(kp))=jq                                                     1346
      jp(2,m(kp))=itr                                                    1347
   12 return                                                             1348
      entry selpar(ib)                                                   1349
      do 13 i=lq,1,-1                                                    1350
      if(n(i).ne.ib) go to 13                                            1351
      jp(1,i)=0                                                          1352
      go to 14                                                           1353
   13 continue                                                           1354
   14 return                                                             1355
      entry itrpar(iarg)                                                 1356
      itr=iarg                                                           1357
      return                                                             1358
      entry setmpr(iarg)                                                 1359
      mpr=iarg                                                           1360
      return                                                             1361
      entry setbta(arg)                                                  1362
      beta=arg                                                           1363
      return                                                             1364
      entry setfrq(arg)                                                  1365
      mtr=1.0/amax1(arg,0.01)+.1                                         1366
      return                                                             1367
      end                                                                1368
      subroutine speed(is)                                               1369
      integer lque(5)                                                    1370
      real freq(5)                                                       1371
      save lque,freq                                                     1372
      data lque /9999,20,20,10,5/                                        1373
      data freq /9.e30,9.e30,0.2,0.2,0.2/                                1374
      j=is                                                               1375
      if(is.lt.1) j=1                                                    1376
      if(is.gt.5) j=5                                                    1377
      call setmpr(lque(j))                                               1378
      call setfrq(freq(j))                                               1379
      return                                                             1380
      end                                                                1381
      subroutine atoscl(n,p,w,x,lx,mm,xm,xs,cm,z)                        1382
      integer p,lx(p),mm(n,p)                                            1383
      real w(n),x(n,p),z(n,p),xm(p),xs(p),cm(*)                          1384
      double precision s,t,sw                                            1385
      sw=0.d0                                                            1386
      do 1 j=1,n                                                         1387
      sw=sw+w(j)                                                         1388
    1 continue                                                           1389
      ip=0                                                               1390
      nct=ip                                                             1391
      do 12 i=1,p                                                        1392
      if(lx(i) .ne. 0) go to 2                                           1393
      xm(i)=0.0                                                          1394
      xs(i)=xm(i)                                                        1395
      go to 12                                                           1396
    2 if(lx(i) .ge. 0) go to 8                                           1397
      nc=0                                                               1398
      xm(i)=ip                                                           1399
      j=1                                                                1400
      nct=nct+1                                                          1401
    3 j0=j                                                               1402
      if(j .ge. n) go to 5                                               1403
    4 if(x(mm(j+1,i),i).gt.x(mm(j,i),i)) go to 5                         1404
      j=j+1                                                              1405
      if(j.ge.n) go to 5                                                 1406
      go to 4                                                            1407
    5 ip=ip+1                                                            1408
      cm(ip)=x(mm(j,i),i)                                                1409
      nc=nc+1                                                            1410
      do 6 k=j0,j                                                        1411
      z(mm(k,i),i)=nc                                                    1412
    6 continue                                                           1413
      j=j+1                                                              1414
      if(j.gt.n) go to 7                                                 1415
      go to 3                                                            1416
    7 xs(i)=nc                                                           1417
      go to 12                                                           1418
    8 s=0.d0                                                             1419
      t=s                                                                1420
      do 9 j=1,n                                                         1421
      s=s+w(j)*x(j,i)                                                    1422
    9 continue                                                           1423
      s=s/sw                                                             1424
      xm(i)=s                                                            1425
      do 10 j=1,n                                                        1426
      z(j,i)=x(j,i)-s                                                    1427
      t=t+w(j)*z(j,i)**2                                                 1428
   10 continue                                                           1429
      xs(i)=1.0                                                          1430
      if(t.le.0.d0) go to 12                                             1431
      t=dsqrt(t/sw)                                                      1432
      xs(i)=t                                                            1433
      t=1.d0/t                                                           1434
      do 11 j=1,n                                                        1435
      z(j,i)=t*z(j,i)                                                    1436
   11 continue                                                           1437
   12 continue                                                           1438
      n2=2*p+1                                                           1439
      if(nct .ne. 0) go to 14                                            1440
      do 13 i=1,n2                                                       1441
      cm(i)=0.0                                                          1442
   13 continue                                                           1443
      return                                                             1444
   14 n2p1=n2+1                                                          1445
      i=ip                                                               1446
      go to 16                                                           1447
   15 i=i+(-1)                                                           1448
   16 if((-1)*((i)-(1)).gt.0) go to 17                                   1449
      cm(i+n2)=cm(i)                                                     1450
      go to 15                                                           1451
   17 j=0                                                                1452
      i=2                                                                1453
      go to 19                                                           1454
   18 i=i+(2)                                                            1455
   19 if((2)*((i)-(n2)).gt.0) go to 22                                   1456
      j=j+1                                                              1457
      if(lx(j) .ge. 0) go to 20                                          1458
      cm(i)=xm(j)+n2p1                                                   1459
      cm(i+1)=cm(i)+xs(j)-1.0                                            1460
      go to 18                                                           1461
   20 cm(i)=0.0                                                          1462
      cm(i+1)=cm(i)                                                      1463
      go to 18                                                           1464
   22 cm(1)=nct                                                          1465
      call stfmrs(1)                                                     1466
      call stcmrs(1)                                                     1467
      return                                                             1468
      end                                                                1469
      subroutine orgpl (xm,xs,nk,tb,cm)                                  1470
      real xm(*),xs(*),tb(5,nk),cm(*)                                    1471
      do 1 m=1,nk                                                        1472
      j=abs(tb(2,m))+.1                                                  1473
      if(cm(2*j).gt.0.0) go to 1                                         1474
      tb(3,m)=xm(j)+xs(j)*tb(3,m)                                        1475
    1 continue                                                           1476
      do 4 m=1,nk                                                        1477
      if(tb(1,m).eq.0.0) go to 4                                         1478
      scl=1.0                                                            1479
      ip=m                                                               1480
    2 if(ip.le.0) go to 3                                                1481
      j=abs(tb(2,ip))+.1                                                 1482
      if(cm(2*j).eq.0.0) scl=scl*xs(j)                                   1483
      ip=tb(4,ip)+.1                                                     1484
      go to 2                                                            1485
    3 tb(1,m)=tb(1,m)/scl                                                1486
    4 continue                                                           1487
      return                                                             1488
      end                                                                1489
      subroutine anova (n,x,y,w,nk,it,tb,cm,lp,lv,t,d)                   1490
      integer lp(3,*),lv(*)                                              1491
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),t(n,nk)                       1492
      double precision d(nk,*),s,u,sw,yv,wn,yb                           1493
      if(it.le.0) return                                                 1494
      nkp1=nk+1                                                          1495
      nkp2=nkp1+1                                                        1496
      nkp3=nkp2+1                                                        1497
      lm=nkp3+nk                                                         1498
      sw=0.d0                                                            1499
      wn=sw                                                              1500
      s=wn                                                               1501
      u=s                                                                1502
      do 1 i=1,n                                                         1503
      sw=sw+w(i)                                                         1504
      wn=wn+w(i)**2                                                      1505
      s=s+w(i)*y(i)                                                      1506
    1 continue                                                           1507
      s=s/sw                                                             1508
      yb=s                                                               1509
      wn=sw**2/wn                                                        1510
      do 2 i=1,n                                                         1511
      u=u+w(i)*(y(i)-s)**2                                               1512
    2 continue                                                           1513
      yv=u/sw                                                            1514
      eft=1.0                                                            1515
      do 3 m=1,nk                                                        1516
      if(tb(1,m).ne.0.0) eft=eft+tb(5,m)                                 1517
    3 continue                                                           1518
      ni=0                                                               1519
      do 9 m=1,nk                                                        1520
      if(tb(1,m).eq.0.0) go to 9                                         1521
      ni=ni+1                                                            1522
      s=0.d0                                                             1523
      do 4 j=1,n                                                         1524
      t(j,ni)=phi(m,j,n,x,tb,cm)                                         1525
      s=s+w(j)*t(j,ni)                                                   1526
    4 continue                                                           1527
      s=s/sw                                                             1528
      do 5 j=1,n                                                         1529
      t(j,ni)=t(j,ni)-s                                                  1530
    5 continue                                                           1531
      do 7 i=1,ni                                                        1532
      s=0.d0                                                             1533
      do 6 j=1,n                                                         1534
      s=s+w(j)*t(j,i)*t(j,ni)                                            1535
    6 continue                                                           1536
      d(i,ni)=s                                                          1537
    7 continue                                                           1538
      s=0.d0                                                             1539
      do 8 j=1,n                                                         1540
      s=s+w(j)*t(j,ni)*(y(j)-yb)                                         1541
    8 continue                                                           1542
      d(ni,nkp1)=s                                                       1543
      d(ni,nkp2)=tb(1,m)                                                 1544
      lv(ni)=m                                                           1545
    9 continue                                                           1546
      if(ni .ne. 0) go to 10                                             1547
c     write(it,26)                                                       1548
      return                                                             1549
   10 do 11 m=1,ni                                                       1550
      t(m,1)=lv(m)                                                       1551
   11 continue                                                           1552
c     write(it,24) ni                                                    1553
      call coll(nk,tb,lp,lv,lp(1,nkp1))                                  1554
      m=1                                                                1555
   12 if(lp(1,m).eq.0) go to 13                                          1556
      m=m+1                                                              1557
      go to 12                                                           1558
   13 na=m-1                                                             1559
      m=1                                                                1560
      nim1=ni-1                                                          1561
      if(na .ne. 1) go to 14                                             1562
      k2=lp(2,m)                                                         1563
      i2=lp(1,m)+k2-1                                                    1564
      efm=eft-1.0                                                        1565
      u=yv/(1.d0-1.d0/wn)**2                                             1566
      s=sqrt(varf(nk,d,d(1,nkp2),sw,1,ni))                               1567
c     write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)                     1568
      return                                                             1569
   14 do 23 m=1,na                                                       1570
      k2=lp(2,m)                                                         1571
      l=lp(1,m)                                                          1572
      i2=l+k2-1                                                          1573
      ll=k2-1                                                            1574
      np=ni                                                              1575
      do 19 im=1,ni                                                      1576
      i=t(im,1)+.1                                                       1577
      if(nord(i,tb) .eq. l) go to 15                                     1578
      t(im,2)=0.0                                                        1579
      go to 19                                                           1580
   15 k=0                                                                1581
      do 16 j=1,l                                                        1582
      if(jf(i,lv(ll+j),tb).eq.1) go to 16                                1583
      k=1                                                                1584
      go to 17                                                           1585
   16 continue                                                           1586
   17 if(k .ne. 1) go to 18                                              1587
      t(im,2)=0.0                                                        1588
      go to 19                                                           1589
   18 t(im,2)=1.0                                                        1590
      np=np-1                                                            1591
   19 continue                                                           1592
   20 k=0                                                                1593
      do 21 i=1,nim1                                                     1594
      if(t(i,2) .le. t(i+1,2)) go to 21                                  1595
      k=1                                                                1596
      call exch(nk,ni,i,d,t,t(1,2))                                      1597
   21 continue                                                           1598
      if(k.eq.0) go to 22                                                1599
      go to 20                                                           1600
   22 call lsf(nk,np,nkp1,0.d0,d,d(1,lm),s,u,d(1,nkp3),1)                1601
      efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)                                 1602
      u=(u/sw+yv)/(1.d0-(eft-efm)/wn)**2                                 1603
      s=sqrt(varf(nk,d,d(1,nkp2),sw,np+1,ni))                            1604
c     write(it,25) m,s,u,lp(3,m),efm,(lv(i),i=k2,i2)                     1605
   23 continue                                                           1606
      return                                                             1607
   24 format(/' anova decomposition on',i3,' basis functions:',/  '  fun 1608
     1. std. dev.     -gcv    #bsfns  #efprms  variable(s)')             1609
   25 format(' ',i3,' ',2g12.4,'  ',i2,'      ',f4.1,'  ',20i4)          1610
   26 format(/' estimated optimal model = response mean.')               1611
      end                                                                1612
      subroutine anoval (n,x,y,w,nk,il,it,az,tb,cm,lp,lv,sc,d)           1613
      integer lp(3,*),lv(*)                                              1614
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*)                       1615
      double precision d(nk,*),sw,yv,wn,yb                               1616
      if(it.le.0) return                                                 1617
      sw=0.d0                                                            1618
      yb=sw                                                              1619
      yv=yb                                                              1620
      wn=yv                                                              1621
      do 1 i=1,n                                                         1622
      sw=sw+w(i)                                                         1623
      wn=wn+w(i)**2                                                      1624
      yb=yb+w(i)*y(i)                                                    1625
    1 continue                                                           1626
      yb=yb/sw                                                           1627
      wn=sw**2/wn                                                        1628
      do 2 i=1,n                                                         1629
      yv=yv+w(i)*(y(i)-yb)**2                                            1630
    2 continue                                                           1631
      yv=yv/sw                                                           1632
      eft=1.0                                                            1633
      ni=0                                                               1634
      do 3 m=1,nk                                                        1635
      if(tb(1,m).eq.0.0) go to 3                                         1636
      ni=ni+1                                                            1637
      eft=eft+tb(5,m)                                                    1638
    3 continue                                                           1639
      if(ni .ne. 0) go to 4                                              1640
c     write(it,14)                                                       1641
      return                                                             1642
    4 continue
c     write(it,12) ni                                                    1643
      call coll(nk,tb,lp,lv,lp(1,nk+1))                                  1644
      m=1                                                                1645
    5 if(lp(1,m).eq.0) go to 6                                           1646
      m=m+1                                                              1647
      go to 5                                                            1648
    6 na=m-1                                                             1649
      m=1                                                                1650
      if(na .ne. 1) go to 7                                              1651
      k2=lp(2,m)                                                         1652
      i2=lp(1,m)+k2-1                                                    1653
      efm=eft-1.0                                                        1654
      u=yv/(1.d0-1.d0/wn)**2                                             1655
c     write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)                       1656
      return                                                             1657
    7 ip=nk+4                                                            1658
      do 11 m=1,na                                                       1659
      k2=lp(2,m)                                                         1660
      l=lp(1,m)                                                          1661
      i2=l+k2-1                                                          1662
      ll=k2-1                                                            1663
      call cptb(nk,tb,sc(1,ip))                                          1664
      do 10 i=1,nk                                                       1665
      if(tb(1,i).eq.0.0) go to 10                                        1666
      if(nord(i,tb).ne.l) go to 10                                       1667
      k=0                                                                1668
      do 8 j=1,l                                                         1669
      if(jf(i,lv(ll+j),tb).eq.1) go to 8                                 1670
      k=1                                                                1671
      go to 9                                                            1672
    8 continue                                                           1673
    9 if(k.eq.1) go to 10                                                1674
      call setz(i,sc(1,ip))                                              1675
   10 continue                                                           1676
      a0=az                                                              1677
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,u,sc,d)                 1678
      efm=efp(lp(1,m),lv(lp(2,m)),nk,tb)                                 1679
      u=u/(1.d0-(eft-efm)/wn)**2                                         1680
c     write(it,13) m,u,lp(3,m),efm,(lv(i),i=k2,i2)                       1681
   11 continue                                                           1682
      return                                                             1683
   12 format(/' logit anova decomposition on',i3,' basis functions:',/   1684
     1'  fun.    -gcv    #bsfns  #efprms  variable(s)')                  1685
   13 format(' ',i3,' ',g12.4,'   ',i2,'     ',f4.1,'    ',20i4)         1686
   14 format(/' estimated optimal model = response mean.')               1687
      end                                                                1688
      subroutine cptb(nk,tb,ub)                                          1689
      real tb(5,nk),ub(5,*)                                              1690
      do 2 m=1,nk                                                        1691
      do 1 k=1,5                                                         1692
      ub(k,m)=tb(k,m)                                                    1693
    1 continue                                                           1694
    2 continue                                                           1695
      return                                                             1696
      entry setz(l,ub)                                                   1697
      ub(1,l)=0.0                                                        1698
      return                                                             1699
      end                                                                1700
      subroutine fun (l,jv,n,x,nk,tb,cm,jl,kv,t,js)                      1701
      integer jv(l),kv(2,jl),js(*)                                       1702
      real x(n,l),tb(5,nk),cm(*),t(n)                                    1703
      double precision s                                                 1704
      do 8 i=1,n                                                         1705
      s=0.d0                                                             1706
      do 7 m=1,nk                                                        1707
      if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7                             1708
      if(nordc(1,m,tb,cm).ne.l) go to 7                                  1709
      k=0                                                                1710
      do 1 j=1,l                                                         1711
      if(jf(m,jv(j),tb).eq.1) go to 1                                    1712
      k=1                                                                1713
      go to 2                                                            1714
    1 continue                                                           1715
    2 if(k.eq.1) go to 7                                                 1716
      phi=1.0                                                            1717
      ip=m                                                               1718
    3 if(ip.le.0) go to 6                                                1719
      u=tb(2,ip)                                                         1720
      j=abs(u)+.1                                                        1721
      if(cm(2*j) .eq. 0.0) go to 4                                       1722
      ip=tb(4,ip)+.1                                                     1723
      go to 3                                                            1724
    4 do 5 k=1,l                                                         1725
      if(j.eq.jv(k)) j=k                                                 1726
    5 continue                                                           1727
      phi=phi*amax1(0.0,sign(1.0,u)*(x(i,j)-tb(3,ip)))                   1728
      ip=tb(4,ip)+.1                                                     1729
      go to 3                                                            1730
    6 s=s+tb(1,m)*phi                                                    1731
    7 continue                                                           1732
      t(i)=s                                                             1733
    8 continue                                                           1734
      return                                                             1735
      end                                                                1736
      subroutine cubic (n,p,x,y,w,nk,it,tb,cm,kp,kv,lp,lv,bz,tc,t,z,sc,j 1737
     1s,d)                                                               1738
      integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*)                      1739
      real x(n,p),y(n),w(n),tb(5,nk),cm(*),tc(*),t(n,nk),z(2,p),sc(n)    1740
      double precision d(nk,*),s,u,sw,yb,wn,yv                           1741
      data big /9.9e30/                                                  1742
      yb=0.d0                                                            1743
      sw=yb                                                              1744
      wn=sw                                                              1745
      yv=wn                                                              1746
      do 1 i=1,n                                                         1747
      sw=sw+w(i)                                                         1748
      wn=wn+w(i)**2                                                      1749
      yb=yb+w(i)*y(i)                                                    1750
    1 continue                                                           1751
      yb=yb/sw                                                           1752
      wn=sw**2/wn                                                        1753
      do 2 i=1,n                                                         1754
      yv=yv+w(i)*(y(i)-yb)**2                                            1755
    2 continue                                                           1756
      yv=yv/sw                                                           1757
      ni=0                                                               1758
      do 3 m=1,nk                                                        1759
      if(tb(1,m).ne.0.0) ni=ni+1                                         1760
    3 continue                                                           1761
      if(ni .ne. 0) go to 4                                              1762
      bz=yb                                                              1763
      u=yv/(1.0-1.0/wn)**2                                               1764
c     if(it.gt.0) write(it,34) ni,u                                      1765
      return                                                             1766
    4 nkp1=nk+1                                                          1767
      nkp2=nk+2                                                          1768
      nkp3=nk+3                                                          1769
      lm=nkp3+nk                                                         1770
      do 6 i=1,p                                                         1771
      xl=big                                                             1772
      xr=-xl                                                             1773
      do 5 j=1,n                                                         1774
      xl=amin1(xl,x(j,i))                                                1775
      xr=amax1(xr,x(j,i))                                                1776
    5 continue                                                           1777
      z(1,i)=xl                                                          1778
      z(2,i)=xr                                                          1779
    6 continue                                                           1780
      ll=1                                                               1781
      la=ll                                                              1782
      l1=la                                                              1783
      lt=0                                                               1784
    7 if(kp(1,ll).lt.0) go to 20                                         1785
      do 8 i=1,n                                                         1786
      sc(i)=1.0                                                          1787
    8 continue                                                           1788
      if(kp(1,ll) .le. 0) go to 12                                       1789
      jl=kp(1,ll)                                                        1790
      do 11 il=1,jl                                                      1791
      k=kp(2,ll)+il-1                                                    1792
      jj=kv(1,k)                                                         1793
      j=iabs(jj)                                                         1794
      kk=kv(2,k)                                                         1795
      do 10 i=1,n                                                        1796
      if(sc(i).eq.0.0) go to 10                                          1797
      ic=x(i,j)+.1                                                       1798
      sc(i)=cm(ic+kk)                                                    1799
      if(jj .ge. 0) go to 10                                             1800
      if(sc(i) .ne. 0.0) go to 9                                         1801
      sc(i)=1.0                                                          1802
      go to 10                                                           1803
    9 sc(i)=0.0                                                          1804
   10 continue                                                           1805
   11 continue                                                           1806
      go to 13                                                           1807
   12 if(kp(3,ll) .gt. 0) go to 13                                       1808
      ll=ll+1                                                            1809
      go to 7                                                            1810
   13 if(kp(3,ll) .gt. 0) go to 15                                       1811
      lt=lt+1                                                            1812
      kp(5,ll)=0                                                         1813
      do 14 i=1,n                                                        1814
      t(i,lt)=sc(i)                                                      1815
   14 continue                                                           1816
      go to 19                                                           1817
   15 kp3=kp(3,ll)                                                       1818
      kp(5,ll)=la                                                        1819
      do 18 m=1,kp3                                                      1820
      l=lp(1,l1)                                                         1821
      nt=lp(3,l1)                                                        1822
      call knts(l,nt,lv(lp(2,l1)),kp(1,ll),kv(1,kp(2,ll)),  nk,tb,cm,tc( 1823
     1la),js)                                                            1824
      call side(l,nt,lv(lp(2,l1)),z,tc(la))                              1825
      do 17 jp=1,nt                                                      1826
      lt=lt+1                                                            1827
      do 16 i=1,n                                                        1828
      t(i,lt)=sc(i)                                                      1829
   16 continue                                                           1830
      call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),t(1,lt))                  1831
   17 continue                                                           1832
      l1=l1+1                                                            1833
      la=la+nt*(5*l+1)                                                   1834
   18 continue                                                           1835
   19 ll=ll+1                                                            1836
      go to 7                                                            1837
   20 do 26 j=1,lt                                                       1838
      s=0.d0                                                             1839
      u=s                                                                1840
      do 21 i=1,n                                                        1841
      s=s+w(i)*t(i,j)                                                    1842
   21 continue                                                           1843
      s=s/sw                                                             1844
      d(j,nkp2)=s                                                        1845
      do 22 i=1,n                                                        1846
      t(i,j)=t(i,j)-s                                                    1847
   22 continue                                                           1848
      s=0.d0                                                             1849
      do 23 i=1,n                                                        1850
      s=s+w(i)*(y(i)-yb)*t(i,j)                                          1851
   23 continue                                                           1852
      d(j,nkp1)=s                                                        1853
      do 25 k=1,j                                                        1854
      s=0.d0                                                             1855
      do 24 i=1,n                                                        1856
      s=s+w(i)*t(i,k)*t(i,j)                                             1857
   24 continue                                                           1858
      d(k,j)=s                                                           1859
   25 continue                                                           1860
   26 continue                                                           1861
      call lsf(nk,lt,nkp1,yb,d,d(1,lm),s,u,d(1,nkp3),1)                  1862
      eft=1.0                                                            1863
      do 27 i=1,nk                                                       1864
      if(tb(1,i).ne.0.0) eft=eft+tb(5,i)                                 1865
   27 continue                                                           1866
      u=(u/sw+yv)/(1.0-eft/wn)**2                                        1867
      bz=s                                                               1868
      ll=1                                                               1869
      l1=ll                                                              1870
      le=la-1                                                            1871
      la=0                                                               1872
      lt=la                                                              1873
   28 if(kp(1,ll).lt.0) go to 33                                         1874
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 29                  1875
      ll=ll+1                                                            1876
      go to 28                                                           1877
   29 if(kp(3,ll) .gt. 0) go to 30                                       1878
      le=le+1                                                            1879
      kp(3,ll)=-le                                                       1880
      lt=lt+1                                                            1881
      tc(le)=d(lt,lm)                                                    1882
      ll=ll+1                                                            1883
      go to 28                                                           1884
   30 kp3=kp(3,ll)                                                       1885
      do 32 m=1,kp3                                                      1886
      nt=lp(3,l1)                                                        1887
      la=la+5*lp(1,l1)*nt                                                1888
      do 31 i=1,nt                                                       1889
      lt=lt+1                                                            1890
      tc(i+la)=d(lt,lm)                                                  1891
   31 continue                                                           1892
      la=la+nt                                                           1893
      l1=l1+1                                                            1894
   32 continue                                                           1895
      ll=ll+1                                                            1896
      go to 28                                                           1897
   33 continue
c     if(it.gt.0) write(it,34) lt,u                                      1898
      return                                                             1899
   34 format(/' piecewise cubic fit on',i3,' basis functions, gcv =',g12 1900
     1.4)                                                                1901
      end                                                                1902
      subroutine cfun (l,jv,n,x,nf,lp,lv,tc,t,sc,jw)                     1903
      integer jv(l),lp(3,*),lv(*),jw(l)                                  1904
      real x(n,l),tc(*),t(n),sc(n)                                       1905
      do 1 i=1,n                                                         1906
      t(i)=0.0                                                           1907
    1 continue                                                           1908
      la=1                                                               1909
      do 10 l1=1,nf                                                      1910
      if(lp(1,l1).ne.l) go to 9                                          1911
      l2=lp(2,l1)-1                                                      1912
      do 3 j=1,l                                                         1913
      m=0                                                                1914
      do 2 k=1,l                                                         1915
      if(jv(j).eq.lv(k+l2)) m=1                                          1916
    2 continue                                                           1917
      if(m.eq.0) go to 9                                                 1918
    3 continue                                                           1919
      nt=lp(3,l1)                                                        1920
      lb=la+5*l*nt-1                                                     1921
      do 8 j=1,nt                                                        1922
      do 5 k=1,l                                                         1923
      do 4 i=1,l                                                         1924
      if(lv(k+l2).eq.jv(i)) jw(k)=i                                      1925
    4 continue                                                           1926
    5 continue                                                           1927
      do 6 i=1,n                                                         1928
      sc(i)=1.0                                                          1929
    6 continue                                                           1930
      call que(j,l,nt,jw,n,x,tc(la),sc)                                  1931
      do 7 i=1,n                                                         1932
      t(i)=t(i)+tc(lb+j)*sc(i)                                           1933
    7 continue                                                           1934
    8 continue                                                           1935
      go to 11                                                           1936
    9 la=la+lp(3,l1)*(5*lp(1,l1)+1)                                      1937
   10 continue                                                           1938
   11 return                                                             1939
      end                                                                1940
      subroutine orgpc (xm,xs,lp,lv,tc)                                  1941
      integer lp(3,*),lv(*)                                              1942
      real xm(*),xs(*),tc(*)                                             1943
      la=1                                                               1944
      l1=la                                                              1945
    1 if(lp(1,l1).eq.0) go to 3                                          1946
      l=lp(1,l1)                                                         1947
      nt=lp(3,l1)                                                        1948
      lb=la+5*l*nt-1                                                     1949
      do 2 j=1,nt                                                        1950
      call scpc(xm,xs,j,l,nt,lv(lp(2,l1)),tc(la),tc(lb+j))               1951
    2 continue                                                           1952
      la=lb+nt+1                                                         1953
      l1=l1+1                                                            1954
      go to 1                                                            1955
    3 return                                                             1956
      end                                                                1957
      subroutine pair (jv,n,x,nk,tb,cm,jl,kv,f,sc,js)                    1958
      integer jv(2),kv(2,jl),js(*)                                       1959
      real x(n,*),tb(5,nk),cm(*),f(n),sc(n)                              1960
      call fun(2,jv,n,x,nk,tb,cm,jl,kv,f,js)                             1961
      do 2 k=1,2                                                         1962
      call fun(1,jv(k),n,x(1,k),nk,tb,cm,jl,kv,sc,js)                    1963
      do 1 i=1,n                                                         1964
      f(i)=f(i)+sc(i)                                                    1965
    1 continue                                                           1966
    2 continue                                                           1967
      return                                                             1968
      end                                                                1969
      subroutine cpair (jv,n,x,nf,lp,lv,tc,f,sc)                         1970
      integer jv(2),lp(3,*),lv(*),jw(2)                                  1971
      real x(n,*),tc(*),f(n),sc(n,2)                                     1972
      call cfun(2,jv,n,x,nf,lp,lv,tc,f,sc,jw)                            1973
      do 2 k=1,2                                                         1974
      call cfun(1,jv(k),n,x(1,k),nf,lp,lv,tc,sc,sc(1,2),jw)              1975
      do 1 i=1,n                                                         1976
      f(i)=f(i)+sc(i,1)                                                  1977
    1 continue                                                           1978
    2 continue                                                           1979
      return                                                             1980
      end                                                                1981
      subroutine logitl (n,x,y,w,nk,il,az,tb,cm,sc,d)                    1982
      integer kp(5,*),kv(2,*),lp(3,*),lv(*)                              1983
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,*),tc(*),ss(n)           1984
      double precision d(nk,*),a,b,s,sw,yb                               1985
      data niter,wm,thr /30,0.0001,0.0001/                               1986
      do 2 i=1,n                                                         1987
      k=0                                                                1988
      do 1 m=1,nk                                                        1989
      if(tb(1,m).eq.0.0) go to 1                                         1990
      k=k+1                                                              1991
      sc(i,k)=phi(m,i,n,x,tb,cm)                                         1992
    1 continue                                                           1993
    2 continue                                                           1994
      if(k .ne. 0) go to 3                                               1995
      az=alog(az/(1.0-az))                                               1996
      return                                                             1997
    3 mk=k                                                               1998
      a=az                                                               1999
      jnt=1                                                              2000
      go to 19                                                           2001
      entry logitc (n,x,y,w,nk,il,cm,kp,kv,lp,lv,bz,tc,sc,ss,d)          2002
      ll=1                                                               2003
      la=ll                                                              2004
      l1=la                                                              2005
      lt=0                                                               2006
    4 if(kp(1,ll).lt.0) go to 17                                         2007
      do 5 i=1,n                                                         2008
      ss(i)=1.0                                                          2009
    5 continue                                                           2010
      if(kp(1,ll) .le. 0) go to 9                                        2011
      jl=kp(1,ll)                                                        2012
      do 8 il=1,jl                                                       2013
      k=kp(2,ll)+il-1                                                    2014
      jj=kv(1,k)                                                         2015
      j=iabs(jj)                                                         2016
      kk=kv(2,k)                                                         2017
      do 7 i=1,n                                                         2018
      if(ss(i).eq.0.0) go to 7                                           2019
      ic=x(i,j)+.1                                                       2020
      ss(i)=cm(ic+kk)                                                    2021
      if(jj .ge. 0) go to 7                                              2022
      if(ss(i) .ne. 0.0) go to 6                                         2023
      ss(i)=1.0                                                          2024
      go to 7                                                            2025
    6 ss(i)=0.0                                                          2026
    7 continue                                                           2027
    8 continue                                                           2028
      go to 10                                                           2029
    9 if(kp(3,ll) .gt. 0) go to 10                                       2030
      ll=ll+1                                                            2031
      go to 4                                                            2032
   10 if(kp(3,ll) .gt. 0) go to 12                                       2033
      lt=lt+1                                                            2034
      do 11 i=1,n                                                        2035
      sc(i,lt)=ss(i)                                                     2036
   11 continue                                                           2037
      go to 16                                                           2038
   12 kp3=kp(3,ll)                                                       2039
      do 15 m=1,kp3                                                      2040
      l=lp(1,l1)                                                         2041
      nt=lp(3,l1)                                                        2042
      do 14 jp=1,nt                                                      2043
      lt=lt+1                                                            2044
      do 13 i=1,n                                                        2045
      sc(i,lt)=ss(i)                                                     2046
   13 continue                                                           2047
      call que(jp,l,nt,lv(lp(2,l1)),n,x,tc(la),sc(1,lt))                 2048
   14 continue                                                           2049
      l1=l1+1                                                            2050
      la=la+nt*(5*l+1)                                                   2051
   15 continue                                                           2052
   16 ll=ll+1                                                            2053
      go to 4                                                            2054
   17 if(lt .ne. 0) go to 18                                             2055
      bz=alog(bz/(1.0-bz))                                               2056
      return                                                             2057
   18 mk=lt                                                              2058
      a=bz                                                               2059
      jnt=2                                                              2060
   19 mkp1=mk+1                                                          2061
      mkp2=mk+2                                                          2062
      mkp3=mk+3                                                          2063
      mkp4=mk+4                                                          2064
      iter=0                                                             2065
      if(jnt .ne. 1) go to 21                                            2066
      k=0                                                                2067
      do 20 m=1,nk                                                       2068
      if(tb(1,m).eq.0.0) go to 20                                        2069
      k=k+1                                                              2070
      d(k,mkp3)=tb(1,m)                                                  2071
   20 continue                                                           2072
      go to 27                                                           2073
   21 ll=1                                                               2074
      l1=ll                                                              2075
      la=0                                                               2076
      lt=la                                                              2077
   22 if(kp(1,ll).lt.0) go to 27                                         2078
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 23                  2079
      ll=ll+1                                                            2080
      go to 22                                                           2081
   23 if(kp(3,ll) .gt. 0) go to 24                                       2082
      lt=lt+1                                                            2083
      d(lt,mkp3)=tc(-kp(3,ll))                                           2084
      ll=ll+1                                                            2085
      go to 22                                                           2086
   24 kp3=kp(3,ll)                                                       2087
      do 26 m=1,kp3                                                      2088
      nt=lp(3,l1)                                                        2089
      la=la+5*lp(1,l1)*nt                                                2090
      do 25 i=1,nt                                                       2091
      lt=lt+1                                                            2092
      d(lt,mkp3)=tc(i+la)                                                2093
   25 continue                                                           2094
      la=la+nt                                                           2095
      l1=l1+1                                                            2096
   26 continue                                                           2097
      ll=ll+1                                                            2098
      go to 22                                                           2099
   27 iter=iter+1                                                        2100
      b=0.d0                                                             2101
      sw=b                                                               2102
      yb=sw                                                              2103
      do 29 i=1,n                                                        2104
      s=a                                                                2105
      do 28 m=1,mk                                                       2106
      s=s+d(m,mkp3)*sc(i,m)                                              2107
   28 continue                                                           2108
      sc(i,mkp3)=s                                                       2109
      pp=1.0/(1.0+exp(-sc(i,mkp3)))                                      2110
      ww=amax1(pp*(1.0-pp),wm)                                           2111
      sc(i,mkp3)=sc(i,mkp3)+(y(i)-pp)/ww                                 2112
      if(il.eq.2) ww=ww**2                                               2113
      ww=ww*w(i)                                                         2114
      sc(i,mkp2)=ww                                                      2115
      sw=sw+ww                                                           2116
      yb=yb+ww*sc(i,mkp3)                                                2117
      if(iter.gt.1) b=b+abs(pp-sc(i,mkp1))                               2118
      sc(i,mkp1)=pp                                                      2119
   29 continue                                                           2120
      if(iter.gt.niter.or.(iter.gt.1.and.b/n.lt.thr)) go to 37           2121
      yb=yb/sw                                                           2122
      do 36 m=1,mk                                                       2123
      b=0.d0                                                             2124
      do 30 i=1,n                                                        2125
      b=b+sc(i,mkp2)*sc(i,m)                                             2126
   30 continue                                                           2127
      b=b/sw                                                             2128
      mm1=m-1                                                            2129
      l=1                                                                2130
      go to 32                                                           2131
   31 l=l+1                                                              2132
   32 if((l).gt.(mm1)) go to 34                                          2133
      s=0.d0                                                             2134
      do 33 i=1,n                                                        2135
      s=s+sc(i,mkp2)*(sc(i,m)-b)*sc(i,l)                                 2136
   33 continue                                                           2137
      d(l,m)=s                                                           2138
      go to 31                                                           2139
   34 a=0.d0                                                             2140
      s=a                                                                2141
      do 35 i=1,n                                                        2142
      ww=sc(i,mkp2)                                                      2143
      pp=sc(i,m)-b                                                       2144
      s=s+ww*pp**2                                                       2145
      a=a+ww*pp*sc(i,mkp3)                                               2146
   35 continue                                                           2147
      d(m,m)=s                                                           2148
      d(m,mkp1)=a                                                        2149
      d(m,mkp2)=b                                                        2150
   36 continue                                                           2151
      call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)                2152
      go to 27                                                           2153
   37 if(jnt .ne. 1) go to 39                                            2154
      az=a                                                               2155
      k=0                                                                2156
      do 38 m=1,nk                                                       2157
      if(tb(1,m).eq.0.0) go to 38                                        2158
      k=k+1                                                              2159
      tb(1,m)=d(k,mkp3)                                                  2160
   38 continue                                                           2161
      go to 45                                                           2162
   39 bz=a                                                               2163
      ll=1                                                               2164
      l1=ll                                                              2165
      la=0                                                               2166
      lt=la                                                              2167
   40 if(kp(1,ll).lt.0) go to 45                                         2168
      if(kp(1,ll) .ne. 0 .or. kp(3,ll) .gt. 0) go to 41                  2169
      ll=ll+1                                                            2170
      go to 40                                                           2171
   41 if(kp(3,ll) .gt. 0) go to 42                                       2172
      lt=lt+1                                                            2173
      tc(-kp(3,ll))=d(lt,mkp3)                                           2174
      ll=ll+1                                                            2175
      go to 40                                                           2176
   42 kp3=kp(3,ll)                                                       2177
      do 44 m=1,kp3                                                      2178
      nt=lp(3,l1)                                                        2179
      la=la+5*lp(1,l1)*nt                                                2180
      do 43 i=1,nt                                                       2181
      lt=lt+1                                                            2182
      tc(i+la)=d(lt,mkp3)                                                2183
   43 continue                                                           2184
      la=la+nt                                                           2185
      l1=l1+1                                                            2186
   44 continue                                                           2187
      ll=ll+1                                                            2188
      go to 40                                                           2189
   45 return                                                             2190
      end                                                                2191
      subroutine varimp (n,p,x,y,w,nk,il,it,az,tb,cm,vip,sc,d)           2192
      integer p                                                          2193
      real x(n,p),y(n),w(n),tb(5,nk),cm(*),vip(p),sc(n,*)                2194
      double precision d(nk,*),sw,yb,yv,wn                               2195
      sw=0.d0                                                            2196
      yb=sw                                                              2197
      yv=yb                                                              2198
      wn=yv                                                              2199
      do 1 i=1,n                                                         2200
      sw=sw+w(i)                                                         2201
      wn=wn+w(i)**2                                                      2202
      yb=yb+w(i)*y(i)                                                    2203
    1 continue                                                           2204
      yb=yb/sw                                                           2205
      do 2 i=1,n                                                         2206
      yv=yv+w(i)*(y(i)-yb)**2                                            2207
    2 continue                                                           2208
      yv=yv/sw                                                           2209
      wn=sw**2/wn                                                        2210
      ip=nk+4                                                            2211
      call varz(0,nk,tb,sc(1,ip),cst,nd)                                 2212
      if(cst .ne. 1.0) go to 3                                           2213
      g0=0.0                                                             2214
      if(il.gt.0) g0=yv                                                  2215
      go to 4                                                            2216
    3 a0=az                                                              2217
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g0,sc,d)                2218
    4 cst=1.d0/(1.d0-cst/wn)**2                                          2219
      if(il .ne. 0) go to 5                                              2220
      g0=(g0+yv)*cst                                                     2221
      go to 6                                                            2222
    5 g0=g0*cst                                                          2223
    6 do 12 j=1,p                                                        2224
      call varz(j,nk,tb,sc(1,ip),cst,nd)                                 2225
      if(nd .ne. 0) go to 7                                              2226
      vip(j)=g0                                                          2227
      go to 12                                                           2228
    7 if(cst .ne. 1.0) go to 8                                           2229
      g=0.0                                                              2230
      if(il.gt.0) g=yv                                                   2231
      go to 9                                                            2232
    8 a0=az                                                              2233
      call vp(n,x,y,w,nk,il,yb,sw,a0,sc(1,ip),cm,g,sc,d)                 2234
    9 cst=1.d0/(1.d0-cst/wn)**2                                          2235
      if(il .ne. 0) go to 10                                             2236
      g=(g+yv)*cst                                                       2237
      go to 11                                                           2238
   10 g=g*cst                                                            2239
   11 vip(j)=g                                                           2240
   12 continue                                                           2241
      if(it .le. 0) go to 13                                             2242
c     write(it,17)                                                       2243
      call numprt(it,p,vip)                                              2244
   13 a0=0.0                                                             2245
      do 14 j=1,p                                                        2246
      vip(j)=sqrt(amax1(0.0,vip(j)-g0))                                  2247
      a0=amax1(a0,vip(j))                                                2248
   14 continue                                                           2249
      if(a0.le.0.0) return                                               2250
      do 15 j=1,p                                                        2251
      vip(j)=100.0*vip(j)/a0                                             2252
   15 continue                                                           2253
      if(it .le. 0) go to 16                                             2254
c     write(it,18)                                                       2255
      call numprt(it,p,vip)                                              2256
   16 return                                                             2257
   17 format(/,' -gcv removing each variable:')                          2258
   18 format(/,' relative variable importance:')                         2259
      end                                                                2260
      subroutine numprt(it,n,a)                                          2261
      real a(*)                                                          2262
      i2=0                                                               2263
    1 if(i2.ge.n) go to 2                                                2264
      i1=i2+1                                                            2265
      i2=i2+6                                                            2266
      if(i2.gt.n) i2=n                                                   2267
c     write(it,'(/,'' '',6(''    '',i4,''    ''))') (i,i=i1,i2)          2268
c     write(it,'('' '',6g12.4)') (a(i),i=i1,i2)                          2269
      go to 1                                                            2270
    2 return                                                             2271
      end                                                                2272
      subroutine varz(j,nk,tb,ub,cst,nd)                                 2273
      real tb(5,nk),ub(5,nk)                                             2274
      do 2 m=1,nk                                                        2275
      do 1 k=1,5                                                         2276
      ub(k,m)=tb(k,m)                                                    2277
    1 continue                                                           2278
    2 continue                                                           2279
      nd=0                                                               2280
      if(j .le. 0) go to 4                                               2281
      do 3 m=1,nk                                                        2282
      if(ub(1,m).eq.0.0) go to 3                                         2283
      if(jf(m,j,ub) .eq. 0) go to 3                                      2284
      ub(1,m)=0.0                                                        2285
      nd=nd+1                                                            2286
    3 continue                                                           2287
    4 cst=1.0                                                            2288
      do 5 m=1,nk                                                        2289
      if(ub(1,m).ne.0.0) cst=cst+ub(5,m)                                 2290
    5 continue                                                           2291
      return                                                             2292
      end                                                                2293
      subroutine vp (n,x,y,w,nk,il,yb,sw,az,tb,cm,gof,sc,d)              2294
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)                      2295
      double precision d(nk,*),s,t,yb,sw                                 2296
      if(il .ne. 0) go to 1                                              2297
      call lstsqr(n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)                       2298
      return                                                             2299
    1 call logitl(n,x,y,w,nk,il,az,tb,cm,sc,d)                           2300
      t=0.d0                                                             2301
      do 3 i=1,n                                                         2302
      s=az                                                               2303
      k=0                                                                2304
      do 2 m=1,nk                                                        2305
      if(tb(1,m).eq.0.0) go to 2                                         2306
      k=k+1                                                              2307
      s=s+tb(1,m)*sc(i,k)                                                2308
    2 continue                                                           2309
      a=s                                                                2310
      pp=1.0/(1.0+exp(-a))                                               2311
      t=t+w(i)*(y(i)-pp)**2                                              2312
    3 continue                                                           2313
      gof=t/sw                                                           2314
      return                                                             2315
      end                                                                2316
      subroutine lstsqr (n,x,y,w,nk,yb,sw,tb,cm,gof,sc,d)                2317
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n,nk)                      2318
      double precision d(nk,*),a,b,s,yb,sw                               2319
      do 2 i=1,n                                                         2320
      k=0                                                                2321
      do 1 m=1,nk                                                        2322
      if(tb(1,m).eq.0.0) go to 1                                         2323
      k=k+1                                                              2324
      sc(i,k)=phi(m,i,n,x,tb,cm)                                         2325
    1 continue                                                           2326
    2 continue                                                           2327
      mk=k                                                               2328
      mkp1=mk+1                                                          2329
      mkp2=mk+2                                                          2330
      mkp3=mk+3                                                          2331
      mkp4=mk+4                                                          2332
      do 9 m=1,mk                                                        2333
      b=0.d0                                                             2334
      do 3 i=1,n                                                         2335
      b=b+w(i)*sc(i,m)                                                   2336
    3 continue                                                           2337
      b=b/sw                                                             2338
      mm1=m-1                                                            2339
      l=1                                                                2340
      go to 5                                                            2341
    4 l=l+1                                                              2342
    5 if((l).gt.(mm1)) go to 7                                           2343
      s=0.d0                                                             2344
      do 6 i=1,n                                                         2345
      s=s+w(i)*(sc(i,m)-b)*sc(i,l)                                       2346
    6 continue                                                           2347
      d(l,m)=s                                                           2348
      go to 4                                                            2349
    7 a=0.d0                                                             2350
      s=a                                                                2351
      do 8 i=1,n                                                         2352
      ww=w(i)                                                            2353
      pp=sc(i,m)-b                                                       2354
      s=s+ww*pp**2                                                       2355
      a=a+ww*pp*y(i)                                                     2356
    8 continue                                                           2357
      d(m,m)=s                                                           2358
      d(m,mkp1)=a                                                        2359
      d(m,mkp2)=b                                                        2360
    9 continue                                                           2361
      call lsf(nk,mk,mkp1,yb,d,d(1,mkp3),a,s,d(1,mkp4),1)                2362
      gof=s/sw                                                           2363
      return                                                             2364
      end                                                                2365
      function efp(l,jv,nk,tb)                                           2366
      integer jv(l)                                                      2367
      real tb(5,nk)                                                      2368
      efp=0.0                                                            2369
      do 3 m=1,nk                                                        2370
      if(tb(1,m).eq.0.0) go to 3                                         2371
      if(nord(m,tb).ne.l) go to 3                                        2372
      k=0                                                                2373
      do 1 j=1,l                                                         2374
      if(jf(m,jv(j),tb).eq.1) go to 1                                    2375
      k=1                                                                2376
      go to 2                                                            2377
    1 continue                                                           2378
    2 if(k.eq.1) go to 3                                                 2379
      efp=efp+tb(5,m)                                                    2380
    3 continue                                                           2381
      return                                                             2382
      end                                                                2383
      function elg(jv,l,lx,tb,cm)                                        2384
      real tb(5,*),cm(*)                                                 2385
      integer lx(*)                                                      2386
      logical elg                                                        2387
      data ic /0/                                                        2388
      elg=.false.                                                        2389
      kx=iabs(lx(jv))                                                    2390
      if(kx.eq.0) return                                                 2391
      if(l .ne. 0) go to 1                                               2392
      elg=.true.                                                         2393
      return                                                             2394
    1 if((kx .ne. 2) .and. (kx .ne. 3)) go to 2                          2395
      if(nnord(l,tb).gt.0) return                                        2396
    2 ip=l                                                               2397
    3 if(ip.le.0) go to 4                                                2398
      jl=abs(tb(2,ip))+.1                                                2399
      ip=tb(4,ip)+.1                                                     2400
      go to 3                                                            2401
    4 k=iabs(lx(jl))                                                     2402
      call isnstr(jl,jb)                                                 2403
      if((k.eq.2.or.k.eq.3).and.jb.eq.0) return                          2404
      if(ic .ne. 1) go to 5                                              2405
      if(lx(jv).lt.0.and.nordc(1,l,tb,cm).gt.0) return                   2406
      if(lx(jv).gt.0.and.nordc(2,l,tb,cm).gt.0) return                   2407
      go to 6                                                            2408
    5 if(ic .ne. 2) go to 6                                              2409
      if(lx(jv).gt.0.and.nordc(1,l,tb,cm).ge.2) return                   2410
    6 ip=l                                                               2411
    7 if(ip.le.0) go to 8                                                2412
      jl=abs(tb(2,ip))+.1                                                2413
      call intalw(jv,jl,k)                                               2414
      if(k.eq.0) return                                                  2415
      ip=tb(4,ip)+.1                                                     2416
      go to 7                                                            2417
    8 elg=.true.                                                         2418
      return                                                             2419
      entry stelg(i1)                                                    2420
      ic=i1                                                              2421
      return                                                             2422
      end                                                                2423
      function phi(m,i,n,x,tb,cm)                                        2424
      real tb(5,*),x(n,*),cm(*)                                          2425
      phi=1.0                                                            2426
      ip=m                                                               2427
    1 if(ip.le.0) go to 7                                                2428
      t=tb(2,ip)                                                         2429
      j=abs(t)+.1                                                        2430
      if(cm(2*j) .le. 0.0) go to 4                                       2431
      u=cm(int(x(i,j)+.1)+int(tb(3,ip)+.1))                              2432
      if(t .ge. 0.0) go to 5                                             2433
      if(u .ne. 0.0) go to 2                                             2434
      u=1.0                                                              2435
      go to 5                                                            2436
    2 u=0.0                                                              2437
      go to 5                                                            2438
    4 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,ip)))                         2439
    5 if(u .gt. 0.0) go to 6                                             2440
      phi=0.0                                                            2441
      return                                                             2442
    6 phi=phi*u                                                          2443
      ip=tb(4,ip)+.1                                                     2444
      go to 1                                                            2445
    7 return                                                             2446
      end                                                                2447
      function nord(m,tb)                                                2448
      real tb(5,*)                                                       2449
      ip=m                                                               2450
      nord=0                                                             2451
    1 if(ip.le.0) go to 2                                                2452
      nord=nord+1                                                        2453
      ip=tb(4,ip)+.1                                                     2454
      go to 1                                                            2455
    2 return                                                             2456
      end                                                                2457
      function jf(m,j,tb)                                                2458
      real tb(5,*)                                                       2459
      ip=m                                                               2460
      jf=0                                                               2461
    1 if(ip.le.0) go to 2                                                2462
      jp=abs(tb(2,ip))+.1                                                2463
      if(jp.eq.j) jf=1                                                   2464
      ip=tb(4,ip)+.1                                                     2465
      go to 1                                                            2466
    2 return                                                             2467
      end                                                                2468
      subroutine lsf(nk,m,mkp1,yb,d,a,a0,gf,dp,k1)                       2469
      double precision a(m),d(nk,*),dp(nk,*),eps,yb,a0,gf,s,t,sl,dps     2470
      data big,eps /9.9e30,1.d-05/                                       2471
      mkp2=mkp1+1                                                        2472
      gf=big                                                             2473
      if(d(m,m).le.0.d0) return                                          2474
      sl=1.d0+eps                                                        2475
      do 2 k=k1,m                                                        2476
      do 1 i=1,k                                                         2477
      dp(i,k)=d(i,k)                                                     2478
    1 continue                                                           2479
      dp(k,k)=dp(k,k)*sl                                                 2480
    2 continue                                                           2481
      do 3 k=1,m                                                         2482
      a(k)=d(k,mkp1)                                                     2483
    3 continue                                                           2484
      info=k1                                                            2485
      call spofa(dp,nk,m,info)                                           2486
      if(info.ne.0) return                                               2487
      call sposl(dp,nk,m,a)                                              2488
      s=yb                                                               2489
      t=0.d0                                                             2490
      do 4 i=1,m                                                         2491
      s=s-a(i)*d(i,mkp2)                                                 2492
      t=t-a(i)*(d(i,mkp1)+eps*d(i,i)*a(i))                               2493
    4 continue                                                           2494
      a0=s                                                               2495
      gf=t                                                               2496
      return                                                             2497
      entry seteps(dps)                                                  2498
      eps=dps                                                            2499
      return                                                             2500
      end                                                                2501
      subroutine coll(nk,tb,lp,lv,jv)                                    2502
      integer lp(3,*),lv(*),jv(*)                                        2503
      real tb(5,nk)                                                      2504
      mo=0                                                               2505
      do 1 m=1,nk                                                        2506
      if(tb(1,m).ne.0.0) mo=max0(mo,nord(m,tb))                          2507
    1 continue                                                           2508
      if(mo .ne. 0) go to 2                                              2509
      lp(1,1)=0                                                          2510
      return                                                             2511
    2 l1=1                                                               2512
      l2=l1                                                              2513
      do 11 mt=1,mo                                                      2514
      l10=l1                                                             2515
      do 10 m=1,nk                                                       2516
      if(tb(1,m).eq.0.0.or.nord(m,tb).ne.mt) go to 10                    2517
      call jfv(m,tb,jv)                                                  2518
      jg=0                                                               2519
      l1m1=l1-1                                                          2520
      i=l10                                                              2521
      go to 4                                                            2522
    3 i=i+1                                                              2523
    4 if((i).gt.(l1m1)) go to 8                                          2524
      k=lp(2,i)-1                                                        2525
      ig=0                                                               2526
      do 5 j=1,mt                                                        2527
      if(jv(j).eq.lv(k+j)) go to 5                                       2528
      ig=1                                                               2529
      go to 6                                                            2530
    5 continue                                                           2531
    6 if(ig .ne. 0) go to 3                                              2532
      jg=1                                                               2533
      lp(3,i)=lp(3,i)+1                                                  2534
    8 if(jg .ne. 0) go to 10                                             2535
      lp(1,l1)=mt                                                        2536
      lp(2,l1)=l2                                                        2537
      lp(3,l1)=1                                                         2538
      k=l2-1                                                             2539
      do 9 i=1,mt                                                        2540
      lv(i+k)=jv(i)                                                      2541
    9 continue                                                           2542
      l1=l1+1                                                            2543
      l2=l2+mt                                                           2544
   10 continue                                                           2545
   11 continue                                                           2546
      lp(1,l1)=0                                                         2547
      return                                                             2548
      end                                                                2549
      subroutine jfv(m,tb,jv)                                            2550
      integer jv(*)                                                      2551
      real tb(5,*)                                                       2552
      ip=m                                                               2553
      j=0                                                                2554
    1 if(ip.le.0) go to 2                                                2555
      j=j+1                                                              2556
      jv(j)=abs(tb(2,ip))+.1                                             2557
      ip=tb(4,ip)+.1                                                     2558
      go to 1                                                            2559
    2 if(j.eq.1) return                                                  2560
      j=j-1                                                              2561
    3 l=0                                                                2562
      do 4 i=1,j                                                         2563
      if(jv(i) .le. jv(i+1)) go to 4                                     2564
      k=jv(i)                                                            2565
      jv(i)=jv(i+1)                                                      2566
      jv(i+1)=k                                                          2567
      l=1                                                                2568
    4 continue                                                           2569
      if(l.eq.0) go to 5                                                 2570
      go to 3                                                            2571
    5 return                                                             2572
      end                                                                2573
      subroutine side(l,nt,jv,xe,x)                                      2574
      integer jv(l)                                                      2575
      real xe(2,*),x(nt,*)                                               2576
      l2=l+l                                                             2577
      l3=l2+l                                                            2578
      l4=l3+l                                                            2579
      do 7 k=1,l                                                         2580
      xl=xe(1,jv(k))                                                     2581
      xr=xe(2,jv(k))                                                     2582
      do 6 j=1,nt                                                        2583
      z=x(j,k)                                                           2584
      if(z .gt. xl) go to 1                                              2585
      x(j,k+l)=xl                                                        2586
      x(j,k+l2)=x(j,k+l)                                                 2587
      x(j,k+l3)=0.0                                                      2588
      x(j,k+l4)=x(j,k+l3)                                                2589
      go to 6                                                            2590
    1 dl=z-xl                                                            2591
      dr=xr-z                                                            2592
      x1=xl                                                              2593
      x2=xr                                                              2594
      do 3 m=1,nt                                                        2595
      a=x(m,k)                                                           2596
      if(a.eq.z) go to 3                                                 2597
      dx=a-z                                                             2598
      if(dx .ge. 0.0 .or. -dx .ge. dl) go to 2                           2599
      dl=-dx                                                             2600
      x1=a                                                               2601
    2 if(dx .le. 0.0 .or. dx .ge. dr) go to 3                            2602
      dr=dx                                                              2603
      x2=a                                                               2604
    3 continue                                                           2605
      x1=0.5*(x1+z)                                                      2606
      x2=0.5*(x2+z)                                                      2607
      if(x(j,k+l) .le. 0.0) go to 4                                      2608
      x(j,k+l)=x1                                                        2609
      x(j,k+l2)=x2                                                       2610
      go to 5                                                            2611
    4 x(j,k+l)=x2                                                        2612
      x(j,k+l2)=x1                                                       2613
    5 call pr(x(j,k+l),x(j,k),x(j,k+l2),x(j,k+l3),x(j,k+l4))             2614
    6 continue                                                           2615
    7 continue                                                           2616
      return                                                             2617
      end                                                                2618
      subroutine que(jp,l,nt,jv,n,x,tc,t)                                2619
      integer jv(l)                                                      2620
      real x(n,*),tc(nt,*),t(n)                                          2621
      l2=l+l                                                             2622
      l3=l2+l                                                            2623
      l4=l3+l                                                            2624
      do 3 i=1,n                                                         2625
      if(t(i).eq.0.0) go to 3                                            2626
      q=1.0                                                              2627
      do 1 k=1,l                                                         2628
      j=jv(k)                                                            2629
      q=q*cue(x(i,j),tc(jp,k+l),tc(jp,k),tc(jp,k+l2),tc(jp,k+l3),tc(jp,k 2630
     1+l4))                                                              2631
      if(q.eq.0.0) go to 2                                               2632
    1 continue                                                           2633
    2 t(i)=q                                                             2634
    3 continue                                                           2635
      return                                                             2636
      end                                                                2637
      subroutine scpc(xm,xs,jp,l,nt,jv,tc,b)                             2638
      integer jv(l)                                                      2639
      real xm(*),xs(*),tc(nt,*)                                          2640
      double precision g,h,q                                             2641
      l2=l+l                                                             2642
      l3=l2+l                                                            2643
      l4=l3+l                                                            2644
      q=1.d0                                                             2645
      do 1 k=1,l                                                         2646
      j=jv(k)                                                            2647
      g=xm(j)                                                            2648
      h=xs(j)                                                            2649
      q=q*h                                                              2650
      tc(jp,k+l)=g+h*tc(jp,k+l)                                          2651
      tc(jp,k)=g+h*tc(jp,k)                                              2652
      tc(jp,k+l2)=g+h*tc(jp,k+l2)                                        2653
      tc(jp,k+l3)=tc(jp,k+l3)/h                                          2654
      tc(jp,k+l4)=tc(jp,k+l4)/h**2                                       2655
    1 continue                                                           2656
      b=b/q                                                              2657
      return                                                             2658
      end                                                                2659
      subroutine update(il,n,m,kr,x,y,w,sw,yb,tb,cm,sc,bl,d,dy,db)       2660
      real x(n,*),y(n),w(n),tb(5,*),cm(*),sc(n,*),bl(n)                  2661
      double precision d(n,*),dy(*),db(*),b,s,yb,sw,dv,eps,v,q           2662
      data eps /1.d-4/                                                   2663
      kp=kr+1                                                            2664
      b=0.d0                                                             2665
      t=tb(2,m)                                                          2666
      j=abs(t)+.1                                                        2667
      if(il .ne. 1) go to 3                                              2668
      tk=tb(3,m)                                                         2669
      do 2 i=1,n                                                         2670
      h=bl(i)                                                            2671
      if(h .gt. 0.0) go to 1                                             2672
      sc(i,m)=0.0                                                        2673
      go to 2                                                            2674
    1 sc(i,m)=h*(x(i,j)-tk)                                              2675
      b=b+w(i)*sc(i,m)                                                   2676
    2 continue                                                           2677
      go to 17                                                           2678
    3 if(cm(2*j) .le. 0.0) go to 12                                      2679
      k=tb(3,m)+.1                                                       2680
      nw=0                                                               2681
      n0=nw                                                              2682
      do 11 i=1,n                                                        2683
      h=bl(i)                                                            2684
      if(h .gt. 0.0) go to 4                                             2685
      sc(i,m)=0.0                                                        2686
      go to 11                                                           2687
    4 u=cm(int(x(i,j)+.1)+k)                                             2688
      if(w(i) .le. 0.0) go to 5                                          2689
      nw=nw+1                                                            2690
      if(u.eq.0.0) n0=n0+1                                               2691
    5 if(t .ge. 0.0) go to 8                                             2692
      if(u .ne. 0.0) go to 6                                             2693
      sc(i,m)=h                                                          2694
      go to 10                                                           2695
    6 sc(i,m)=0.0                                                        2696
      go to 10                                                           2697
    8 if(u .gt. 0.0) go to 9                                             2698
      sc(i,m)=0.0                                                        2699
      go to 10                                                           2700
    9 sc(i,m)=h                                                          2701
   10 b=b+w(i)*sc(i,m)                                                   2702
   11 continue                                                           2703
      if(n0.eq.0.or.n0.eq.nw) return                                     2704
      go to 17                                                           2705
   12 tk=tb(3,m)                                                         2706
      sg=sign(1.0,t)                                                     2707
      do 16 i=1,n                                                        2708
      h=bl(i)                                                            2709
      if(h .gt. 0.0) go to 13                                            2710
      sc(i,m)=0.0                                                        2711
      go to 16                                                           2712
   13 u=amax1(0.0,sg*(x(i,j)-tk))                                        2713
      if(u .gt. 0.0) go to 14                                            2714
      sc(i,m)=0.0                                                        2715
      go to 15                                                           2716
   14 sc(i,m)=h*u                                                        2717
   15 b=b+w(i)*sc(i,m)                                                   2718
   16 continue                                                           2719
   17 b=b/sw                                                             2720
      s=0.d0                                                             2721
      v=s                                                                2722
      do 18 j=1,kr                                                       2723
      db(j)=0.d0                                                         2724
   18 continue                                                           2725
      do 21 i=1,n                                                        2726
      d(i,kp)=sc(i,m)-b                                                  2727
      if(sc(i,m).le.0.0) go to 21                                        2728
      q=w(i)*sc(i,m)                                                     2729
      s=s+q*(sc(i,m)-b)                                                  2730
      v=v+q*(y(i)-yb)                                                    2731
      j=1                                                                2732
      go to 20                                                           2733
   19 j=j+1                                                              2734
   20 if((j).gt.(kr)) go to 21                                           2735
      db(j)=db(j)+q*d(i,j)                                               2736
      go to 19                                                           2737
   21 continue                                                           2738
      if(s.le.0.d0) return                                               2739
      dv=s                                                               2740
      j=1                                                                2741
      go to 23                                                           2742
   22 j=j+1                                                              2743
   23 if((j).gt.(kr)) go to 24                                           2744
      s=s-db(j)**2                                                       2745
      go to 22                                                           2746
   24 if(s.lt.eps*dv) return                                             2747
      j=1                                                                2748
      go to 26                                                           2749
   25 j=j+1                                                              2750
   26 if((j).gt.(kr)) go to 28                                           2751
      do 27 i=1,n                                                        2752
      d(i,kp)=d(i,kp)-db(j)*d(i,j)                                       2753
   27 continue                                                           2754
      go to 25                                                           2755
   28 s=1.d0/dsqrt(s)                                                    2756
      j=1                                                                2757
      go to 30                                                           2758
   29 j=j+1                                                              2759
   30 if((j).gt.(kr)) go to 31                                           2760
      v=v-db(j)*dy(j)                                                    2761
      go to 29                                                           2762
   31 dy(kp)=v*s                                                         2763
      do 32 i=1,n                                                        2764
      d(i,kp)=d(i,kp)*s                                                  2765
   32 continue                                                           2766
      kr=kp                                                              2767
      return                                                             2768
      end                                                                2769
      subroutine pr(um,u,up,p,r)                                         2770
      s=1.0                                                              2771
      if(um.gt.up) s=-1.0                                                2772
      p=s*(2.0*up+um-3.0*u)/(up-um)**2                                   2773
      r=s*(2.0*u-up-um)/(up-um)**3                                       2774
      return                                                             2775
      end                                                                2776
      function cue(x,um,u,up,p,r)                                        2777
      s=1.0                                                              2778
      if(um.gt.up) s=-1.0                                                2779
      y=s*x                                                              2780
      if(y .gt. s*um) go to 1                                            2781
      cue=0.0                                                            2782
      return                                                             2783
    1 if(y .lt. s*up) go to 2                                            2784
      cue=y-s*u                                                          2785
      return                                                             2786
    2 cue=p*(x-um)**2+r*(x-um)**3                                        2787
      return                                                             2788
      end                                                                2789
      function varf(nk,d,a,sw,k1,k2)                                     2790
      double precision d(nk,*),a(nk),sw,s,t,u                            2791
      s=0.d0                                                             2792
      do 4 i=k1,k2                                                       2793
      t=0.d0                                                             2794
      do 3 j=k1,k2                                                       2795
      if(j .gt. i) go to 1                                               2796
      u=d(j,i)                                                           2797
      go to 2                                                            2798
    1 u=d(i,j)                                                           2799
    2 t=t+a(j)*u                                                         2800
    3 continue                                                           2801
      s=s+a(i)*t                                                         2802
    4 continue                                                           2803
      varf=s/sw                                                          2804
      return                                                             2805
      end                                                                2806
      subroutine exch(nk,m,k,d,a,b)                                      2807
      real a(m),b(m)                                                     2808
      double precision d(nk,m),t                                         2809
      l=k+1                                                              2810
      km1=k-1                                                            2811
      kp2=k+2                                                            2812
      r=a(k)                                                             2813
      a(k)=a(l)                                                          2814
      a(l)=r                                                             2815
      r=b(k)                                                             2816
      b(k)=b(l)                                                          2817
      b(l)=r                                                             2818
      do 1 j=1,2                                                         2819
      i=nk+j                                                             2820
      t=d(k,i)                                                           2821
      d(k,i)=d(l,i)                                                      2822
      d(l,i)=t                                                           2823
    1 continue                                                           2824
      t=d(k,k)                                                           2825
      d(k,k)=d(l,l)                                                      2826
      d(l,l)=t                                                           2827
      j=1                                                                2828
      go to 3                                                            2829
    2 j=j+1                                                              2830
    3 if((j).gt.(km1)) go to 4                                           2831
      t=d(j,k)                                                           2832
      d(j,k)=d(j,l)                                                      2833
      d(j,l)=t                                                           2834
      go to 2                                                            2835
    4 j=kp2                                                              2836
      go to 6                                                            2837
    5 j=j+1                                                              2838
    6 if((j).gt.(m)) go to 7                                             2839
      t=d(k,j)                                                           2840
      d(k,j)=d(l,j)                                                      2841
      d(l,j)=t                                                           2842
      go to 5                                                            2843
    7 return                                                             2844
      end                                                                2845
      function jft(m,j,tb)                                               2846
      real tb(5,*)                                                       2847
      k=1                                                                2848
      go to 2                                                            2849
    1 k=k+1                                                              2850
    2 if((k).gt.(m)) go to 4                                             2851
      if(int(abs(tb(2,k))+.1) .ne. j) go to 1                            2852
      jft=1                                                              2853
      return                                                             2854
    4 jft=0                                                              2855
      return                                                             2856
      end                                                                2857
      subroutine coefpr (it,nk,az,tb,cm,xs)                              2858
      real tb(5,*),cm(*),xs(*),a(6)                                      2859
      i2=0                                                               2860
    1 if(i2.ge.nk) go to 4                                               2861
      if(i2 .ne. 0) go to 2                                              2862
      i1=0                                                               2863
      i2=min0(5,nk)                                                      2864
      l2=i2+1                                                            2865
      a(1)=az                                                            2866
      call org(1,i2,tb,cm,xs,a(2))                                       2867
      go to 3                                                            2868
    2 i1=i2+1                                                            2869
      i2=i2+6                                                            2870
      if(i2.gt.nk) i2=nk                                                 2871
      l2=i2-i1+1                                                         2872
      call org(i1,i2,tb,cm,xs,a)                                         2873
    3 continue
c     write(it,'(/,'' bsfn:'',6(''    '',i4,''    ''))') (i,i=i1,i2)     2874
c     write(it,'('' coef:'',6g12.4)') (a(i),i=1,l2)                      2875
      go to 1                                                            2876
    4 return                                                             2877
      end                                                                2878
      subroutine org(m1,m2,tb,cm,xs,a)                                   2879
      real xs(*),tb(5,*),cm(*),a(*)                                      2880
      k=0                                                                2881
      do 4 m=m1,m2                                                       2882
      k=k+1                                                              2883
      if(tb(1,m) .ne. 0.0) go to 1                                       2884
      a(k)=0.0                                                           2885
      go to 4                                                            2886
    1 s=1.0                                                              2887
      ip=m                                                               2888
    2 if(ip.le.0) go to 3                                                2889
      j=abs(tb(2,ip))+.1                                                 2890
      if(cm(2*j).eq.0.0) s=s*xs(j)                                       2891
      ip=tb(4,ip)+.1                                                     2892
      go to 2                                                            2893
    3 a(k)=tb(1,m)/s                                                     2894
    4 continue                                                           2895
      return                                                             2896
      end                                                                2897
      subroutine hulset (n,x,big,nh,xh,y)                                2898
      real x(n,*),y(n),xh(3,nh)                                          2899
      do 5 j=1,n                                                         2900
      k=0                                                                2901
      x1=x(j,1)                                                          2902
      x2=x(j,2)                                                          2903
      do 3 i=1,nh                                                        2904
      a=xh(1,i)                                                          2905
      b=xh(2,i)                                                          2906
      sg=xh(3,i)                                                         2907
      if(a .lt. big) go to 1                                             2908
      s=x1-b                                                             2909
      go to 2                                                            2910
    1 s=x2-a*x1-b                                                        2911
    2 if(s*sg .ge. 0.0) go to 3                                          2912
      k=1                                                                2913
      go to 4                                                            2914
    3 continue                                                           2915
    4 if(k.eq.1) y(j)=big                                                2916
    5 continue                                                           2917
      return                                                             2918
      end                                                                2919
      subroutine cvxhul (n,x1,x2,big,nh,xh)                              2920
      real x1(n),x2(n),xh(3,*)                                           2921
      data eps /1.e-3/                                                   2922
      x0=big                                                             2923
      y0=x0                                                              2924
      xm=-big                                                            2925
      ym=xm                                                              2926
      xb=0.0                                                             2927
      yb=xb                                                              2928
      do 1 i=1,n                                                         2929
      xb=xb+x1(i)                                                        2930
      yb=yb+x2(i)                                                        2931
      x0=amin1(x0,x1(i))                                                 2932
      y0=amin1(y0,x2(i))                                                 2933
      xm=amax1(xm,x1(i))                                                 2934
      ym=amax1(ym,x2(i))                                                 2935
    1 continue                                                           2936
      x0=x0-eps*(xm-x0)                                                  2937
      y0=y0-eps*(ym-y0)                                                  2938
      xb=xb/n                                                            2939
      yb=yb/n                                                            2940
      nh=0                                                               2941
      a0=0.0                                                             2942
      lq=1                                                               2943
    2 am=big                                                             2944
      kq=4                                                               2945
      k=0                                                                2946
      do 15 i=1,n                                                        2947
      x=x1(i)                                                            2948
      y=x2(i)                                                            2949
      if(x .ne. x0) go to 5                                              2950
      if(y.eq.y0) go to 15                                               2951
      if(y .le. y0) go to 3                                              2952
      iq=2                                                               2953
      go to 10                                                           2954
    3 iq=4                                                               2955
      go to 10                                                           2956
    5 if(x .le. x0) go to 8                                              2957
      if(y .lt. y0) go to 6                                              2958
      iq=1                                                               2959
      go to 10                                                           2960
    6 iq=4                                                               2961
      go to 10                                                           2962
    8 if(y .le. y0) go to 9                                              2963
      iq=2                                                               2964
      go to 10                                                           2965
    9 iq=3                                                               2966
   10 if(iq.gt.kq) go to 15                                              2967
      if(iq.lt.lq) go to 15                                              2968
      if((iq .ne. 1) .and. (iq .ne. 3)) go to 11                         2969
      a=abs((y-y0)/(x-x0))                                               2970
      go to 12                                                           2971
   11 a=abs((x-x0)/(y-y0))                                               2972
   12 if(iq.eq.lq.and.a.lt.a0) go to 15                                  2973
      if(iq .ge. kq) go to 13                                            2974
      kq=iq                                                              2975
      am=a                                                               2976
      k=i                                                                2977
      go to 15                                                           2978
   13 if(a .ge. am) go to 14                                             2979
      am=a                                                               2980
      k=i                                                                2981
      go to 15                                                           2982
   14 if(a .ne. am) go to 15                                             2983
      if((x-x0)**2+(y-y0)**2.gt.(x1(k)-x0)**2+(x2(k)-y0)**2) k=i         2984
   15 continue                                                           2985
      if(k .ne. 0) go to 16                                              2986
      a0=0.0                                                             2987
      lq=1                                                               2988
      go to 2                                                            2989
   16 if(nh .ne. 0) go to 17                                             2990
      k0=k                                                               2991
      go to 20                                                           2992
   17 if(x1(k) .ne. x0) go to 18                                         2993
      a=big                                                              2994
      b=x0                                                               2995
      s=xb-b                                                             2996
      go to 19                                                           2997
   18 a=(x2(k)-y0)/(x1(k)-x0)                                            2998
      b=y0-a*x0                                                          2999
      s=yb-a*xb-b                                                        3000
   19 xh(1,nh)=a                                                         3001
      xh(2,nh)=b                                                         3002
      xh(3,nh)=sign(1.0,s)                                               3003
      if(k.eq.k0) go to 21                                               3004
   20 nh=nh+1                                                            3005
      x0=x1(k)                                                           3006
      y0=x2(k)                                                           3007
      lq=kq                                                              3008
      a0=am                                                              3009
      go to 2                                                            3010
   21 return                                                             3011
      end                                                                3012
      subroutine knts (l,nt,jv,jl,kv,nk,tb,cm,x,js)                      3013
      integer jv(l),kv(2,jl),js(*)                                       3014
      real tb(5,nk),cm(*),x(nt,l)                                        3015
      l1=0                                                               3016
      do 7 m=1,nk                                                        3017
      if(icf(m,tb,cm,jl,kv,js).eq.0) go to 7                             3018
      if(nordc(1,m,tb,cm).ne.l) go to 7                                  3019
      k=0                                                                3020
      do 1 j=1,l                                                         3021
      if(jf(m,jv(j),tb).eq.1) go to 1                                    3022
      k=1                                                                3023
      go to 2                                                            3024
    1 continue                                                           3025
    2 if(k.eq.1) go to 7                                                 3026
      ip=m                                                               3027
      l1=l1+1                                                            3028
    3 if(ip.le.0) go to 7                                                3029
      t=tb(2,ip)                                                         3030
      j=abs(t)+.1                                                        3031
      if(cm(2*j) .eq. 0.0) go to 4                                       3032
      ip=tb(4,ip)+.1                                                     3033
      go to 3                                                            3034
    4 k=1                                                                3035
    5 if(jv(k).eq.j) go to 6                                             3036
      k=k+1                                                              3037
      go to 5                                                            3038
    6 x(l1,k)=tb(3,ip)                                                   3039
      x(l1,l+k)=sign(1.0,t)                                              3040
      ip=tb(4,ip)+.1                                                     3041
      go to 3                                                            3042
    7 continue                                                           3043
      return                                                             3044
      end                                                                3045
      subroutine sclato (n,p,x,xm,xs,cm,z)                               3046
      integer p                                                          3047
      real x(n,p),xm(p),xs(p),cm(*),z(n,p)                               3048
      do 4 j=1,p                                                         3049
      j1=cm(2*j)+.1                                                      3050
      if(j1 .ne. 0) go to 2                                              3051
      if(xs(j).le.0.0) go to 4                                           3052
      do 1 i=1,n                                                         3053
      z(i,j)=xs(j)*x(i,j)+xm(j)                                          3054
    1 continue                                                           3055
      go to 4                                                            3056
    2 j1=j1-1                                                            3057
      do 3 i=1,n                                                         3058
      l=x(i,j)+.1                                                        3059
      z(i,j)=cm(l+j1)                                                    3060
    3 continue                                                           3061
    4 continue                                                           3062
      call stfmrs(0)                                                     3063
      call stcmrs(0)                                                     3064
      return                                                             3065
      end                                                                3066
      function ncat(kp)                                                  3067
      integer kp(5,*)                                                    3068
      ncat=0                                                             3069
      ll=1                                                               3070
    1 if(kp(1,ll).lt.0) go to 2                                          3071
      if(kp(1,ll).gt.0.and.kp(3,ll).le.0) ncat=ncat+1                    3072
      ll=ll+1                                                            3073
      go to 1                                                            3074
    2 return                                                             3075
      end                                                                3076
      subroutine catv (jl,kp,kv,nv,jv)                                   3077
      integer kp(5,*),kv(2,*),jv(jl,*)                                   3078
      nv=0                                                               3079
      ll=1                                                               3080
    1 if(kp(1,ll).lt.0) go to 20                                         3081
      if(kp(3,ll) .le. 0) go to 2                                        3082
      ll=ll+1                                                            3083
      go to 1                                                            3084
    2 if(kp(1,ll) .eq. jl) go to 3                                       3085
      ll=ll+1                                                            3086
      go to 1                                                            3087
    3 jg=0                                                               3088
      j=1                                                                3089
      go to 5                                                            3090
    4 j=j+1                                                              3091
    5 if((j).gt.(nv)) go to 9                                            3092
      ig=0                                                               3093
      do 6 i=1,jl                                                        3094
      if(jv(i,j).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 6                    3095
      ig=1                                                               3096
      go to 7                                                            3097
    6 continue                                                           3098
    7 if(ig .ne. 0) go to 4                                              3099
      jg=1                                                               3100
    9 if(jg .ne. 1) go to 10                                             3101
      ll=ll+1                                                            3102
      go to 1                                                            3103
   10 l1=ll+1                                                            3104
      ig=0                                                               3105
   11 if(kp(1,l1).lt.0) go to 17                                         3106
      k1=kp(1,l1)                                                        3107
      if((k1 .ne. jl) .and. (kp(3,l1) .le. 0)) go to 12                  3108
      l1=l1+1                                                            3109
      go to 11                                                           3110
   12 do 15 i=1,k1                                                       3111
      k=iabs(kv(1,kp(2,l1)+i-1))                                         3112
      do 13 j=1,jl                                                       3113
      l=iabs(kv(1,kp(2,ll)+j-1))                                         3114
      if(l .ne. k) go to 13                                              3115
      ig=1                                                               3116
      go to 14                                                           3117
   13 continue                                                           3118
   14 if(ig.eq.1) go to 16                                               3119
   15 continue                                                           3120
   16 if(ig.eq.1) go to 17                                               3121
      l1=l1+1                                                            3122
      go to 11                                                           3123
   17 if(ig .ne. 1) go to 18                                             3124
      ll=ll+1                                                            3125
      go to 1                                                            3126
   18 nv=nv+1                                                            3127
      do 19 i=1,jl                                                       3128
      jv(i,nv)=iabs(kv(1,kp(2,ll)+i-1))                                  3129
   19 continue                                                           3130
      ll=ll+1                                                            3131
      go to 1                                                            3132
   20 return                                                             3133
      end                                                                3134
      function cvlv (m,jl,jv,lv,nk,kp,kv,tb,cm,tc)                       3135
      integer jv(jl),lv(jl),kp(5,*),kv(2,*)                              3136
      real tb(5,nk),cm(*),tc(*)                                          3137
      if(m .ne. 1) go to 1                                               3138
      cvlv=cvll(jl,jv,lv,nk,tb,cm)                                       3139
      go to 2                                                            3140
    1 cvlv=cvlq(jl,jv,lv,kp,kv,cm,tc)                                    3141
    2 return                                                             3142
      end                                                                3143
      function cvll (jl,jv,lv,nk,tb,cm)                                  3144
      integer jv(jl),lv(jl),iv(2)                                        3145
      real tb(5,nk),cm(*)                                                3146
      cvll=0.0                                                           3147
      if(jl.gt.2) return                                                 3148
      do 11 m=1,nk                                                       3149
      if(tb(1,m).eq.0.0) go to 11                                        3150
      if(nordc(1,m,tb,cm).gt.0) go to 11                                 3151
      if(nord(m,tb).ne.jl) go to 11                                      3152
      call jfv(m,tb,iv)                                                  3153
      ig=0                                                               3154
      do 1 i=1,jl                                                        3155
      if(jv(i).eq.iv(i)) go to 1                                         3156
      ig=1                                                               3157
      go to 2                                                            3158
    1 continue                                                           3159
    2 if(ig.eq.1) go to 11                                               3160
      phi=1.0                                                            3161
      ip=m                                                               3162
    3 if(ip.le.0) go to 10                                               3163
      t=tb(2,ip)                                                         3164
      j=abs(t)+.1                                                        3165
      i=1                                                                3166
      go to 5                                                            3167
    4 i=i+1                                                              3168
    5 if((i).gt.(jl)) go to 6                                            3169
      if(jv(i).eq.j) go to 6                                             3170
      go to 4                                                            3171
    6 u=cm(lv(i)+int(tb(3,ip)+.1))                                       3172
      if(t .ge. 0.0) go to 8                                             3173
      if(u .ne. 0.0) go to 7                                             3174
      u=1.0                                                              3175
      go to 8                                                            3176
    7 u=0.0                                                              3177
    8 if(u .ne. 0.0) go to 9                                             3178
      phi=0.0                                                            3179
      go to 10                                                           3180
    9 ip=tb(4,ip)+.1                                                     3181
      go to 3                                                            3182
   10 if(phi.gt.0.0) cvll=cvll+tb(1,m)                                   3183
   11 continue                                                           3184
      return                                                             3185
      end                                                                3186
      function cvlq (jl,jv,lv,kp,kv,cm,tc)                               3187
      integer jv(jl),lv(jl),kp(5,*),kv(2,*)                              3188
      real cm(*),tc(*)                                                   3189
      ll=1                                                               3190
      cvlq=0.0                                                           3191
    1 if(kp(1,ll).lt.0) go to 12                                         3192
      if(kp(3,ll) .le. 0) go to 2                                        3193
      ll=ll+1                                                            3194
      go to 1                                                            3195
    2 if(kp(1,ll) .eq. jl) go to 3                                       3196
      ll=ll+1                                                            3197
      go to 1                                                            3198
    3 ig=0                                                               3199
      do 4 i=1,jl                                                        3200
      if(jv(i).eq.iabs(kv(1,kp(2,ll)+i-1))) go to 4                      3201
      ig=1                                                               3202
      go to 5                                                            3203
    4 continue                                                           3204
    5 if(ig .ne. 1) go to 6                                              3205
      ll=ll+1                                                            3206
      go to 1                                                            3207
    6 kt=1                                                               3208
      do 9 j=1,jl                                                        3209
      k=kp(2,ll)+j-1                                                     3210
      jt=cm(lv(j)+kv(2,k))+.1                                            3211
      if(kv(1,k) .ge. 0) go to 8                                         3212
      if(jt .ne. 0) go to 7                                              3213
      jt=1                                                               3214
      go to 8                                                            3215
    7 jt=0                                                               3216
    8 if(jt .ne. 0) go to 9                                              3217
      kt=0                                                               3218
      go to 10                                                           3219
    9 continue                                                           3220
   10 if(kt .ne. 1) go to 11                                             3221
      cvlq=cvlq+tc(-kp(3,ll))                                            3222
   11 ll=ll+1                                                            3223
      go to 1                                                            3224
   12 return                                                             3225
      end                                                                3226
      subroutine ccoll (nk,tb,cm,kp,kv,lp,lv,jv)                         3227
      integer kp(5,*),kv(2,*),lp(3,*),lv(*),jv(*)                        3228
      real tb(5,*),cm(*)                                                 3229
      call collc(nk,tb,cm,kp,kv,jv)                                      3230
      call purcat(nk,tb,cm,kp,kv,li,jv)                                  3231
      ll=li+1                                                            3232
      l1=1                                                               3233
      l2=l1                                                              3234
    1 if(kp(1,ll).lt.0) go to 2                                          3235
      kp(4,ll)=l1                                                        3236
      call collf(nk,tb,cm,kp(1,ll),kv(1,kp(2,ll)),l1,l2,lp,lv,jv)        3237
      kp(3,ll)=l1-kp(4,ll)                                               3238
      ll=ll+1                                                            3239
      go to 1                                                            3240
    2 lp(1,l1)=0                                                         3241
      return                                                             3242
      end                                                                3243
      subroutine collc (nk,tb,cm,kp,kv,jv)                               3244
      integer kp(5,*),kv(2,*),jv(*)                                      3245
      real tb(5,nk),cm(*)                                                3246
      kp(1,1)=0                                                          3247
      kp(2,1)=1                                                          3248
      l1=2                                                               3249
      l2=1                                                               3250
      mc=0                                                               3251
      do 1 m=1,nk                                                        3252
      if(tb(1,m).ne.0.0) mc=max0(mc,nordc(2,m,tb,cm))                    3253
    1 continue                                                           3254
      mt=1                                                               3255
      go to 3                                                            3256
    2 mt=mt+1                                                            3257
    3 if((mt).gt.(mc)) go to 18                                          3258
      l10=l1                                                             3259
      do 17 m=1,nk                                                       3260
      if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.mt) go to 17              3261
      call jfvc(2,m,tb,cm,nv,jv,jv(mt+1))                                3262
      jg=0                                                               3263
      l1m1=l1-1                                                          3264
      i=l10                                                              3265
      go to 5                                                            3266
    4 i=i+1                                                              3267
    5 if((i).gt.(l1m1)) go to 15                                         3268
      k=kp(2,i)-1                                                        3269
      ig=0                                                               3270
      do 6 j=1,mt                                                        3271
      if(iabs(jv(j)).eq.iabs(kv(1,k+j))) go to 6                         3272
      ig=1                                                               3273
      go to 7                                                            3274
    6 continue                                                           3275
    7 if(ig .ne. 0) go to 13                                             3276
      do 12 j=1,mt                                                       3277
      m1=kv(2,k+j)                                                       3278
      m2=jv(mt+j)                                                        3279
      jj=iabs(jv(j))                                                     3280
      nc=int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1                           3281
      kk=jv(j)*kv(1,k+j)                                                 3282
      do 10 jk=1,nc                                                      3283
      z=cm(jk+m2)                                                        3284
      if(kk .ge. 0) go to 9                                              3285
      if(z .ne. 0.0) go to 8                                             3286
      z=1.0                                                              3287
      go to 9                                                            3288
    8 z=0.0                                                              3289
    9 if(cm(jk+m1).eq.z) go to 10                                        3290
      ig=1                                                               3291
      go to 11                                                           3292
   10 continue                                                           3293
   11 if(ig.eq.1) go to 13                                               3294
   12 continue                                                           3295
   13 if(ig .ne. 0) go to 4                                              3296
      jg=1                                                               3297
   15 if(jg .ne. 0) go to 17                                             3298
      kp(1,l1)=mt                                                        3299
      kp(2,l1)=l2                                                        3300
      k=l2-1                                                             3301
      do 16 i=1,mt                                                       3302
      kv(1,i+k)=jv(i)                                                    3303
      kv(2,i+k)=jv(i+mt)                                                 3304
   16 continue                                                           3305
      l1=l1+1                                                            3306
      l2=l2+mt                                                           3307
   17 continue                                                           3308
      go to 2                                                            3309
   18 kp(1,l1)=-1                                                        3310
      return                                                             3311
      end                                                                3312
      function nordc (l,m,tb,cm)                                         3313
      real tb(5,*),cm(*)                                                 3314
      ip=m                                                               3315
      nordc=0                                                            3316
    1 if(ip.le.0) go to 4                                                3317
      j=abs(tb(2,ip))+.1                                                 3318
      if(l .ne. 1) go to 2                                               3319
      if(cm(2*j).eq.0.0) nordc=nordc+1                                   3320
      go to 3                                                            3321
    2 if(cm(2*j).gt.0.0) nordc=nordc+1                                   3322
    3 ip=tb(4,ip)+.1                                                     3323
      go to 1                                                            3324
    4 return                                                             3325
      end                                                                3326
      subroutine jfvc (l,m,tb,cm,nv,jv,jp)                               3327
      integer jv(*),jp(*)                                                3328
      real tb(5,*),cm(*)                                                 3329
      ip=m                                                               3330
      nv=0                                                               3331
    1 if(ip.le.0) go to 5                                                3332
      j=abs(tb(2,ip))+.1                                                 3333
      if(l .ne. 1) go to 3                                               3334
      if(cm(2*j) .le. 0.0) go to 4                                       3335
      ip=tb(4,ip)+.1                                                     3336
      go to 1                                                            3337
    3 if(cm(2*j) .ne. 0.0) go to 4                                       3338
      ip=tb(4,ip)+.1                                                     3339
      go to 1                                                            3340
    4 nv=nv+1                                                            3341
      jv(nv)=j                                                           3342
      if(l.ne.1.and.tb(2,ip).lt.0.0) jv(nv)=-j                           3343
      if(l.ne.1) jp(nv)=tb(3,ip)+.1                                      3344
      ip=tb(4,ip)+.1                                                     3345
      go to 1                                                            3346
    5 if(nv.le.1) return                                                 3347
      j=nv-1                                                             3348
    6 ll=0                                                               3349
      do 7 i=1,j                                                         3350
      if(iabs(jv(i)) .le. iabs(jv(i+1))) go to 7                         3351
      ll=1                                                               3352
      k=jv(i)                                                            3353
      jv(i)=jv(i+1)                                                      3354
      jv(i+1)=k                                                          3355
      if(l .eq. 1) go to 7                                               3356
      k=jp(i)                                                            3357
      jp(i)=jp(i+1)                                                      3358
      jp(i+1)=k                                                          3359
    7 continue                                                           3360
      if(ll.eq.0) go to 8                                                3361
      go to 6                                                            3362
    8 return                                                             3363
      end                                                                3364
      subroutine purcat (nk,tb,cm,kp,kv,li,jv)                           3365
      integer kp(5,*),kv(2,*),jv(*)                                      3366
      real tb(5,nk),cm(*)                                                3367
      lm=1                                                               3368
    1 if(kp(1,lm).lt.0) go to 2                                          3369
      lm=lm+1                                                            3370
      go to 1                                                            3371
    2 ll=1                                                               3372
      li=0                                                               3373
    3 if(kp(1,ll).lt.0) go to 20                                         3374
      jl=kp(1,ll)                                                        3375
      if(jl .gt. 0) go to 4                                              3376
      ll=ll+1                                                            3377
      go to 3                                                            3378
    4 ifg=0                                                              3379
      jfg=ifg                                                            3380
      do 6 m=1,nk                                                        3381
      if(icf(m,tb,cm,jl,kv(1,kp(2,ll)),jv).eq.0) go to 6                 3382
      if(nord(m,tb) .ne. jl) go to 5                                     3383
      ifg=1                                                              3384
      go to 6                                                            3385
    5 jfg=1                                                              3386
    6 continue                                                           3387
      if(ifg .ne. 0) go to 9                                             3388
      if(jfg .ne. 0) go to 8                                             3389
c     write(6,7)                                                         3390
    7 format (' bug in purcat - term not found.')                        3391
      stop                                                               3392
    8 ll=ll+1                                                            3393
      go to 3                                                            3394
    9 li=li+1                                                            3395
      j=lm                                                               3396
      go to 11                                                           3397
   10 j=j+(-1)                                                           3398
   11 if((-1)*((j)-(li)).gt.0) go to 13                                  3399
      do 12 i=1,5                                                        3400
      kp(i,j+1)=kp(i,j)                                                  3401
   12 continue                                                           3402
      go to 10                                                           3403
   13 lm=lm+1                                                            3404
      ll=ll+1                                                            3405
      do 14 i=1,5                                                        3406
      kp(i,li)=kp(i,ll)                                                  3407
   14 continue                                                           3408
      kp(3,li)=0                                                         3409
      kp(4,li)=1                                                         3410
      kp(5,li)=0                                                         3411
      if(jfg .ne. 1) go to 15                                            3412
      ll=ll+1                                                            3413
      go to 3                                                            3414
   15 j=ll+1                                                             3415
      go to 17                                                           3416
   16 j=j+1                                                              3417
   17 if((j).gt.(lm)) go to 19                                           3418
      do 18 i=1,5                                                        3419
      kp(i,j-1)=kp(i,j)                                                  3420
   18 continue                                                           3421
      go to 16                                                           3422
   19 lm=lm-1                                                            3423
      go to 3                                                            3424
   20 return                                                             3425
      end                                                                3426
      subroutine collf (nk,tb,cm,jl,kv,l1,l2,lp,lv,jv)                   3427
      integer kv(2,*),lp(3,*),lv(*),jv(*)                                3428
      real tb(5,*),cm(*)                                                 3429
      mo=0                                                               3430
      do 1 m=1,nk                                                        3431
      if(icf(m,tb,cm,jl,kv,jv).ne.0) mo=max0(mo,nordc(1,m,tb,cm))        3432
    1 continue                                                           3433
      if(mo.eq.0) return                                                 3434
      do 10 mt=1,mo                                                      3435
      l10=l1                                                             3436
      do 9 m=1,nk                                                        3437
      if(icf(m,tb,cm,jl,kv,jv).eq.0) go to 9                             3438
      if(nordc(1,m,tb,cm).ne.mt) go to 9                                 3439
      call jfvc(1,m,tb,cm,nv,jv,jv)                                      3440
      jg=0                                                               3441
      l1m1=l1-1                                                          3442
      i=l10                                                              3443
      go to 3                                                            3444
    2 i=i+1                                                              3445
    3 if((i).gt.(l1m1)) go to 7                                          3446
      k=lp(2,i)-1                                                        3447
      ig=0                                                               3448
      do 4 j=1,mt                                                        3449
      if(jv(j).eq.lv(k+j)) go to 4                                       3450
      ig=1                                                               3451
      go to 5                                                            3452
    4 continue                                                           3453
    5 if(ig .ne. 0) go to 2                                              3454
      jg=1                                                               3455
      lp(3,i)=lp(3,i)+1                                                  3456
    7 if(jg .ne. 0) go to 9                                              3457
      lp(1,l1)=mt                                                        3458
      lp(2,l1)=l2                                                        3459
      lp(3,l1)=1                                                         3460
      k=l2-1                                                             3461
      do 8 i=1,mt                                                        3462
      lv(i+k)=jv(i)                                                      3463
    8 continue                                                           3464
      l1=l1+1                                                            3465
      l2=l2+mt                                                           3466
    9 continue                                                           3467
   10 continue                                                           3468
      return                                                             3469
      end                                                                3470
      function icf (m,tb,cm,jl,kv,jv)                                    3471
      integer kv(2,jl),jv(*)                                             3472
      real tb(5,*),cm(*)                                                 3473
      icf=0                                                              3474
      if(tb(1,m).eq.0.0.or.nordc(2,m,tb,cm).ne.jl) return                3475
      if(jl .ne. 0) go to 1                                              3476
      icf=1                                                              3477
      return                                                             3478
    1 call jfvc(2,m,tb,cm,nv,jv,jv(jl+1))                                3479
      do 2 j=1,jl                                                        3480
      if(iabs(jv(j)).ne.iabs(kv(1,j))) return                            3481
    2 continue                                                           3482
      do 6 j=1,jl                                                        3483
      l1=kv(2,j)                                                         3484
      l2=jv(jl+j)                                                        3485
      k=2*iabs(jv(j))                                                    3486
      kk=jv(j)*kv(1,j)                                                   3487
      nc=int(cm(k+1)+.1)-int(cm(k)+.1)+1                                 3488
      do 5 i=1,nc                                                        3489
      z=cm(i+l2)                                                         3490
      if(kk .ge. 0) go to 4                                              3491
      if(z .ne. 0.0) go to 3                                             3492
      z=1.0                                                              3493
      go to 4                                                            3494
    3 z=0.0                                                              3495
    4 if(cm(i+l1).ne.z) return                                           3496
    5 continue                                                           3497
    6 continue                                                           3498
      icf=1                                                              3499
      return                                                             3500
      end                                                                3501
      function icat (x,j,cm)                                             3502
      real cm(*)                                                         3503
      j0=cm(2*j)+.1                                                      3504
      j1=j0                                                              3505
      j2=cm(2*j+1)+.1                                                    3506
    1 if(j2.eq.j1+1) go to 5                                             3507
      k=(j1+j2)/2                                                        3508
      if(cm(k) .ne. x) go to 2                                           3509
      icat=k-j0+1                                                        3510
      return                                                             3511
    2 if(cm(k) .ge. x) go to 3                                           3512
      j1=k                                                               3513
      go to 1                                                            3514
    3 j2=k                                                               3515
      go to 1                                                            3516
    5 if(x .ne. cm(j1)) go to 6                                          3517
      icat=j1-j0+1                                                       3518
      go to 8                                                            3519
    6 if(x .ne. cm(j2)) go to 7                                          3520
      icat=j2-j0+1                                                       3521
      go to 8                                                            3522
    7 icat=0                                                             3523
    8 return                                                             3524
      end                                                                3525
      subroutine csp  (jp,nc,m,n,x,y,w,nk,tb,cm,kcp,yb,d,kr,ntt,sw,me,mk 3526
     1p2,nop,sc,db,sp,mm)                                                3527
      integer mm(nc,2)                                                   3528
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(n)                         3529
      double precision yb,sw,d(nk,*),db(n,*),sp(mkp2,*),a,b,s,eps,dv,dy  3530
      data eps,big /1.d-4,9.9e30/                                        3531
      nop=0                                                              3532
      if(nc .gt. 1) go to 1                                              3533
      tb(1,m)=big                                                        3534
      return                                                             3535
    1 mk=mkp2-2                                                          3536
      mkp1=mk+1                                                          3537
      n1=nc+1                                                            3538
      do 3 j=1,n1                                                        3539
      do 2 i=1,mkp2                                                      3540
      sp(i,j)=0.d0                                                       3541
    2 continue                                                           3542
    3 continue                                                           3543
      do 4 j=1,nc                                                        3544
      mm(j,2)=0                                                          3545
    4 continue                                                           3546
      do 7 i=1,n                                                         3547
      h=sc(i)                                                            3548
      if(h.le.0.0.or.w(i).le.0.0) go to 7                                3549
      wh=w(i)*h                                                          3550
      k=x(i,jp)+.1                                                       3551
      mm(k,2)=mm(k,2)+1                                                  3552
      sp(mkp2,k)=sp(mkp2,k)+wh                                           3553
      sp(mkp1,k)=sp(mkp1,k)+wh*(y(i)-yb)                                 3554
      sp(m,k)=sp(m,k)+wh*h                                               3555
      j=1                                                                3556
      go to 6                                                            3557
    5 j=j+1                                                              3558
    6 if((j).gt.(kr)) go to 7                                            3559
      sp(j,k)=sp(j,k)+wh*db(i,j)                                         3560
      go to 5                                                            3561
    7 continue                                                           3562
      do 8 j=1,nc                                                        3563
      mm(j,1)=j                                                          3564
    8 continue                                                           3565
      bof0=big                                                           3566
      ns=0                                                               3567
      jj=nc                                                              3568
      nrt=0                                                              3569
      k1=1                                                               3570
    9 bof1=big                                                           3571
      js=0                                                               3572
      do 18 j=1,jj                                                       3573
      k=mm(j,1)                                                          3574
      if(mm(k,2).eq.0) go to 18                                          3575
      nr=nrt+mm(k,2)                                                     3576
      if(nr.le.me.or.ntt-nr.le.me) go to 18                              3577
      dy=sp(mkp1,n1)+sp(mkp1,k)                                          3578
      a=sp(mkp2,n1)+sp(mkp2,k)                                           3579
      dv=sp(m,n1)+sp(m,k)-a**2/sw                                        3580
      if(dv .le. 0.d0) go to 17                                          3581
      i=1                                                                3582
      go to 11                                                           3583
   10 i=i+1                                                              3584
   11 if((i).gt.(kr)) go to 12                                           3585
      d(i,2)=sp(i,n1)+sp(i,k)                                            3586
      go to 10                                                           3587
   12 a=0.d0                                                             3588
      b=a                                                                3589
      i=1                                                                3590
      go to 14                                                           3591
   13 i=i+1                                                              3592
   14 if((i).gt.(kr)) go to 15                                           3593
      s=d(i,2)                                                           3594
      a=a+s*d(i,1)                                                       3595
      b=b+s**2                                                           3596
      go to 13                                                           3597
   15 b=dv-b                                                             3598
      if(b .le. eps*dv) go to 17                                         3599
      nop=nop+1                                                          3600
      b=-(dy-a)**2/b                                                     3601
      if(b .ge. bof1) go to 16                                           3602
      bof1=b                                                             3603
      js=j                                                               3604
   16 if(b .ge. bof0) go to 17                                           3605
      bof0=b                                                             3606
      ns=jj                                                              3607
   17 if(nc.eq.2) go to 19                                               3608
   18 continue                                                           3609
   19 if(js.eq.0) go to 23                                               3610
      k=mm(js,1)                                                         3611
      mm(js,1)=mm(jj,1)                                                  3612
      mm(jj,1)=k                                                         3613
      sp(mkp1,n1)=sp(mkp1,n1)+sp(mkp1,k)                                 3614
      sp(mkp2,n1)=sp(mkp2,n1)+sp(mkp2,k)                                 3615
      nrt=nrt+mm(k,2)                                                    3616
      sp(m,n1)=sp(m,n1)+sp(m,k)                                          3617
      i=1                                                                3618
      go to 21                                                           3619
   20 i=i+1                                                              3620
   21 if((i).gt.(kr)) go to 22                                           3621
      sp(i,n1)=sp(i,n1)+sp(i,k)                                          3622
      go to 20                                                           3623
   22 jj=jj-1                                                            3624
      if(jj.le.2) go to 23                                               3625
      go to 9                                                            3626
   23 tb(1,m)=bof0                                                       3627
      tb(3,m)=kcp                                                        3628
      do 24 j=1,nc                                                       3629
      cm(j+kcp)=0.0                                                      3630
   24 continue                                                           3631
      if(ns.eq.0) return                                                 3632
      do 25 j=ns,nc                                                      3633
      cm(mm(j,1)+kcp)=1.0                                                3634
   25 continue                                                           3635
      return                                                             3636
      end                                                                3637
      subroutine rspnpr (it,il,n,y,w,m)                                  3638
      integer m(n)                                                       3639
      real y(n),w(n),wm(2)                                               3640
      if(it.le.0) return                                                 3641
      if(il .ne. 1) go to 2                                              3642
      wm(1)=0.0                                                          3643
      wm(2)=wm(1)                                                        3644
      do 1 i=1,n                                                         3645
      k=y(i)+1.1                                                         3646
      wm(k)=wm(k)+w(i)                                                   3647
    1 continue                                                           3648
      wt=wm(1)+wm(2)                                                     3649
      wm(1)=wm(1)/wt                                                     3650
      wm(2)=wm(2)/wt                                                     3651
c     write(it,'(/,'' binary (0/1) response:  mass(0) ='',g12.4,         3652
c    1                           ''   mass(1) ='',g12.4)') wm(1),wm(2)   3653
      return                                                             3654
    2 continue
c     write(it,'(/,'' ordinal response:'')')                             3655
c     write(it,'(''      min         n/4         n/2        3n/4         3656
c    1 max'')')                                                          3657
      do 3 i=1,n                                                         3658
      m(i)=i                                                             3659
    3 continue                                                           3660
      call psort(y,m,1,n)                                                3661
c     write(it,'('' '',5g12.4)') y(m(1)),y(m(n/4)),y(m(n/2)),y(m(n-n/4)) 3662
c    1,y(m(n))                                                           3663
      return                                                             3664
      end                                                                3665
      subroutine ordpr (it,n,p,x,lx,m)                                   3666
      integer p,lx(p),m(n,p)                                             3667
      real x(n,p)                                                        3668
      if(it.le.0) return                                                 3669
      no=0                                                               3670
      do 1 j=1,p                                                         3671
      if(lx(j).gt.0) no=no+1                                             3672
    1 continue                                                           3673
      if(no.eq.0) return                                                 3674
c     write(it,'(/,'' there are'',i3,'' ordinal predictor variables.'',/ 3675
c    1)') no                                                             3676
c     write(it,'(''  var     min         n/4         n/2        3n/4     3677
c    1     max'')')                                                      3678
      n1=n/4                                                             3679
      n2=n/2                                                             3680
      n3=n-n1                                                            3681
      do 2 j=1,p                                                         3682
      if(lx(j).le.0) go to 2                                             3683
c     write(it,'('' '',i3,'' '',5g12.4)') j,  x(m(1,j),j),x(m(n1,j),j),x 3684
c    1(m(n2,j),j),x(m(n3,j),j),x(m(n,j),j)                               3685
    2 continue                                                           3686
      return                                                             3687
      end                                                                3688
      subroutine catpr(it,n,p,x,cm,mm)                                   3689
      integer p,mm(*)                                                    3690
      real x(n,p),cm(*)                                                  3691
      if(it.le.0) return                                                 3692
      nct=cm(1)+.1                                                       3693
      if(nct.eq.0) return                                                3694
      n2=2*p+1                                                           3695
      np=0                                                               3696
c     write(it,'(/,'' there are'',i3,'' categorical predictor variables. 3697
c    1'')') nct                                                          3698
      i=2                                                                3699
      go to 2                                                            3700
    1 i=i+(2)                                                            3701
    2 if((2)*((i)-(n2)).gt.0) go to 6                                    3702
      np=np+1                                                            3703
      j1=cm(i)+.1                                                        3704
      if(j1.eq.0) go to 1                                                3705
      j2=cm(i+1)+.1                                                      3706
      nv=j2-j1+1                                                         3707
      do 3 j=1,nv                                                        3708
      mm(j)=0                                                            3709
    3 continue                                                           3710
      do 4 j=1,n                                                         3711
      ic=x(j,np)+.1                                                      3712
      mm(ic)=mm(ic)+1                                                    3713
    4 continue                                                           3714
c     write(it,'(/,'' categorical variable'',i3,'' has'',i3,'' values.'' 3715
c    1)') np,nv                                                          3716
c     write(it,'(''  value     internal code     counts'')')             3717
      k=0                                                                3718
      do 5 j=j1,j2                                                       3719
      k=k+1                                                              3720
c     write(it,'(f6.0,i13,i15)') cm(j),k,mm(k)                           3721
    5 continue                                                           3722
      go to 1                                                            3723
    6 return                                                             3724
      end                                                                3725
      subroutine holl (jp,cm,t,h)                                        3726
      real cm(*)                                                         3727
      character*28 h                                                     3728
      j1=cm(2*jp)+.1                                                     3729
      j2=cm(2*jp+1)+.1                                                   3730
      j2=j2-j1+1                                                         3731
      if(j2 .le. 28) go to 1                                             3732
      h='   cat. factor > 28 values  '                                   3733
      return                                                             3734
    1 h='                            '                                   3735
      j1=(28-j2)/2                                                       3736
      j2=j1+j2-1                                                         3737
      k=t+.1                                                             3738
      do 3 j=j1,j2                                                       3739
      if(cm(k+j-j1+1) .le. 0.0) go to 2                                  3740
      h(j:j)='1'                                                         3741
      go to 3                                                            3742
    2 h(j:j)='0'                                                         3743
    3 continue                                                           3744
      return                                                             3745
      end                                                                3746
      subroutine slova (nk,it,tb,ni,lp,lv)                               3747
      integer lp(3,*),lv(*)                                              3748
      real tb(5,nk)                                                      3749
c     write(it,4) ni                                                     3750
      call coll(nk,tb,lp,lv,lp(1,nk+1))                                  3751
      m=1                                                                3752
    1 if(lp(1,m).eq.0) go to 2                                           3753
      m=m+1                                                              3754
      go to 1                                                            3755
    2 na=m-1                                                             3756
      do 3 m=1,na                                                        3757
      k2=lp(2,m)                                                         3758
      i2=lp(1,m)+k2-1                                                    3759
c     write(it,5) m,lp(3,m),(lv(i),i=k2,i2)                              3760
    3 continue                                                           3761
      return                                                             3762
    4 format(/,' sliced anova decomposition on',i3,' basis functions:',/ 3763
     1,  '   fun      #bsfns      variable(s)')                          3764
    5 format('  ',i3,'         ',i2,'       ',20i4)                      3765
      end                                                                3766
      subroutine reducq (flg,x,nk,tb,cm,tc,kp,kv,lp,lv,r,td,sc,fc)       3767
      integer kp(5,*),kv(2,*),lp(3,*),lv(*)                              3768
      real x(*),tb(5,nk),cm(*),tc(*),r(*),td(2,nk),sc(2,*),fc(*)         3769
      ll=1                                                               3770
      la=ll                                                              3771
      l1=la                                                              3772
      laa=0                                                              3773
      do 1 m=1,nk                                                        3774
      td(1,m)=0.0                                                        3775
    1 continue                                                           3776
    2 if(kp(1,ll).lt.0) go to 9                                          3777
      nv=0                                                               3778
      if(kp(1,ll) .le. 0) go to 4                                        3779
      jl=kp(1,ll)                                                        3780
      do 3 il=1,jl                                                       3781
      k=kp(2,ll)+il-1                                                    3782
      nv=nv+1                                                            3783
      sc(1,nv)=kv(1,k)                                                   3784
      sc(2,nv)=kv(2,k)                                                   3785
    3 continue                                                           3786
      go to 5                                                            3787
    4 if(kp(3,ll) .gt. 0) go to 5                                        3788
      ll=ll+1                                                            3789
      go to 2                                                            3790
    5 if(kp(3,ll) .gt. 0) go to 6                                        3791
      m=match(nv,sc,nk,tb,cm,r,0)                                        3792
      td(1,m)=tc(-kp(3,ll))                                              3793
      ll=ll+1                                                            3794
      go to 2                                                            3795
    6 kp3=kp(3,ll)                                                       3796
      do 8 k=1,kp3                                                       3797
      l=lp(1,l1)                                                         3798
      nt=lp(3,l1)                                                        3799
      laa=laa+5*l*nt                                                     3800
      do 7 jp=1,nt                                                       3801
      call gtrm(1,jp,l,nt,lv(lp(2,l1)),flg,x,nk,tb,tc(la),sc(1,nv+1),fc) 3802
      m=match(nv+l,sc,nk,tb,cm,r,0)                                      3803
      td(1,m)=tc(jp+laa)                                                 3804
      call std(m,flg,x,l,sc(1,nv+1),fc,nk,tb,r,td)                       3805
    7 continue                                                           3806
      laa=laa+nt                                                         3807
      l1=l1+1                                                            3808
      la=la+nt*(5*l+1)                                                   3809
    8 continue                                                           3810
      ll=ll+1                                                            3811
      go to 2                                                            3812
    9 return                                                             3813
      end                                                                3814
      subroutine gtrm (il,jp,l,nt,jv,flg,x,nk,tb,tc,te,fc)               3815
      integer jv(l)                                                      3816
      real x(*),tb(5,nk),tc(nt,*),te(2,*),fc(*)                          3817
      l2=l+l                                                             3818
      nf=0                                                               3819
      l3=l2+l                                                            3820
      l4=l3+l                                                            3821
      do 1 k=1,l                                                         3822
      j=jv(k)                                                            3823
      jj=j                                                               3824
      if(tc(jp,k+l).gt.tc(jp,k+l2)) jj=-jj                               3825
      te(1,k)=jj                                                         3826
      te(2,k)=tc(jp,k)                                                   3827
      if(il.eq.2) go to 1                                                3828
      if(x(j).eq.flg) go to 1                                            3829
      nf=nf+1                                                            3830
      fc(nf)=cue(x(j),tc(jp,k+l),tc(jp,k),tc(jp,k+l2),tc(jp,k+l3),tc(jp, 3831
     1k+l4))                                                             3832
    1 continue                                                           3833
      return                                                             3834
      end                                                                3835
      function match (nv,te,nk,tb,cm,r,iz)                               3836
      real te(2,nv),tb(5,nk),cm(*),r(*)                                  3837
      match=0                                                            3838
      do 15 m=1,nk                                                       3839
      if(tb(1,m).eq.0.0) go to 15                                        3840
      if(nord(m,tb).ne.nv) go to 15                                      3841
      jg=0                                                               3842
      do 13 j=1,nv                                                       3843
      t=te(1,j)                                                          3844
      u=te(2,j)                                                          3845
      jp=abs(t)+.1                                                       3846
      jp2=2*jp                                                           3847
      jp21=jp2+1                                                         3848
      jpp=jp                                                             3849
      if(t.lt.0.0) jpp=-jpp                                              3850
      ig=0                                                               3851
      ip=m                                                               3852
    1 if(ip.le.0) go to 12                                               3853
      t=tb(2,ip)                                                         3854
      jq=abs(t)+.1                                                       3855
      jqq=jq                                                             3856
      if(t.lt.0.0) jqq=-jqq                                              3857
      if(jp .eq. jq) go to 2                                             3858
      ip=tb(4,ip)+.1                                                     3859
      go to 1                                                            3860
    2 if(cm(jp2) .ne. 0.0) go to 4                                       3861
      if(jpp .ne. jqq .or. ieq(tb(3,ip),u,r(jp)) .ne. 1) go to 3         3862
      ig=1                                                               3863
      go to 12                                                           3864
    3 ip=tb(4,ip)+.1                                                     3865
      go to 1                                                            3866
    4 nc=int(cm(jp21)+.1)-int(cm(jp2)+.1)+1                              3867
      i1=u+.1                                                            3868
      i2=tb(3,ip)+.1                                                     3869
      kg=0                                                               3870
      do 9 i=1,nc                                                        3871
      j1=cm(i1+i)                                                        3872
      j2=cm(i2+i)                                                        3873
      if(jpp .ge. 0) go to 6                                             3874
      if(j1 .ne. 0) go to 5                                              3875
      j1=1                                                               3876
      go to 6                                                            3877
    5 j1=0                                                               3878
    6 if(jqq .ge. 0) go to 8                                             3879
      if(j2 .ne. 0) go to 7                                              3880
      j2=1                                                               3881
      go to 8                                                            3882
    7 j2=0                                                               3883
    8 if(j1 .eq. j2) go to 9                                             3884
      kg=1                                                               3885
      go to 10                                                           3886
    9 continue                                                           3887
   10 if(kg .ne. 0) go to 11                                             3888
      ig=1                                                               3889
      go to 12                                                           3890
   11 ip=tb(4,ip)+.1                                                     3891
      go to 1                                                            3892
   12 if(ig .ne. 0) go to 13                                             3893
      jg=1                                                               3894
      go to 14                                                           3895
   13 continue                                                           3896
   14 if(jg .ne. 0) go to 15                                             3897
      match=m                                                            3898
      go to 16                                                           3899
   15 continue                                                           3900
   16 if(match.gt.0.or.iz.ne.0) return                                   3901
c     write(6,17)                                                        3902
   17 format (' bug in match - term not found.')                         3903
      do 19 j=1,nv                                                       3904
c     write(6,18)j,te(1,j),te(2,j)                                       3905
   18 format (' te(',i2,')=',2g12.4)                                     3906
   19 continue                                                           3907
      do 21 j=1,nk                                                       3908
c     write(6,20)j,(tb(i,j),i=1,4)                                       3909
   20 format (' tb(',i2,')=',4g12.4)                                     3910
   21 continue                                                           3911
      stop                                                               3912
      end                                                                3913
      subroutine std (m,flg,x,l,te,fc,nk,tb,r,td)                        3914
      real x(*),te(2,*),fc(*),tb(5,nk),r(*),td(2,nk)                     3915
      ip=m                                                               3916
    1 if(ip.le.0) go to 4                                                3917
      t=tb(2,ip)                                                         3918
      j=abs(t)+.1                                                        3919
      jj=j                                                               3920
      if(t.lt.0.0) jj=-jj                                                3921
      u=tb(3,ip)                                                         3922
      k=0                                                                3923
      ig=0                                                               3924
      do 2 i=1,l                                                         3925
      t=te(1,i)                                                          3926
      jp=abs(t)+.1                                                       3927
      if(x(jp).ne.flg) k=k+1                                             3928
      if(t.lt.0.0) jp=-jp                                                3929
      if(jj .ne. jp .or. ieq(te(2,i),u,r(j)) .ne. 1) go to 2             3930
      ig=1                                                               3931
      go to 3                                                            3932
    2 continue                                                           3933
    3 if(ig.eq.1.and.x(j).ne.flg) td(2,ip)=fc(k)                         3934
      ip=tb(4,ip)+.1                                                     3935
      go to 1                                                            3936
    4 return                                                             3937
      end                                                                3938
      subroutine reducl (flg,x,nk,az,tb,cm,bz,td,r,azn,tbn,bzn,sc)       3939
      real x(*),tb(5,nk),cm(*),td(2,nk),r(*),tbn(5,nk),sc(*)             3940
      azn=az                                                             3941
      do 2 m=1,nk                                                        3942
      do 1 i=1,5                                                         3943
      tbn(i,m)=tb(i,m)                                                   3944
    1 continue                                                           3945
    2 continue                                                           3946
      bzn=bz                                                             3947
      do 9 m=1,nk                                                        3948
      t=tb(2,m)                                                          3949
      j=abs(t)+.1                                                        3950
      if(x(j).eq.flg) go to 9                                            3951
      if(cm(2*j) .le. 0.0) go to 7                                       3952
      k=icat(x(j),j,cm)                                                  3953
      if(k .ne. 0) go to 3                                               3954
      u=0.0                                                              3955
      go to 4                                                            3956
    3 u=cm(k+int(tb(3,m)+.1))                                            3957
    4 if(t .ge. 0.0) go to 6                                             3958
      if(u .ne. 0.0) go to 5                                             3959
      u=1.0                                                              3960
      go to 6                                                            3961
    5 u=0.0                                                              3962
    6 td(2,m)=u                                                          3963
      go to 8                                                            3964
    7 u=amax1(0.0,sign(1.0,t)*(x(j)-tb(3,m)))                            3965
    8 sc(m)=u                                                            3966
    9 continue                                                           3967
      m=nk                                                               3968
      go to 11                                                           3969
   10 m=m+(-1)                                                           3970
   11 if((-1)*((m)-(1)).gt.0) go to 21                                   3971
      ip=tbn(4,m)+.1                                                     3972
      t=tbn(2,m)                                                         3973
      j=abs(t)+.1                                                        3974
      if(x(j) .ne. flg) go to 15                                         3975
      if(tbn(1,m) .eq. 0.0) go to 10                                     3976
      iq=ip                                                              3977
   12 if(iq.le.0) go to 10                                               3978
      t=tbn(2,iq)                                                        3979
      j=abs(t)+.1                                                        3980
      if(x(j) .eq. flg) go to 13                                         3981
      tbn(1,m)=tbn(1,m)*sc(iq)                                           3982
      td(1,m)=td(1,m)*td(2,iq)                                           3983
   13 iq=tbn(4,iq)+.1                                                    3984
      go to 12                                                           3985
   15 k=m+1                                                              3986
      go to 17                                                           3987
   16 k=k+1                                                              3988
   17 if((k).gt.(nk)) go to 18                                           3989
      if(int(tbn(4,k)+.1).eq.m) tbn(4,k)=tbn(4,m)                        3990
      go to 16                                                           3991
   18 if(tbn(1,m).eq.0.0) go to 10                                       3992
      if(ip .ne. 0) go to 19                                             3993
      azn=azn+tbn(1,m)*sc(m)                                             3994
      bzn=bzn+td(1,m)*td(2,m)                                            3995
      go to 10                                                           3996
   19 tbn(1,ip)=tbn(1,ip)+tbn(1,m)*sc(m)                                 3997
      td(1,ip)=td(1,ip)+td(1,m)*td(2,m)                                  3998
      go to 10                                                           3999
   21 no=nk                                                              4000
      m=nk                                                               4001
      go to 23                                                           4002
   22 m=m+(-1)                                                           4003
   23 if((-1)*((m)-(1)).gt.0) go to 31                                   4004
      t=tb(2,m)                                                          4005
      j=abs(t)+.1                                                        4006
      if(x(j).eq.flg) go to 22                                           4007
      k=m+1                                                              4008
      go to 25                                                           4009
   24 k=k+1                                                              4010
   25 if((k).gt.(no)) go to 30                                           4011
      td(1,k-1)=td(1,k)                                                  4012
      do 26 i=1,5                                                        4013
      tbn(i,k-1)=tbn(i,k)                                                4014
   26 continue                                                           4015
      i=k+1                                                              4016
      go to 28                                                           4017
   27 i=i+1                                                              4018
   28 if((i).gt.(no)) go to 24                                           4019
      if(int(tbn(4,i)+.1).eq.k) tbn(4,i)=k-1                             4020
      go to 27                                                           4021
   30 no=no-1                                                            4022
      go to 22                                                           4023
   31 m=no+1                                                             4024
      go to 33                                                           4025
   32 m=m+1                                                              4026
   33 if((m).gt.(nk)) go to 34                                           4027
      tbn(1,m)=0.0                                                       4028
      go to 32                                                           4029
   34 m=no                                                               4030
      go to 36                                                           4031
   35 m=m+(-1)                                                           4032
   36 if((-1)*((m)-(2)).gt.0) go to 39                                   4033
      if(tbn(1,m).eq.0.0) go to 35                                       4034
      nv=0                                                               4035
      ip=m                                                               4036
   37 if(ip.le.0) go to 38                                               4037
      nv=nv+1                                                            4038
      sc(2*nv-1)=tbn(2,ip)                                               4039
      sc(2*nv)=tbn(3,ip)                                                 4040
      ip=tbn(4,ip)+.1                                                    4041
      go to 37                                                           4042
   38 k=match(nv,sc,m-1,tbn,cm,r,1)                                      4043
      if(k.eq.0) go to 35                                                4044
      tbn(1,k)=tbn(1,k)+tbn(1,m)                                         4045
      td(1,k)=td(1,k)+td(1,m)                                            4046
      tbn(1,m)=0.0                                                       4047
      go to 35                                                           4048
   39 return                                                             4049
      end                                                                4050
      subroutine qslice (p,nk,tb,cm,td,kp,kv,lp,lv,tc,r,sc,js)           4051
      integer p,kp(5,*),kv(2,*),lp(3,*),lv(*),js(*)                      4052
      real tb(5,nk),cm(*),td(2,*),tc(*),r(p,2),sc(2,p)                   4053
      do 1 j=1,p                                                         4054
      sc(1,j)=r(j,2)                                                     4055
      sc(2,j)=sc(1,j)+r(j,1)                                             4056
    1 continue                                                           4057
      ll=1                                                               4058
      la=ll                                                              4059
      l1=la                                                              4060
    2 if(kp(1,ll).lt.0) go to 5                                          4061
      if(kp(3,ll) .gt. 0) go to 3                                        4062
      kp(5,ll)=0                                                         4063
      ll=ll+1                                                            4064
      go to 2                                                            4065
    3 kp3=kp(3,ll)                                                       4066
      kp(5,ll)=la                                                        4067
      do 4 m=1,kp3                                                       4068
      l=lp(1,l1)                                                         4069
      nt=lp(3,l1)                                                        4070
      call knts(l,nt,lv(lp(2,l1)),kp(1,ll),kv(1,kp(2,ll)),nk,tb,cm,tc(la 4071
     1),js)                                                              4072
      call side(l,nt,lv(lp(2,l1)),sc,tc(la))                             4073
      l1=l1+1                                                            4074
      la=la+nt*(5*l+1)                                                   4075
    4 continue                                                           4076
      ll=ll+1                                                            4077
      go to 2                                                            4078
    5 le=la-1                                                            4079
      ll=1                                                               4080
      la=ll                                                              4081
      l1=la                                                              4082
      laa=0                                                              4083
    6 if(kp(1,ll).lt.0) go to 13                                         4084
      nv=0                                                               4085
      if(kp(1,ll) .le. 0) go to 8                                        4086
      jl=kp(1,ll)                                                        4087
      do 7 il=1,jl                                                       4088
      k=kp(2,ll)+il-1                                                    4089
      nv=nv+1                                                            4090
      sc(1,nv)=kv(1,k)                                                   4091
      sc(2,nv)=kv(2,k)                                                   4092
    7 continue                                                           4093
      go to 9                                                            4094
    8 if(kp(3,ll) .gt. 0) go to 9                                        4095
      ll=ll+1                                                            4096
      go to 6                                                            4097
    9 if(kp(3,ll) .gt. 0) go to 10                                       4098
      m=match(nv,sc,nk,tb,cm,r,0)                                        4099
      le=le+1                                                            4100
      kp(3,ll)=-le                                                       4101
      tc(le)=td(1,m)                                                     4102
      ll=ll+1                                                            4103
      go to 6                                                            4104
   10 kp3=kp(3,ll)                                                       4105
      do 12 k=1,kp3                                                      4106
      l=lp(1,l1)                                                         4107
      nt=lp(3,l1)                                                        4108
      laa=laa+5*l*nt                                                     4109
      do 11 jp=1,nt                                                      4110
      call gtrm(2,jp,l,nt,lv(lp(2,l1)),dum,dum,nk,tb,tc(la),sc(1,nv+1),d 4111
     1um)                                                                4112
      m=match(nv+l,sc,nk,tb,cm,r,0)                                      4113
      tc(jp+laa)=td(1,m)                                                 4114
   11 continue                                                           4115
      laa=laa+nt                                                         4116
      l1=l1+1                                                            4117
      la=la+nt*(5*l+1)                                                   4118
   12 continue                                                           4119
      ll=ll+1                                                            4120
      go to 6                                                            4121
   13 return                                                             4122
      end                                                                4123
      function ieq(a,b,r)                                                4124
      ieq=0                                                              4125
      if(abs((a-b)/r).lt.1.e-5) ieq=1                                    4126
      return                                                             4127
      end                                                                4128
      function lcm (p,nk,tb,cm)                                          4129
      integer p                                                          4130
      real tb(5,nk),cm(*)                                                4131
      ix=0                                                               4132
      do 1 m=1,nk                                                        4133
      j=abs(tb(2,m))+.1                                                  4134
      if(cm(2*j).eq.0.0) go to 1                                         4135
      if(int(tb(3,m)+.1) .le. ix) go to 1                                4136
      ix=tb(3,m)+.1                                                      4137
      jj=j                                                               4138
    1 continue                                                           4139
      if(ix .le. 0) go to 2                                              4140
      lcm=ix+int(cm(2*jj+1)+.1)-int(cm(2*jj)+.1)+1                       4141
      return                                                             4142
    2 lcm=2*p+1                                                          4143
      do 3 j=1,p                                                         4144
      if(cm(2*j).eq.0.0) go to 3                                         4145
      lcm=lcm+int(cm(2*j+1)+.1)-int(cm(2*j)+.1)+1                        4146
    3 continue                                                           4147
      return                                                             4148
      end                                                                4149
      function newb (m,tb)                                               4150
      real tb(5,m)                                                       4151
      newb=0                                                             4152
      mm1=m-1                                                            4153
      do 1 k=1,mm1                                                       4154
      if(ieq(tb(2,k),tb(2,m),1.0).eq.0) go to 1                          4155
      if(ieq(tb(3,k),tb(3,m),1.0).eq.0) go to 1                          4156
      if(ieq(tb(4,k),tb(4,m),1.0).eq.0) go to 1                          4157
      newb=1                                                             4158
      go to 2                                                            4159
    1 continue                                                           4160
    2 return                                                             4161
      end                                                                4162
      subroutine sscp (n,m,sc,y,w,yb,yv,sw,d,da)                         4163
      real sc(n,*),y(n),w(n)                                             4164
      double precision d(m,m),da(*),yb,yv,sw,s                           4165
      mm1=m-1                                                            4166
      do 6 k=1,mm1                                                       4167
      s=0.d0                                                             4168
      do 1 i=1,n                                                         4169
      s=s+w(i)*sc(i,k)                                                   4170
    1 continue                                                           4171
      s=s/sw                                                             4172
      da(k)=s                                                            4173
      do 2 i=1,n                                                         4174
      sc(i,k)=sc(i,k)-s                                                  4175
    2 continue                                                           4176
      do 4 j=1,k                                                         4177
      s=0.d0                                                             4178
      do 3 i=1,n                                                         4179
      s=s+w(i)*sc(i,j)*sc(i,k)                                           4180
    3 continue                                                           4181
      d(j,k)=s                                                           4182
    4 continue                                                           4183
      s=0.d0                                                             4184
      do 5 i=1,n                                                         4185
      s=s+w(i)*sc(i,k)*(y(i)-yb)                                         4186
    5 continue                                                           4187
      d(k,m)=s                                                           4188
    6 continue                                                           4189
      d(m,m)=sw*yv                                                       4190
      return                                                             4191
      end                                                                4192
      subroutine lsf1 (d,m,xb,yb,al,rss,a,a0,dp)                         4193
      double precision d(m,m),xb(*),yb,al,rss,a(*),a0,dp(*),eps,s        4194
      data eps /1.d-4/                                                   4195
      mm1=m-1                                                            4196
      do 1 i=1,mm1                                                       4197
      dp(i)=d(i,i)                                                       4198
      d(i,i)=d(i,i)*(1.d0+al)                                            4199
    1 continue                                                           4200
      do 5 i=1,mm1                                                       4201
      if(dp(i).le.0.d0) go to 5                                          4202
      im1=i-1                                                            4203
      s=dp(i)                                                            4204
      j=1                                                                4205
      go to 3                                                            4206
    2 j=j+1                                                              4207
    3 if((j).gt.(im1)) go to 4                                           4208
      if(d(j,j).lt.0.d0) s=s+dp(j)*d(j,i)**2                             4209
      go to 2                                                            4210
    4 if((d(i,i)-al*s)/dp(i).lt.eps) go to 5                             4211
      call sweep(d,m,i,-1.d0,dp(m))                                      4212
    5 continue                                                           4213
      rss=0.d0                                                           4214
      a0=yb                                                              4215
      do 6 i=1,mm1                                                       4216
      a(i)=0.d0                                                          4217
      if(d(i,i).ge.0.d0) go to 6                                         4218
      a(i)=d(i,m)                                                        4219
      a0=a0-a(i)*xb(i)                                                   4220
      rss=rss+dp(i)*a(i)**2                                              4221
    6 continue                                                           4222
      rss=d(m,m)-al*rss                                                  4223
      return                                                             4224
      end                                                                4225
      subroutine bkstp (d,m,xb,yb,al,rss,a,a0,k,dp)                      4226
      double precision d(m,m),xb(*),yb,al                                4227
      double precision a(*),a0,rss,dp(*),s                               4228
      data big /9.9e30/                                                  4229
      mm1=m-1                                                            4230
      rss=big                                                            4231
      k=0                                                                4232
      do 4 i=1,mm1                                                       4233
      if(d(i,i).ge.0.d0) go to 4                                         4234
      s=0.d0                                                             4235
      do 3 j=1,mm1                                                       4236
      if(d(j,j).ge.0.d0) go to 3                                         4237
      if(j.eq.i) go to 3                                                 4238
      if(j .ge. i) go to 1                                               4239
      a0=d(j,i)                                                          4240
      go to 2                                                            4241
    1 a0=d(i,j)                                                          4242
    2 s=s+dp(j)*(d(j,m)-a0*d(i,m)/d(i,i))**2                             4243
    3 continue                                                           4244
      s=d(m,m)-d(i,m)**2/d(i,i)-al*s                                     4245
      if(s .gt. rss) go to 4                                             4246
      rss=s                                                              4247
      k=i                                                                4248
    4 continue                                                           4249
      if(k.gt.0) call sweep(d,m,k,1.d0,dp(m))                            4250
      a0=yb                                                              4251
      rss=0.d0                                                           4252
      do 5 i=1,mm1                                                       4253
      a(i)=0.d0                                                          4254
      if(d(i,i).ge.0.d0) go to 5                                         4255
      a(i)=d(i,m)                                                        4256
      a0=a0-a(i)*xb(i)                                                   4257
      rss=rss+dp(i)*a(i)**2                                              4258
    5 continue                                                           4259
      rss=d(m,m)-al*rss                                                  4260
      return                                                             4261
      end                                                                4262
      subroutine sweep (a,m,k,fl,u)                                      4263
      double precision a(m,m),u(m),fl,c                                  4264
      c=a(k,k)                                                           4265
      do 1 i=1,k                                                         4266
      u(i)=a(i,k)                                                        4267
      a(i,k)=0.d0                                                        4268
    1 continue                                                           4269
      do 2 i=k,m                                                         4270
      u(i)=a(k,i)                                                        4271
      a(k,i)=0.d0                                                        4272
    2 continue                                                           4273
      u(k)=fl                                                            4274
      do 4 i=1,m                                                         4275
      do 3 j=i,m                                                         4276
      a(i,j)=a(i,j)-u(i)*u(j)/c                                          4277
    3 continue                                                           4278
    4 continue                                                           4279
      return                                                             4280
      end                                                                4281
      subroutine array(p,n,i,j)                                          4282
      integer p                                                          4283
      i=mod(p,n)                                                         4284
      if(i.eq.0) i=n                                                     4285
      j=(p-i)/n+1                                                        4286
      return                                                             4287
      end                                                                4288
      subroutine cvmars  (ix,n,p,x,y,w,nk,ms,df,fv,mi,lx,it,xm,xs,tb,cm, 4289
     1sc,db,d,mm,wt,cv)                                                  4290
      integer p,mm(n,*),lx(p)                                            4291
      real x(n,p),y(n),w(n),xm(p),xs(p),tb(5,nk),cm(*),sc(*),wt(n,2),cv( 4292
     1nk,4)                                                              4293
      double precision db(n,*),d(nk,*)                                   4294
      data eps,big,dfs,cvm,im /1.e-6,9.9e30,2*0.0,0/                     4295
c     if(it.gt.0) write(it,'(/,'' sample reuse to estimate df:'')')      4296
      if(ix .le. 0) go to 1                                              4297
      nr=ix                                                              4298
      nd=nr                                                              4299
c     if(it.gt.0) write(it,'('' '',i3,'' - fold cross-validation.'',/)') 4300
c    1 ix                                                                4301
      go to 2                                                            4302
    1 nr=1                                                               4303
      nd=-ix                                                             4304
c     if(it.gt.0) write(it,'('' independent test set - every'',i4,'' obs 4305
c    1ervations.'',/)') nd                                               4306
    2 do 3 i=1,n                                                         4307
      wt(i,1)=w(i)                                                       4308
      wt(i,2)=i                                                          4309
    3 continue                                                           4310
      do 4 i=1,n                                                         4311
      call rnms(r,1)                                                     4312
      k=(n-i+1)*r+i                                                      4313
      t=wt(i,2)                                                          4314
      wt(i,2)=wt(k,2)                                                    4315
      wt(k,2)=t                                                          4316
    4 continue                                                           4317
      do 6 i=1,3                                                         4318
      do 5 m=1,nk                                                        4319
      cv(m,i)=0.0                                                        4320
    5 continue                                                           4321
    6 continue                                                           4322
      cv0=0.0                                                            4323
      sw=cv0                                                             4324
      wn=sw                                                              4325
      yv=wn                                                              4326
      fc=yv                                                              4327
      do 14 ir=1,nr                                                      4328
      i=ir                                                               4329
    7 if(i.gt.n) go to 8                                                 4330
      wt(int(wt(i,2)+.1),1)=0.0                                          4331
      i=i+nd                                                             4332
      go to 7                                                            4333
    8 call marsgo (n,p,x,y,wt,nk,ms,df,fv,mi,lx,0,xm,xs,az,tb,cm,sc,db,d 4334
     1,mm)                                                               4335
      yv1=sc(3)                                                          4336
      yv=yv+yv1                                                          4337
      wn1=sc(2)                                                          4338
      wn=wn+wn1                                                          4339
      fc=fc+sc(1)                                                        4340
      mk=sc((nk+1)**2+4)+.1                                              4341
      i=ir                                                               4342
    9 if(i.gt.n) go to 10                                                4343
      k=wt(i,2)+.1                                                       4344
      wt(k,1)=w(k)                                                       4345
      sw=sw+w(k)                                                         4346
      call cvmod(k,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv(1,3))                   4347
      i=i+nd                                                             4348
      go to 9                                                            4349
   10 do 13 m=1,nk                                                       4350
      am=sc(m+4)                                                         4351
      cv(m,2)=cv(m,2)+am                                                 4352
      am1=yv1                                                            4353
      if(m.gt.1) am1=sc(m+3)                                             4354
      if(am1/yv1 .le. eps) go to 11                                      4355
      r=sqrt(am/am1)                                                     4356
      go to 12                                                           4357
   11 r=1.0                                                              4358
   12 cv(m,1)=cv(m,1)+((wn1-1.0)*(1.0-r)/(m-r*(m-1))-1.0)/sc(1)          4359
   13 continue                                                           4360
   14 continue                                                           4361
      do 15 m=1,nk                                                       4362
      cv(m,1)=cv(m,1)/nr                                                 4363
      cv(m,2)=cv(m,2)/nr                                                 4364
      cv(m,3)=cv(m,3)/sw                                                 4365
   15 continue                                                           4366
      fc=fc/nr                                                           4367
      yv=yv/nr                                                           4368
      wn=wn/nr                                                           4369
      cv0=cv0/sw                                                         4370
c     if(it.gt.0) write(it,21)                                           4371
      im=0                                                               4372
      cvm=cv0                                                            4373
      dmx=-big                                                           4374
      cvl=cv(nk,1)                                                       4375
      m=nk                                                               4376
      go to 17                                                           4377
   16 m=m+(-1)                                                           4378
   17 if((-1)*((m)-(1)).gt.0) go to 19                                   4379
      if(cv(m,1).le.dmx) go to 16                                        4380
      dmx=cv(m,1)                                                        4381
      dfu=0.5*(cvl+cv(m,1))                                              4382
      cvl=cv(m,1)                                                        4383
      if(cv(m,3) .gt. cvm) go to 18                                      4384
      cvm=cv(m,3)                                                        4385
      df=dfu                                                             4386
      im=m                                                               4387
   18 gcv=cv(m,2)/(1.0-((dfu*fc+1.0)*m+1.0)/wn)**2                       4388
c     if(it.gt.0) write(it,22) m,dfu,cv(m,2),gcv,cv(m,3)                 4389
      go to 16                                                           4390
   19 if(cv0 .gt. cvm) go to 20                                          4391
      cvm=cv0                                                            4392
      df=dmx                                                             4393
      im=0                                                               4394
   20 dfs=df                                                             4395
      gcv=yv/(1.0-1.0/wn)**2                                             4396
c     if(it.gt.0) write(it,22) 0,dmx,yv,gcv,cv0                          4397
c     if(it.gt.0) write(it,'(/,'' estimated optimal df('',i3,'') ='',f7. 4398
c    12,               '' with (estimated) pse ='',g12.4)') im,df,cvm    4399
      return                                                             4400
      entry cvinfo(a1,a2,ia1)                                            4401
      a1=dfs                                                             4402
      a2=cvm                                                             4403
      ia1=im                                                             4404
      return                                                             4405
   21 format('  #bsfns     df        asr           gcv           cv')    4406
   22 format(' ',i5,f10.2,3g14.4)                                        4407
      end                                                                4408
      subroutine cvmod (i,n,x,y,w,nk,mk,tb,cm,sc,cv0,cv)                 4409
      real x(n,*),y(n),w(n),tb(5,nk),cm(*),sc(*),cv(nk,2)                4410
      do 8 m=1,mk                                                        4411
      t=tb(2,m)                                                          4412
      j=abs(t)+.1                                                        4413
      if(cm(2*j) .le. 0.0) go to 5                                       4414
      k=x(i,j)+.1                                                        4415
      if(k .ne. 0) go to 1                                               4416
      u=0.0                                                              4417
      go to 2                                                            4418
    1 u=cm(k+int(tb(3,m)+.1))                                            4419
    2 if(t .ge. 0.0) go to 6                                             4420
      if(u .ne. 0.0) go to 3                                             4421
      u=1.0                                                              4422
      go to 6                                                            4423
    3 u=0.0                                                              4424
      go to 6                                                            4425
    5 u=amax1(0.0,sign(1.0,t)*(x(i,j)-tb(3,m)))                          4426
    6 l=tb(4,m)+.1                                                       4427
      if(l .le. 0) go to 7                                               4428
      cv(m,2)=u*cv(l,2)                                                  4429
      go to 8                                                            4430
    7 cv(m,2)=u                                                          4431
    8 continue                                                           4432
      kp=nk+4                                                            4433
      cv0=cv0+w(i)*(y(i)-sc(4))**2                                       4434
      do 10 m=1,nk                                                       4435
      kp=kp+1                                                            4436
      s=sc(kp)                                                           4437
      do 9 l=1,nk                                                        4438
      kp=kp+1                                                            4439
      if(l.le.mk) s=s+sc(kp)*cv(l,2)                                     4440
    9 continue                                                           4441
      cv(m,1)=cv(m,1)+w(i)*(y(i)-s)**2                                   4442
   10 continue                                                           4443
      return                                                             4444
      end                                                                4445
      subroutine nest (n,i,j,nv,vals)                                    4446
      parameter (mlist=200, nlist=2000)                                  4447
      real vals(*),vm(nlist),tb(5,*),cm(*),x(n,*),bl(*)                  4448
      integer m(4,mlist),p,lx(*)                                         4449
      save m,vm                                                          4450
      data il,jl /2*0/                                                   4451
      if((i .ne. 0) .and. (j .ne. 0)) go to 1                            4452
      il=0                                                               4453
      jl=il                                                              4454
      return                                                             4455
    1 if(i.eq.j) return                                                  4456
      ig=0                                                               4457
      if(nv .le. 0) go to 8                                              4458
      k=1                                                                4459
      go to 3                                                            4460
    2 k=k+1                                                              4461
    3 if((k).gt.(il)) go to 4                                            4462
      if(m(1,k).eq.i.or.m(1,k).eq.j) return                              4463
      go to 2                                                            4464
    4 il=il+1                                                            4465
      if(il .le. mlist) go to 5                                          4466
c     write(6,  '('' increase parameter mlist in subroutine nest to grea 4467
c    1ter than'',               i5,'' and recompile.'')') il             4468
      stop                                                               4469
    5 m(1,il)=i                                                          4470
      m(2,il)=j                                                          4471
      m(3,il)=nv                                                         4472
      m(4,il)=jl                                                         4473
      if(jl+nv .le. nlist) go to 6                                       4474
c     write(6,  '('' increase parameter nlist in subroutine nest to grea 4475
c    1ter than'',               i5,'' and recompile.'')') jl+nv          4476
      stop                                                               4477
    6 do 7 k=1,nv                                                        4478
      jl=jl+1                                                            4479
      vm(jl)=vals(k)                                                     4480
    7 continue                                                           4481
      return                                                             4482
    8 k=1                                                                4483
      go to 10                                                           4484
    9 k=k+1                                                              4485
   10 if((k).gt.(il)) go to 12                                           4486
      if(m(1,k) .ne. i .or. m(2,k) .ne. j) go to 9                       4487
      ig=1                                                               4488
   12 if(ig.eq.0) return                                                 4489
      il=il-1                                                            4490
      ll=k                                                               4491
      go to 14                                                           4492
   13 ll=ll+1                                                            4493
   14 if((ll).gt.(il)) go to 16                                          4494
      do 15 l=1,4                                                        4495
      m(l,ll)=m(l,ll+1)                                                  4496
   15 continue                                                           4497
      go to 13                                                           4498
   16 return                                                             4499
      entry nstlst(it)                                                   4500
      if(it.le.0) return                                                 4501
      if(il.eq.0) return                                                 4502
c     write(it,'(/,'' variable nesting:'',/)')                           4503
      do 18 k=1,il                                                       4504
      if(m(3,k) .le. 5) go to 17                                         4505
c     write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') ='' 4506
c    1)')  k,m(1,k),m(2,k)                                               4507
c     write(it,'(100('' '',10f7.1))') (vm(l),l=m(4,k)+1,m(4,k)+m(3,k))   4508
      go to 18                                                           4509
   17 continue
c     write(it,'('' '',i3,'': var('',i3,'') exists for var('',i3,'') ='' 4510
c    1,5f7.1)')  k,m(1,k),m(2,k),(vm(l),l=m(4,k)+1,m(4,k)+m(3,k))        4511
   18 continue                                                           4512
      return                                                             4513
      entry oknest(it,p,lx,cm)                                           4514
      if(it.le.0) return                                                 4515
      l=1                                                                4516
      go to 20                                                           4517
   19 l=l+1                                                              4518
   20 if((l).gt.(il)) go to 24                                           4519
      j1=m(1,l)                                                          4520
      jn=m(2,l)                                                          4521
      jv=m(3,l)                                                          4522
      jp=m(4,l)                                                          4523
c     if(j1.lt.1.or.j1.gt.p) write(it,25) l,j1                           4524
c     if(jn.lt.1.or.jn.gt.p) write(it,26) l,jn                           4525
c     if(lx(jn).ge.0) write(it,27) l,jn,lx(jn)                           4526
      k1=cm(2*jn)+.1                                                     4527
      k2=cm(2*jn+1)+.1                                                   4528
      do 23 k=jp+1,jp+jv                                                 4529
      ig=0                                                               4530
      do 21 kk=k1,k2                                                     4531
      if(vm(k) .ne. cm(kk)) go to 21                                     4532
      ig=1                                                               4533
      go to 22                                                           4534
   21 continue                                                           4535
   22 continue
c     if(ig.eq.0) write(it,28) l,vm(k),jn                                4536
   23 continue                                                           4537
      go to 19                                                           4538
   24 return                                                             4539
   25 format(' nesting entry',i3,', invalid variable',i3,' to be nested. 4540
     1')                                                                 4541
   26 format(' nesting entry',i3,', invalid nesting variable',i3,'.')    4542
   27 format(' nesting entry',i3,', lx(',i3,') =',i2,'. must be < 0.')   4543
   28 format(' nesting entry',i3,', categorical value ',g12.4,/,  ' not  4544
     1among the data values for variable',i3,'.')                        4545
      entry isnstr(j,jb)                                                 4546
      jb=0                                                               4547
      k=1                                                                4548
      go to 30                                                           4549
   29 k=k+1                                                              4550
   30 if((k).gt.(il)) go to 32                                           4551
      if(m(2,k) .ne. j) go to 29                                         4552
      jb=m(1,k)                                                          4553
   32 return                                                             4554
      entry isfac (lm,j,mk,tb,cm,ja)                                     4555
      ja=0                                                               4556
      ig=ja                                                              4557
      l=1                                                                4558
      go to 34                                                           4559
   33 l=l+1                                                              4560
   34 if((l).gt.(il)) go to 36                                           4561
      if(j .ne. m(1,l)) go to 33                                         4562
      ig=1                                                               4563
   36 if(ig.eq.0) return                                                 4564
      jn=m(2,l)                                                          4565
      if(cm(2*jn).eq.0.0) return                                         4566
      jv=m(3,l)                                                          4567
      jp=m(4,l)                                                          4568
      ig=0                                                               4569
      ip=lm                                                              4570
   37 if(ip.le.0) go to 39                                               4571
      j1=abs(tb(2,ip))+.1                                                4572
      if(j1 .ne. jn) go to 38                                            4573
      ig=1                                                               4574
      go to 39                                                           4575
   38 ip=tb(4,ip)+.1                                                     4576
      go to 37                                                           4577
   39 if(ig .eq. 0) go to 45                                             4578
      nc=cm(2*jn+1)-cm(2*jn)+1.1                                         4579
      t=tb(2,ip)                                                         4580
      kp=tb(3,ip)+.1                                                     4581
      do 44 l=1,nc                                                       4582
      lp=l+kp                                                            4583
      if(t .le. 0) go to 40                                              4584
      if(cm(lp).eq.0.0) go to 44                                         4585
      go to 41                                                           4586
   40 if(cm(lp).ne.0.0) go to 44                                         4587
   41 ex=cm(int(cm(2*jn)+.1)+l-1)                                        4588
      ig=0                                                               4589
      do 42 k=jp+1,jp+jv                                                 4590
      if(ex .ne. vm(k)) go to 42                                         4591
      ig=1                                                               4592
      go to 43                                                           4593
   42 continue                                                           4594
   43 if(ig .ne. 0) go to 44                                             4595
      ja=-1                                                              4596
      return                                                             4597
   44 continue                                                           4598
      return                                                             4599
   45 ja=l                                                               4600
      norm=nord(lm,tb)+1                                                 4601
      nc=cm(2*jn+1)-cm(2*jn)+1.1                                         4602
      do 56 lk=1,mk                                                      4603
      if(nord(lk,tb).ne.norm) go to 56                                   4604
      jg=0                                                               4605
      ip=lk                                                              4606
   46 if(ip.le.0) go to 55                                               4607
      t1=tb(2,ip)                                                        4608
      j1=abs(t1)+.1                                                      4609
      if(j1 .ne. jn) go to 54                                            4610
      kp=tb(3,ip)+.1                                                     4611
      kg=0                                                               4612
      do 52 l=1,nc                                                       4613
      lp=l+kp                                                            4614
      lon=cm(lp)+.1                                                      4615
      if(t1 .ge. 0.0) go to 48                                           4616
      if(lon .ne. 0) go to 47                                            4617
      lon=1                                                              4618
      go to 48                                                           4619
   47 lon=0                                                              4620
   48 ex=cm(int(cm(2*jn)+.1)+l-1)                                        4621
      ig=0                                                               4622
      do 49 k=jp+1,jp+jv                                                 4623
      if(ex .ne. vm(k)) go to 49                                         4624
      ig=1                                                               4625
      go to 50                                                           4626
   49 continue                                                           4627
   50 if(lon .ne. 1 .or. ig .ne. 0) go to 51                             4628
      kg=1                                                               4629
      go to 53                                                           4630
   51 if(lon .ne. 0 .or. ig .ne. 1) go to 52                             4631
      kg=1                                                               4632
      go to 53                                                           4633
   52 continue                                                           4634
   53 if(kg .ne. 0) go to 54                                             4635
      jg=1                                                               4636
      go to 55                                                           4637
   54 ip=tb(4,ip)+.1                                                     4638
      go to 46                                                           4639
   55 if(jg.eq.0) go to 56                                               4640
      if(ieqbf(lk,lm,tb,cm) .ne. 1) go to 56                             4641
      ja=-1                                                              4642
      return                                                             4643
   56 continue                                                           4644
      return                                                             4645
      entry cmpnst(ja,n,x,cm,bl)                                         4646
      jn=m(2,ja)                                                         4647
      jv=m(3,ja)                                                         4648
      jp=m(4,ja)                                                         4649
      do 59 l=1,n                                                        4650
      kx=x(l,jn)+.1                                                      4651
      ex=cm(int(cm(2*jn)+.1)+kx-1)                                       4652
      ig=0                                                               4653
      do 57 k=jp+1,jp+jv                                                 4654
      if(ex .ne. vm(k)) go to 57                                         4655
      ig=1                                                               4656
      go to 58                                                           4657
   57 continue                                                           4658
   58 if(ig.eq.1) go to 59                                               4659
      bl(l)=0.0                                                          4660
   59 continue                                                           4661
      return                                                             4662
      entry getnst(ja,cm,j,nv,vals)                                      4663
      j=m(2,ja)                                                          4664
      jv=m(3,ja)                                                         4665
      jp=m(4,ja)                                                         4666
      nv=cm(2*j+1)-cm(2*j)+1.1                                           4667
      do 60 k=1,nv                                                       4668
      vals(k)=0.0                                                        4669
   60 continue                                                           4670
      do 61 l=jp+1,jp+jv                                                 4671
      k=icat(vm(l),j,cm)                                                 4672
      if(k.gt.0) vals(k)=1.0                                             4673
   61 continue                                                           4674
      return                                                             4675
      end                                                                4676
      subroutine blf0(l,ja,n,x,w,cm,sc,nnt,bl)                           4677
      real x(n,*),w(n),cm(*),sc(n,*),bl(n)                               4678
      nnt=0                                                              4679
      call blf(l,n,sc,bl)                                                4680
      if(ja.gt.0) call cmpnst(ja,n,x,cm,bl)                              4681
      do 1 i=1,n                                                         4682
      if(bl(i).gt.0.0.and.w(i).gt.0.0) nnt=nnt+1                         4683
    1 continue                                                           4684
      return                                                             4685
      end                                                                4686
      subroutine blf(l,n,sc,bl)                                          4687
      real sc(n,*),bl(n)                                                 4688
      if(l .gt. 0) go to 2                                               4689
      do 1 i=1,n                                                         4690
      bl(i)=1.0                                                          4691
    1 continue                                                           4692
      go to 4                                                            4693
    2 do 3 i=1,n                                                         4694
      bl(i)=sc(i,l)                                                      4695
    3 continue                                                           4696
    4 return                                                             4697
      end                                                                4698
      subroutine mnspan(ms,alf,nep,nnt,mn,me,mel)                        4699
      parameter(al2=0.693147,al25=1.732868)                              4700
      allf=-alog(1.0-alf)                                                4701
      fmn=-alog(allf/(nep*nnt))/al25                                     4702
      fme=-alog(alf*0.125/nep)/al2                                       4703
      if(ms .le. 0) go to 1                                              4704
      me=ms*fme/fmn+0.5                                                  4705
      mn=ms                                                              4706
      go to 2                                                            4707
    1 me=fme+0.5                                                         4708
      mn=fmn+0.5                                                         4709
    2 me=max0(me,mn,2)                                                   4710
      nst=nnt-2*me-1                                                     4711
      nnr=nst/mn                                                         4712
      nnl=nst-nnr*mn                                                     4713
      nnr=(nnr+1)*mn-nst                                                 4714
      nst=min0(nnl,nnr)                                                  4715
      if(nnl .gt. nnr) go to 3                                           4716
      nnl=1                                                              4717
      go to 4                                                            4718
    3 nnl=-1                                                             4719
    4 nnr=nst/2                                                          4720
      mel=me                                                             4721
      me=me+nnl*nnr                                                      4722
      mel=mel+nnl*nnr                                                    4723
      if(mod(nst,2).ne.0) mel=mel+nnl                                    4724
      return                                                             4725
      end                                                                4726
      function ieqbf(lk,lm,tb,cm)                                        4727
      real tb(5,*),cm(*)                                                 4728
      ipo=lm                                                             4729
      lg=0                                                               4730
    1 if(ipo.le.0) go to 16                                              4731
      to=tb(2,ipo)                                                       4732
      jo=abs(to)+.1                                                      4733
      jg=0                                                               4734
      if(cm(2*jo) .ne. 0.0) go to 2                                      4735
      t=tb(3,ipo)                                                        4736
      ic=0                                                               4737
      go to 3                                                            4738
    2 ko=tb(3,ipo)+.1                                                    4739
      nc=cm(2*jo+1)-cm(2*jo)+1.1                                         4740
      ic=1                                                               4741
    3 ip=lk                                                              4742
    4 if(ip.le.0) go to 14                                               4743
      t1=tb(2,ip)                                                        4744
      j1=abs(t1)+.1                                                      4745
      if(j1 .ne. jo) go to 13                                            4746
      if(ic .ne. 0) go to 6                                              4747
      if(to*t1 .le. 0.0) go to 13                                        4748
      if(ieq(t,tb(3,ip),1.0) .ne. 1) go to 13                            4749
      jg=1                                                               4750
      go to 14                                                           4751
    6 kp=tb(3,ip)+.1                                                     4752
      kg=0                                                               4753
      do 11 l=1,nc                                                       4754
      lo=l+ko                                                            4755
      lp=l+kp                                                            4756
      lon=cm(lo)+.1                                                      4757
      lop=cm(lp)+.1                                                      4758
      if(to .ge. 0.0) go to 8                                            4759
      if(lon .ne. 0) go to 7                                             4760
      lon=1                                                              4761
      go to 8                                                            4762
    7 lon=0                                                              4763
    8 if(t1 .ge. 0.0) go to 10                                           4764
      if(lop .ne. 0) go to 9                                             4765
      lop=1                                                              4766
      go to 10                                                           4767
    9 lop=0                                                              4768
   10 if(lon .eq. lop) go to 11                                          4769
      kg=1                                                               4770
      go to 12                                                           4771
   11 continue                                                           4772
   12 if(kg .ne. 0) go to 13                                             4773
      jg=1                                                               4774
      go to 14                                                           4775
   13 ip=tb(4,ip)+.1                                                     4776
      go to 4                                                            4777
   14 if(jg .ne. 0) go to 15                                             4778
      lg=1                                                               4779
      go to 16                                                           4780
   15 ipo=tb(4,ipo)+.1                                                   4781
      go to 1                                                            4782
   16 if(lg .ne. 0) go to 17                                             4783
      ieqbf=1                                                            4784
      go to 18                                                           4785
   17 ieqbf=0                                                            4786
   18 return                                                             4787
      end                                                                4788
      function ibfext(m,tb,cm)                                           4789
      real tb(5,*),cm(*)                                                 4790
      mm1=m-1                                                            4791
      ibfext=0                                                           4792
      norm=nord(m,tb)                                                    4793
      do 1 l=1,mm1                                                       4794
      if(nord(l,tb).ne.norm) go to 1                                     4795
      if(ieqbf(l,m,tb,cm) .eq. 0) go to 1                                4796
      ibfext=1                                                           4797
      return                                                             4798
    1 continue                                                           4799
      return                                                             4800
      end                                                                4801
      subroutine miss (n,p,x,lx,xm,flg,pn,xn,lxn,xs,xp)                  4802
      integer p,pn,lx(*),lxn(*)                                          4803
      real x(n,*),xm(*),xn(n,*),xs(*),xp(*)                              4804
      double precision s                                                 4805
      pn=p                                                               4806
      xp(1)=p                                                            4807
      do 1 j=2,2*p+1                                                     4808
      xp(j)=0.0                                                          4809
    1 continue                                                           4810
      ss=0.0                                                             4811
      do 7 j=1,p                                                         4812
      lxn(j)=lx(j)                                                       4813
      xs(j)=flg                                                          4814
      if(lx(j).eq.0) go to 7                                             4815
      s=0.d0                                                             4816
      mf=0                                                               4817
      do 4 i=1,n                                                         4818
      xn(i,j)=x(i,j)                                                     4819
      if(x(i,j) .ne. xm(j)) go to 2                                      4820
      mf=mf+1                                                            4821
      go to 4                                                            4822
    2 if(lx(j) .ge. 0) go to 3                                           4823
      ss=x(i,j)                                                          4824
      go to 4                                                            4825
    3 s=s+x(i,j)                                                         4826
    4 continue                                                           4827
      if(mf.eq.0) go to 7                                                4828
      if(mf .ne. n) go to 5                                              4829
      lxn(j)=0                                                           4830
      go to 7                                                            4831
    5 s=s/(n-mf)                                                         4832
      pn=pn+1                                                            4833
      lxn(pn)=-1                                                         4834
      xs(pn)=1.0                                                         4835
      xp(j+1)=pn                                                         4836
      call nest(n,j,pn,1,1.0)                                            4837
      if(lx(j).gt.0) ss=s                                                4838
      xp(j+p+1)=ss                                                       4839
      do 6 i=1,n                                                         4840
      xn(i,pn)=1.0                                                       4841
      if(x(i,j) .ne. xm(j)) go to 6                                      4842
      xn(i,j)=ss                                                         4843
      xn(i,pn)=0.0                                                       4844
    6 continue                                                           4845
    7 continue                                                           4846
      return                                                             4847
      end                                                                4848
      subroutine mkmiss (n,p,x,y,w,xm,pm,nnx,nn,xnn,yn,wn,sc)            4849
      parameter(nlist=500)                                               4850
      integer p,m(nlist)                                                 4851
      real pm(p),xm(p),x(n,p),y(n),w(n),xnn(*),yn(*),wn(*),sc(p,*)       4852
      data tol /0.001/                                                   4853
      if(p .le. nlist) go to 1                                           4854
c     write(6,'('' increase parameter nlist in subroutine mkmiss to '',i 4855
c    15,                      '' and recompile.'')') p                   4856
      stop                                                               4857
    1 do 3 j=1,p                                                         4858
      m(j)=0                                                             4859
      do 2 i=1,n                                                         4860
      if(x(i,j).eq.xm(j)) m(j)=m(j)+1                                    4861
      sc(j,i)=x(i,j)                                                     4862
      yn(i)=y(i)                                                         4863
      wn(i)=w(i)                                                         4864
    2 continue                                                           4865
    3 continue                                                           4866
      nn=n                                                               4867
      jp=0                                                               4868
      km=jp                                                              4869
    4 if(nn.ge.nnx.or.km.gt.p) go to 13                                  4870
      jp=jp+1                                                            4871
      if(jp.gt.p) jp=1                                                   4872
      fin=nn*pm(jp)-m(jp)                                                4873
      if(fin .le. 0.0) go to 5                                           4874
      in=fin+0.5                                                         4875
      go to 6                                                            4876
    5 in=0                                                               4877
    6 in=min0(in,nnx-nn)                                                 4878
      if(in .le. 0) go to 7                                              4879
      km=0                                                               4880
      go to 8                                                            4881
    7 km=km+1                                                            4882
      go to 4                                                            4883
    8 do 11 k=1,in                                                       4884
      call rnms(r,1)                                                     4885
      i=nn*r+1.0                                                         4886
      nnk=nn+k                                                           4887
      do 9 j=1,p                                                         4888
      sc(j,nnk)=sc(j,i)                                                  4889
    9 continue                                                           4890
      sc(jp,nnk)=xm(jp)                                                  4891
      yn(nnk)=yn(i)                                                      4892
      wn(nnk)=wn(i)                                                      4893
      do 10 j=1,p                                                        4894
      if(sc(j,nnk).eq.xm(j)) m(j)=m(j)+1                                 4895
   10 continue                                                           4896
   11 continue                                                           4897
      nn=nn+in                                                           4898
      cvx=-9.9e30                                                        4899
      do 12 j=1,p                                                        4900
      cvx=amax1(cvx,(nn*pm(j)-m(j))/float(nn))                           4901
   12 continue                                                           4902
      if(cvx.lt.tol) go to 13                                            4903
      go to 4                                                            4904
   13 k=0                                                                4905
      do 15 j=1,p                                                        4906
      do 14 i=1,nn                                                       4907
      k=k+1                                                              4908
      xnn(k)=sc(j,i)                                                     4909
   14 continue                                                           4910
   15 continue                                                           4911
      return                                                             4912
      entry smktol(val)                                                  4913
      tol=val                                                            4914
      return                                                             4915
      end                                                                4916
      subroutine xmiss (n,x,xm,xp,xn)                                    4917
      real x(n,*),xm(*),xp(*),xn(n,*)                                    4918
      integer p                                                          4919
      p=xp(1)+.1                                                         4920
      do 3 j=1,p                                                         4921
      k=xp(j+1)+.1                                                       4922
      do 2 i=1,n                                                         4923
      if(x(i,j) .eq. xm(j)) go to 1                                      4924
      xn(i,j)=x(i,j)                                                     4925
      if(k.gt.0) xn(i,k)=1.0                                             4926
      go to 2                                                            4927
    1 xn(i,j)=xp(j+p+1)                                                  4928
      if(k.gt.0) xn(i,k)=0.0                                             4929
    2 continue                                                           4930
    3 continue                                                           4931
      return                                                             4932
      end                                                                4933
      function nnord(m,tb)                                               4934
      real tb(5,*)                                                       4935
      ip=m                                                               4936
      nnord=0                                                            4937
    1 if(ip.le.0) go to 2                                                4938
      call isnstr(int(abs(tb(2,ip))+.1),jb)                              4939
      if(jb.eq.0) nnord=nnord+1                                          4940
      ip=tb(4,ip)+.1                                                     4941
      go to 1                                                            4942
    2 return                                                             4943
      end                                                                4944
      subroutine spofa(a,m,n,info)                                              
      double precision a(m,*),s,t,u                                             
         j1 = info                                                              
         do 30 j = j1, n                                                        
            info=j                                                              
            s = 0.0d0                                                           
            jm1 = j - 1                                                         
            if (jm1 .lt. 1) go to 20                                            
            do 10 k = 1, jm1                                                    
               u=0.0                                                            
               km1=k-1                                                          
               if(km1.le.0) go to 40                                            
               do 50 i=1,km1                                                    
                  u=u+a(i,k)*a(i,j)                                             
   50          continue                                                         
   40          continue                                                         
               t = a(k,j) - u                                                   
               t = t/a(k,k)                                                     
               a(k,j) = t                                                       
               s = s + t*t                                                      
   10       continue                                                            
   20       continue                                                            
            s = a(j,j) - s                                                      
            if (s .le. 0.0d0)  return                                           
            a(j,j) = dsqrt(s)                                                   
   30    continue                                                               
      info=0                                                                    
      return                                                                    
      end                                                                       
      subroutine sposl(a,m,n,b)                                                 
      double precision a(m,*),b(*),t                                            
      do 10 k = 1, n                                                            
         t = 0.0                                                                
         km1=k-1                                                                
         if(km1.le.0) go to 30                                                  
         do 40 i=1,km1                                                          
            t=t+a(i,k)*b(i)                                                     
   40    continue                                                               
   30    continue                                                               
         b(k) = (b(k) - t)/a(k,k)                                               
   10 continue                                                                  
      do 20 kb=1,n                                                              
      k=n+1-kb                                                                  
      b(k)=b(k)/a(k,k)                                                          
      t=-b(k)                                                                   
      km1=k-1                                                                   
      if(km1.le.0) go to 50                                                     
      if(t.eq.0.0) go to 50                                                     
      do 60 i=1,km1                                                             
         b(i)=b(i)+t*a(i,k)                                                     
   60 continue                                                                  
   50 continue                                                                  
   20 continue                                                                  
      return                                                                    
      end                                                                       
      subroutine psort (v,a,ii,jj)                                              
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      real v                                                                    
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       
      subroutine stseed (iseed)                                                 
      real x(1)                                                                 
      double precision u                                                        
      data i /987654321/                                                        
      i=iseed                                                                   
      return                                                                    
      entry rnms (x,n)                                                          
      do 1 j=1,n                                                                
      i=dmod(i*16807.d0,2147483647.d0)                                          
      u=i                                                                       
      u=u*.465661287d-9                                                         
      x(j)=u                                                                    
    1 continue                                                                  
      return                                                                    
      end                                                                       


