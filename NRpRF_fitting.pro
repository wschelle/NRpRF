pro NRpRF_fitting
;+
; :Author: wouter
;
; Non-rigid pRF fitting
;-

;*****Declare variables*****
subs=['V6649','V6682','V6683','V6743','V6744','V7774','V7860','V7882']
nsubs=n_elements(subs)
funrs=['1','2','3','4']
conditions = ['toes', 'ankle', 'knee', 'abdomen', 'shoulder', 'elbow','wrist', 'littleF', 'ringF', 'middleF', 'indexF', 'thumb', 'forehead', 'eye', 'nose', 'lips', 'jaw', 'tongue']

for ii=0,nsubs-1 do begin
  sub=subs[ii]
  subdir='/Fridge/users/wouter/NRPRF/'+sub+'/FSL/'
  filedir1=subdir+funrs[0]+'/'
  filedir2=subdir+funrs[1]+'/'
  filedir3=subdir+funrs[3]+'/'
  filedir4=subdir+funrs[4]+'/'
  surfdir='/Fridge/users/wouter/fs/'+sub+'/surf/caret/'
  expcode='wbm'
  scandur=2100.
  nconds=18.

  ;Load fMRI timeseries
  readmgh,filedir1+sub+'_funcs_'+funrs[0]+'_avg01_lh.mgh',lfuncs1
  readmgh,filedir2+sub+'_funcs_'+funrs[1]+'_avg01_lh.mgh',lfuncs2
  readmgh,filedir3+sub+'_funcs_'+funrs[2]+'_avg01_lh.mgh',lfuncs3
  readmgh,filedir4+sub+'_funcs_'+funrs[3]+'_avg01_lh.mgh',lfuncs4
  readmgh,filedir1+sub+'_funcs_'+funrs[0]+'_avg01_rh.mgh',rfuncs1
  readmgh,filedir2+sub+'_funcs_'+funrs[1]+'_avg01_rh.mgh',rfuncs2
  readmgh,filedir3+sub+'_funcs_'+funrs[2]+'_avg01_rh.mgh',rfuncs3
  readmgh,filedir4+sub+'_funcs_'+funrs[3]+'_avg01_rh.mgh',rfuncs4
  lfuncs1=reform(lfuncs1)
  lfuncs2=reform(lfuncs2)
  lfuncs3=reform(lfuncs3)
  lfuncs3=reform(lfuncs4)
  rfuncs1=reform(rfuncs1)
  rfuncs2=reform(rfuncs2)
  rfuncs3=reform(rfuncs3)
  rfuncs4=reform(rfuncs4)

  nrnodesl=n_elements(lfuncs1[*,0])
  nrnodesr=n_elements(rfuncs1[*,0])
  nrnodes=nrnodesl+nrnodesr
  nscans1=n_elements(lfuncs1[0,*])
  nscans4=n_elements(lfuncs4[0,*])
  funx1=[lfuncs1,rfuncs1]
  funx2=[lfuncs2,rfuncs2]
  funx4=[lfuncs3,rfuncs3]
  funx5=[lfuncs4,rfuncs4]

  ;Concatenate different runs and clean up memory
  funx=transpose([transpose(funx1),transpose(funx2),transpose(funx4),transpose(funx5)])
  delvar,lfuncs1
  delvar,lfuncs2
  delvar,lfuncs3
  delvar,lfuncs4
  delvar,rfuncs1
  delvar,rfuncs2
  delvar,rfuncs3
  delvar,rfuncs4
  delvar,funx1
  delvar,funx2
  delvar,funx3
  delvar,funx4
  nscans=n_elements(funx[0,*])

  ;Make a binary mask
  mask=intarr(nrnodes)
  mask(where(mean(funx,dimension=2) gt 10))=1

  ;Create high-pass filter using a discrete set of cosines with a cut-off at 0.01Hz
  nrfilters = floor(2*(nscans*(scandur/1000.)*0.01));
  filtermatrix=make_filter(nrfilters,nscans)
  
  ;Add additional regressors for the different runs
  sesfactor=fltarr(1,nscans)
  sesfactor[0,0:nscans1-1]=1
  filtermatrix=[filtermatrix,sesfactor]
  sesfactor*=0
  sesfactor[0,nscans1:2*nscans1-1]=1
  filtermatrix=[filtermatrix,sesfactor]
  sesfactor*=0
  sesfactor[0,2*nscans1:2*nscans1+nscans4-1]=1
  filtermatrix=[filtermatrix,sesfactor]
  delvar,sesfactor
  
  ;High-pass filtering of timeseries using linear regression + scaling to percent signal change
  for i=0L,nrnodes-1 do if mask[i] ne 0 then begin
    b=regress(filtermatrix,reform(funx[i,*]),yfit=yfit,const=const,/double)
    funx[i,*]=(funx[i,*]-yfit)/const*100.
  end

  ;Read in experiment onset times
  ti1=read_ascii('/Fridge/users/wouter/WholeBodyMap/timings/'+sub+'timing_ses1.dat')
  ti1=ti1.(0)
  ti2=read_ascii('/Fridge/users/wouter/WholeBodyMap/timings/'+sub+'timing_ses2.dat')
  ti2=ti2.(0)
  ti3=read_ascii('/Fridge/users/wouter/WholeBodyMap/timings/'+sub+'timing_ses3.dat')
  ti3=ti3.(0)
  ti4=read_ascii('/Fridge/users/wouter/WholeBodyMap/timings/'+sub+'timing_ses4.dat')
  ti4=ti4.(0)

  ;Reorganize all conditions to match cortical homunculus
  ti2[1,where(ti2[1,*] eq 9)]=18
  ti2[1,where(ti2[1,*] eq 5)]=17
  ti2[1,where(ti2[1,*] eq 7)]=16
  ti2[1,where(ti2[1,*] eq 8)]=15
  ti2[1,where(ti2[1,*] eq 3)]=14
  ti2[1,where(ti2[1,*] eq 6)]=13
  ti2[1,where(ti2[1,*] eq 1)]=3
  ti2[1,where(ti2[1,*] eq 2)]=2
  ti2[1,where(ti2[1,*] eq 4)]=1
  ti1[1,where(ti1[1,*] eq 7)]=12
  ti1[1,where(ti1[1,*] eq 2)]=11
  ti1[1,where(ti1[1,*] eq 4)]=10
  ti1[1,where(ti1[1,*] eq 8)]=4
  ti1[1,where(ti1[1,*] eq 9)]=7
  ti1[1,where(ti1[1,*] eq 5)]=9
  ti1[1,where(ti1[1,*] eq 6)]=5
  ti1[1,where(ti1[1,*] eq 1)]=6
  ti1[1,where(ti1[1,*] eq 3)]=8
  ti4[1,where(ti4[1,*] eq 9)]=18
  ti4[1,where(ti4[1,*] eq 5)]=17
  ti4[1,where(ti4[1,*] eq 7)]=16
  ti4[1,where(ti4[1,*] eq 8)]=15
  ti4[1,where(ti4[1,*] eq 3)]=14
  ti4[1,where(ti4[1,*] eq 6)]=13
  ti4[1,where(ti4[1,*] eq 1)]=3
  ti4[1,where(ti4[1,*] eq 2)]=2
  ti4[1,where(ti4[1,*] eq 4)]=1
  ti3[1,where(ti3[1,*] eq 7)]=12
  ti3[1,where(ti3[1,*] eq 2)]=11
  ti3[1,where(ti3[1,*] eq 4)]=10
  ti3[1,where(ti3[1,*] eq 8)]=4
  ti3[1,where(ti3[1,*] eq 9)]=7
  ti3[1,where(ti3[1,*] eq 5)]=9
  ti3[1,where(ti3[1,*] eq 6)]=5
  ti3[1,where(ti3[1,*] eq 1)]=6
  ti3[1,where(ti3[1,*] eq 3)]=8

  ti2[0,*]+=nscans1*scandur
  ti3[0,*]+=(2*nscans1*scandur)
  ti4[0,*]+=(2*nscans1*scandur)+(nscans4*scandur)

  ;Make a binary design matrix for each condition
  fmat=fltarr(nconds,nscans*scandur)
  for i=0,44 do fmat[ti1[1,i]-1,ti1[0,i]]=1
  for i=0,44 do fmat[ti1[1,i]-1,ti1[0,i]+1000]=1
  for i=0,44 do fmat[ti2[1,i]-1,ti2[0,i]]=1
  for i=0,44 do fmat[ti2[1,i]-1,ti2[0,i]+1000]=1
  for i=0,35 do fmat[ti4[1,i]-1,ti4[0,i]]=1
  for i=0,35 do fmat[ti4[1,i]-1,ti4[0,i]+1000]=1
  for i=0,35 do fmat[ti5[1,i]-1,ti5[0,i]]=1
  for i=0,35 do fmat[ti5[1,i]-1,ti5[0,i]+1000]=1

  ;Convolve binary design matrix with canonical HRF and interpolate across temporal domain from ms to scans.
  cfmat=bconv(fmat)
  fmatrix=fltarr(nconds,nscans)
  for i=0,nconds-1 do fmatrix[i,*]=interpol(cfmat[i,*],nscans,/spline)

  ;The Non-rigid pRF fitting is performed by the following 2 stages:
  ;First, the postions of conditions in a static Gaussian response field are estimated using a multiple linear regression
  ;Second, non-rigid pRF timeseries are fitted using the Levenberg-Marquardt LS algorithm with the first stage as starting point.
  
  ;Perform a multiple linear regression on fmri timeseries with convolved design matrix
  regord=fltarr(nrnodes,nconds)
  cmrsq=fltarr(nrnodes)
  for i=0L,nrnodes-1 do if mask[i] eq 1 then begin
    regord[i,*]=regress(fmatrix,reform(funx[i,*]),mcorrelation=mcor)
    cmrsq[i]=mcor^2.
  end
  
  ;normalize the regression coefficients
  nreg=regord
  for i=0L,nrnodes-1 do if mask[i] eq 1 then nreg[i,*]/=max(regord[i,*])

  ;Make a static Gaussian Response field
  step=100.
  gcenter=10.
  gsigma=1.
  xax=findgen(gcenter*step)/step
  p0=[0.,1.,gsigma,gcenter]
  gf=gauss1d(xax,p0)
  xr=n_elements(gf)
  
  ;find position for each condition in Gaussian response field
  posrc=nreg*0
  for i=0L,nrnodes-1 do if mask[i] eq 1 then begin
    for j=0,nconds-1 do posrc[i,j]=max(where(nreg[i,j] gt gf))+1
  end
  posrc/=step

  ;Set boundaries for the Levenberg-Marquardt LS fitting procedure
  npar=nconds+4
  parinfo = replicate({fixed:0, limited:[0,0], limits:[0.,0.]}, npar)
  parinfo[4:*].limited(*)=1
  parinfo[4:*].limits(0)=0.
  parinfo[4:*].limits(1)=gcenter
  parinfo[0:1].fixed=1

  ;Define output matrices
  err=fltarr(nscans)+1
  yfitz=funx*0.
  zfitz=fltarr(nrnodes,nconds+3)
  te=100./total(mask)
  count=0L

  ;Levenberg-Marquardt fitting procedure
  for i=0L,nrnodes-1 do if mask[i] eq 1 then begin
    p0=[gcenter,gsigma,median(funx[i,*]),stddev(funx[i,*]),reform(posrc[i,*])]
    zfit=mpfitfun('prftimeorder',fmatrix,reform(funx[i,*]),err,p0,parinfo=parinfo,yfit=yfit,bestnorm=chisq,status=status,/quiet)
    if status gt 0 then begin
      yfitz[i,*]=yfit
      zfitz[i,*]=[zfit[3:*],chisq]
    endif
    count++
    print,te*count
  end

  ;Calculate goodness-of-fit F-statistic
  fval=fltarr(nrnodes)
  for i=0L,nrnodes-1 do if stddev(zfitz[i,*]) ne 0 then fval[i]=(total((mean(funx[i,*])-yfitz[i,*])^2.)/(nconds+1))/(zfitz[i,-1]/((nscans)-(nconds+2)))
  
  ;Add F-statistic to output matrix & save
  zfitz=transpose([transpose(zfitz),transpose(fval)])
  save,zfitz,filename='/Fridge/users/wouter/WholeBodyMap/savs/'+sub+'-zfitz-difpos-FSL.sav'

  ;Save Non-rigid pRF output in freesurfer format and map to average surface
  savemgh,reform(zfitz[0:nrnodesl-1,*]),subdir+'difpos/'+sub+'-zfitz-lh.mgh',hdr=filedir1+sub+'_funcs_'+funrs[0]+'_avg01_lh.mgh'
  savemgh,reform(zfitz[nrnodesl:*,*]),subdir+'difpos/'+sub+'-zfitz-rh.mgh',hdr=filedir1+sub+'_funcs_'+funrs[0]+'_avg01_rh.mgh'
  spawn,'mri_surf2surf --srcsubject '+sub+' --srcsurfval '+subdir+'difpos/'+sub+'-zfitz-lh.mgh --trgsubject avg_wbm --trgsurfval '+subdir+'difpos/'+sub+'-zfitz-lh-avg_wbm.mgh --hemi lh'
  spawn,'mri_surf2surf --srcsubject '+sub+' --srcsurfval '+subdir+'difpos/'+sub+'-zfitz-rh.mgh --trgsubject avg_wbm --trgsurfval '+subdir+'difpos/'+sub+'-zfitz-rh-avg_wbm.mgh --hemi rh'

endfor
end



