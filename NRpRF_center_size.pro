pro NRpRF_center_size

  ;+
  ; :Author: Wouter Schellekens
  ;
  ; Non-rigid pRF center & size calculation
  ;-
  
  ;Declare variables
  rdir='/Fridge/users/wouter/NRPRF/'
  fsdir='/Fridge/users/wouter/fs/avg_wbm/'
  surfdir=fsdir+'surf/caret2/'
  subs=['V6649','V6682','V6683','V6743','V6744','V7774','V7860','V7882']
  nsub=n_elements(subs)
  nconds=18.
  condities = ['toes', 'ankle', 'knee', 'belly', 'shoulder', 'elbow','wrist', 'littleF', 'ringF', 'middleF', 'indexF', 'thumb', 'forehead', 'eye', 'nose', 'lips', 'jaw', 'tongue']
  nscans=856.
  
  ;Define the static Gaussian function
  step=100.
  gcenter=10.
  gsigma=1.
  xax=findgen(gcenter*step)/step
  p0=[0.,1.,gsigma,gcenter]
  gf=gauss1d(xax,p0)
  xr=n_elements(gf)
  
  ;Load pRF fits
  for i=0,nsub-1 do begin
    subdir=rdir+subs[i]+'/FSL/'
    zfitz=loadmgh(subdir+'difpos/'+subs[i]+'-zfitz-lh.mgh')
    nrnodes=n_elements(zfitz[*,0])
    
    ;Set values beyond 1 s.d. to zero
    zfitz2=reform(zfitz[*,2:2+nconds-1])
    hwhm=sqrt(2*alog(2))
    fwhm=2*hwhm
    zfitz3=zfitz2
    low=fwhm
    zfitz3-=low
    zfitz3[where(zfitz2 lt gcenter-hwhm)]=0
    zfitz3/=(gcenter-low)
   
    ;Threshold non-significant vertices
    mask=intarr(nrnodes)
    mask[where(zfitz[*,-1] gt 0)]=1
    fpar=nconds+2
    fthr=f_fdr(reform(zfitz[*,-1]),[fpar-1,nscans-fpar-1])
    m2=mask
    m2[where(reform(zfitz[*,-1]) lt fthr)]=0
    
    ;Calculate pRF center  
    centerpos=fltarr(nrnodes)
    for ii=0L,nrnodes-1 do if m2[ii] eq 1 then begin
      if n_elements(where(zfitz2[ii,*] eq max(zfitz2[ii,*]),/null)) eq 1 then begin
        centerpos[ii]=where(zfitz2[ii,*] eq max(zfitz2[ii,*]))
      endif else begin
        tmp=where(zfitz2[ii,*] eq max(zfitz2[ii,*]))
        tmp2=dblarr(n_elements(tmp))
        for jj=0,n_elements(tmp)-1 do for kk=-1,1,2 do if tmp[jj]+kk ge 0 and tmp[jj]+kk lt nconds then if zfitz2[ii,tmp[jj]+kk] gt tmp2[jj] then tmp2[jj]=zfitz2[ii,tmp[jj]+kk]
        centerpos[ii]=tmp[where(tmp2 eq max(tmp2))]
      endelse
    endif

    ;Calculate pRF size
    zfitz4=zfitz2/gcenter
    sigmapos=fltarr(nrnodes)
    for ii=0L,nrnodes-1 do if m2[ii] eq 1 then sigmapos[ii] = total(zfitz4[ii,where(zfitz3[ii,*] gt 0)])*hwhm
    
    ;Set non-significant vertices to NaN
    centerpos[where(m2 eq 0)]=!values.f_nan
    sigmapos[where(m2 eq 0)]=!values.f_nan
    zfitz3[where(m2 eq 0),*]=!values.f_nan
    
    ;Save it to Freesurfer format and project to average cortical surface
    savemgh,centerpos,subdir+'difpos/'+subs[i]+'-centerpos-lh.mgh'
    spawn,'mri_surf2surf --srcsubject '+subs[i]+' --srcsurfval '+subdir+'difpos/'+subs[i]+'-centerpos-lh.mgh --trgsubject avg_wbm --trgsurfval '+subdir+'difpos/'+subs[i]+'-centerpos-lh-avg_wbm.mgh --hemi lh'
    
    savemgh,sigmapos,subdir+'difpos/'+subs[i]+'-sigmapos-lh.mgh'
    spawn,'mri_surf2surf --srcsubject '+subs[i]+' --srcsurfval '+subdir+'difpos/'+subs[i]+'-sigmapos-lh.mgh --trgsubject avg_wbm --trgsurfval '+subdir+'difpos/'+subs[i]+'-sigmapos-lh-avg_wbm.mgh --hemi lh'

    savemgh,zfitz3,subdir+'difpos/'+subs[i]+'-dxi-lh.mgh'
    spawn,'mri_surf2surf --srcsubject '+subs[i]+' --srcsurfval '+subdir+'difpos/'+subs[i]+'-dxi-lh.mgh --trgsubject avg_wbm --trgsurfval '+subdir+'difpos/'+subs[i]+'-dxi-lh-avg_wbm.mgh --hemi lh'

  endfor

end

