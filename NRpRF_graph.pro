pro NRpRF_graph
  ;+
  ; :Author: Wouter Schellekens
  ;
  ; Non-rigid pRF graph theory
  ;-

;Declare variables
rdir='/Fridge/users/wouter/NRPRF/'
fsdir='/Fridge/users/wouter/fs/avg_wbm/'
surfdir=fsdir+'surf/caret/'
subs=['V6649','V6682','V6683','V6743','V6744','V7774','V7860','V7882']
nsub=n_elements(subs)
nconds=18.
nlab=8.
condities = ['toes', 'ankle', 'knee', 'abdomen', 'shoulder', 'elbow','wrist', 'littleF', 'ringF', 'middleF', 'indexF', 'thumb', 'forehead', 'eye', 'nose', 'lips', 'jaw', 'tongue']
nrnodesl=163842L
sm=fltarr(nsub,nrnodesl,24)
gws=fltarr(nsub,nrnodesl,nconds)

;Load non-rigid pRF fits
for i=0,nsub-1 do begin
  sdir=rdir+subs[i]+'/FSL/difpos/'
  lcp=loadmgh(sdir+subs[i]+'-centerpos-lh-avg_wbm.mgh')
  lsp=loadmgh(sdir+subs[i]+'-sigmapos-lh-avg_wbm.mgh')
  zfz=loadmgh(sdir+subs[i]+'-zfitz-lh-avg_wbm.mgh')
  lgw=loadmgh(sdir+subs[i]+'-dxi-lh-avg_wbm.mgh')
  sm[i,*,0] = lcp
  sm[i,*,1] = lsp
  sm[i,*,2:*] = zfz
  gws[i,*,*] = lgw
end

;Make an array of binary masks for each subject
submask=intarr(nsub,nrnodesl)
for i=0,nsub-1 do submask[i,where(finite(sm[i,*,0]) eq 1,/null)]=1

;Set distance beyond the 1 s.d. of the Gaussian response field to zero
gwc=gws
gwc[where(gwc le .5)]=0
subcen=reform(sm[*,*,0])
subcen[where(finite(subcen) eq 1)]=round(subcen[where(finite(subcen) eq 1)])

;Make the mean pRF for each subject, ROI, and condition
nb=fltarr(nsub,nlab,nconds,nconds)
enb=fltarr(nsub,nlab,nconds,nconds)
for i=0,nsub-1 do for j=0,nlab-1 do for k=0,nconds-1 do if where(labmask eq j+1 and round(subcen[i,*]) eq k and submask[i,*] eq 1,/null) ne !NULL then begin
  nb[i,j,k,*]=mean(gwc[i,where(labmask eq j+1 and subcen[i,*] eq k and submask[i,*] eq 1),*],dimension=2,/nan)
  enb[i,j,k,*]=stddev(gwc[i,where(labmask eq j+1 and subcen[i,*] eq k and submask[i,*] eq 1),*],dimension=2,/nan)/sqrt(n_elements(where(labmask eq j+1 and subcen[i,*] eq k and submask[i,*] eq 1)))
endif else begin
  nb[i,j,k,*]=!values.f_nan
  enb[i,j,k,*]=!values.f_nan
endelse
nb2=mean(nb,dimension=1,/nan)

;Calculate the correlation matrix of mean response fields
nbz=nb
nbz[where(finite(nb) eq 0)]=0
cmat=dblarr(nsub,nlab,nconds,nconds)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do for k=0,nconds-1 do if (stddev(nbz[h,i,j,*]) ne 0 and stddev(nbz[h,i,k,*]) ne 0) then cmat[h,i,j,k]=correlate(reform(nbz[h,i,j,*]),reform(nbz[h,i,k,*]),/double)

;Define the lower-bound cutoff at the 95% confidence interval for each ROI
tmp=cmat
thres=dblarr(nlab)
for i=0,nlab-1 do begin
  m=mean(tmp[*,i,*,*],/nan)
  sd=stddev(tmp[*,i,*,*],/nan)
  se=sd/sqrt(nsub)
  thres[i]=m-1.96*se
end

;Make binary matrix where correlation values exceed the threshold
bcmatid=fix(cmat*0)
for i=0,nsub-1 do for j=0,nlab-1 do for k=0,nconds-1 do bcmatid[i,j,k,where(cmat[i,j,k,*] ge thres[j])]=1

;Zeroize binary correlation matrix where mean pRF of each condition correlates to itself
bcmat=bcmatid
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do bcmat[h,i,j,j]=0

;Apply binary correlation matrix (=threshold) to correlation matrix
tcmat=cmat
tcmat*=bcmat

;weighted connectivity strength
S=dblarr(nsub,nlab,nconds)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do S[h,i,j]=total(tcmat[h,i,j,*])

;total conncectivity strength
Stot=mean(S,dimension=3)

;total possible number of triangles per condition
kn=0.
count=0L
for i=3,nconds do begin
  count++
  kn+=count
end

;the actual triangles per condition
tri=intarr(nsub,nlab,nconds,kn,2)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do if total(bcmat[h,i,j,*]) ge 2 then begin
  count=0
  tmp=where(bcmat[h,i,j,*] eq 1)
  for k=0,n_elements(tmp)-1 do for l=0,n_elements(tmp)-1 do if l gt k then if bcmat[h,i,tmp[k],tmp[l]] eq 1 then begin
    tri[h,i,j,count,*]=[tmp[k],tmp[l]]
    count++
  endif
end

;weighted geometric mean of triangles
gmt=dblarr(nsub,nlab,nconds,kn)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do for k=0,kn-1 do if total(tri[h,i,j,k,*]) gt 0 then gmt[h,i,j,k]=(tcmat[h,i,j,tri[h,i,j,k,0]]*tcmat[h,i,j,tri[h,i,j,k,1]]*tcmat[h,i,tri[h,i,j,k,0],tri[h,i,j,k,1]])^(1/3d)
gmt2=total(gmt,4)/2d

;weighted clustering coefficent
nrtri=intarr(nsub,nlab,nconds)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do nrtri[h,i,j]=n_elements(where(total(tri[h,i,j,*,*],5) gt 0,/null))
wcc=gmt2*0d
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do if nrtri[h,i,j] ne 0 then wcc[h,i,j]=(2d*gmt2[h,i,j])/(total(bcmat[h,i,j,*])*(total(bcmat[h,i,j,*])-1d))
twcc=mean(wcc,dimension=3,/double)

;Find the shortest paths
wmat=1d/tcmat
allpathlengths=dblarr(nsub,nlab,nconds,nconds,3)
allpaths=lonarr(nsub,nlab,nconds,nconds,nconds)
for h=0,nsub-1 do for i=0,nlab-1 do for j=0,nconds-1 do begin
  pathfinder_matrix,[j],where(findgen(nconds) ne j),reform(bcmat[h,i,*,*]),reform(wmat[h,i,*,*]),paths,pathlengths
  allpaths[h,i,j,where(findgen(nconds) ne j),*]=paths
  allpathlengths[h,i,j,where(findgen(nconds) ne j),*]=pathlengths
end
apl=reform(allpathlengths[*,*,*,*,0])
apl[where(apl eq 0)]=!values.f_nan

;Weighted betweenness centrality
bc=dblarr(nsub,nlab,nconds)
nodepath=dblarr(nconds)
for h=0,nsub-1 do for i=0,nlab-1 do begin
  totpath=n_elements(where(allpathlengths[h,i,*,*,0] ne 0))/2d
  nodepath*=0d
  for j=0,nconds-1 do begin
    nodepath[j]=n_elements(where(allpaths[h,i,where(findgen(nconds) ne j),where(findgen(nconds) ne j),j] ne 0,/null))
    bc[h,i,j]=nodepath[j]/totpath
  endfor
end

;modularity Louvain algorithm over average subjects
cmats=dblarr(nlab,nconds,nconds)
for i=0,nlab-1 do cmats[i,*,*]=correlate(reform(nb2[i,*,*]),/double)

;threshold average correlation matrix
tcmats=cmats
for i=0,nlab-1 do begin
  for j=0,nconds-1 do begin
    tcmats[i,j,where(cmats[i,j,*] gt thres[i])]-=thres[i]
    tcmats[i,j,where(cmats[i,j,*] le thres[i])]=0
    tcmats[i,j,j]=0
  endfor
  tcmats[i,*,*]/=max(tcmats[i,*,*])
end

;Binary average correlation matrix
bcmats=fix(cmats*0)
for i=0,nlab-1 do for j=0,nconds-1 do begin
  bcmats[i,j,where(cmats[i,j,*] gt thres[i])]=1
  bcmats[i,j,j]=0
end

modul=lonarr(nlab,nconds,nconds)
for i=0,nlab-1 do modul[i,*,*]=modularity_louvain(reform(bcmats[i,*,*]),reform(tcmats[i,*,*]))

lname=['M1','S1','SMA','PMd','PMv','Ins','iPC','sPC']
cname=['Toe','Ank','Kne','Abd','Sho','Elb','Wri','LiF','RiF','MiF','InF','Thu','For','Eye','Nos','Lip','Jaw','Ton']
tits=strarr(nlab-1,nconds)
for i=0,nlab-1 do for j=0,nconds-1 do tits[i,j]=lname[i]+'-'+cname[j]

;output graph theory metrics
graphw=strarr((nlab-1)*nconds,nsub+1)
graphw[*,0]=reform(tits,(nlab-1)*nconds)
tmp=transpose(S,[1,2,0])
graphw[*,1:*]=string(reform(tmp[0:nlab-1,*,*],(nlab-1)*nconds,nsub))
write_ascii,rdir+'group/graph-connectivity.txt',strcompress(graphw,/remove_all)
save,S,filename=rdir+'group/graph-connectivity.sav'

graphw=strarr((nlab-1)*nconds,nsub+1)
graphw[*,0]=reform(tits,(nlab-1)*nconds)
tmp=transpose(bc,[1,2,0])
graphw[*,1:*]=string(reform(tmp[0:nlab-1,*,*],(nlab-1)*nconds,nsub))
write_ascii,rdir+'group/graph-betweennesscentrality.txt',strcompress(graphw,/remove_all)
save,bc,filename=rdir+'group/graph-betweennesscentrality.sav'

graphw=strarr((nlab-1)*nconds,nsub+1)
graphw[*,0]=reform(tits,(nlab-1)*nconds)
tmp=transpose(wcc,[1,2,0])
graphw[*,1:*]=string(reform(tmp[0:nlab-1,*,*],(nlab-1)*nconds,nsub))
write_ascii,rdir+'group/graph-clustering.txt',strcompress(graphw,/remove_all)
save,wcc,filename=rdir+'group/graph-clustering.sav'



;Make figures

;Graph images
bcmatall=bcmats
points=[[300,100],[400,100],[400,200],[500,300],[400,400],[300,300],[300,400],[100,400],[100,500],[100,600],[200,600],[300,600],[400,800],[400,700],[500,700],[500,600],[500,500],[400,600]]
ct=colortable([[255,0,0],[255,255,0],[0,255,0],[0,0,255],[175,25,85]],ncolors=fix(nconds))
ctbizar=[[0,137,171],[50,128,0],[0,0,0],[202,53,53]]
lnameg=['M1','S1','SMA','PMd','PMv','Insula','iPC','sPC']

linethick=fltarr(nlab,nconds,nconds)
for i=0,nlab-1 do begin
  hg=histogram(tcmats[i,*,*],locations=loc,nbins=4)
  for j=0,nconds-1 do for k=0,nconds-1 do linethick[i,j,k]=loc[min(where(tcmats[i,j,k] le loc))]*9
end
linethick/=2.5

xd=600.
yd=900.
xs=40.
sxd=8*xd+7*xs+100
syd=yd

points2=fltarr(2,nconds,nlab-1)
for i=0,nconds-1 do for j=0,nlab-1 do begin
  points2[0,i,j]=points[0,i]+(j*(xd+xs))+100
  points2[1,i,j]=points[1,i]
end

;Connectivity graph
Sall=mean(S,dimension=1)
hg=histogram(Sall,locations=locc,nbins=5)

;rescale to node sizes to fit image window
Sall/=max(Sall)
Sall*=50d

graphim=intarr(sxd+300,syd)+!values.f_nan
im=image(graphim,margin=0,background_color='white',dimensions=[2500,500])
for h=0,nlab-1 do begin
  for i=0,nconds-1 do for j=0,nconds-1 do if j gt i and bcmatall[h,i,j] eq 1 then pl=polyline([points2[0,i,h],points2[0,j,h]],[points2[1,i,h],points2[1,j,h]],thick=linethick[h,i,j],color='black',/data)
  for i=0,nconds-1 do e = ELLIPSE(reform(points2[0,i,h]), reform(points2[1,i,h]),MAJOR=Sall[h,i],/DATA, FILL_COLOR=reform(ct[i,*]),thick=2)
  tt=text(mean(points2[0,*,h])-200,mean(points2[1,*,h])+400,lnameg[h],font_size=24,font_style=1,/data)
end

for i=5,51,10 do begin
  num=string(locc[i/10],FORMAT='(F5.2)')
  e = ELLIPSE(sxd+220,i*12+syd/5.5,MAJOR=i,/DATA, FILL_COLOR=reform(ctbizar[*,3]),COLOR=reform(ctbizar[*,3]),thick=2)
  t=text(sxd+20,i*12-20+syd/5.5,/data,num,font_size=24)
end

t=text(70,syd/3.3,'Connectivity',/data,font_size=24,font_style=1,orientation=90)

im.save,rdir+'group/graphs/graph-connectivity.png'
im.close

;Clustering graphs
CCall=mean(wcc,dimension=1)
hg=histogram(CCall,locations=loccc,nbins=5)

;rescale to node sizes to fit image window
CCall/=max(CCall)
CCall*=50d

graphim=intarr(sxd+300,syd)+!values.f_nan
im=image(graphim,margin=0,background_color='white',dimensions=[2500,500])
for h=0,nlab-1 do begin
  for i=0,nconds-1 do for j=0,nconds-1 do if j gt i and bcmatall[h,i,j] eq 1 then pl=polyline([points2[0,i,h],points2[0,j,h]],[points2[1,i,h],points2[1,j,h]],thick=linethick[h,i,j],color='black',/data)
  for i=0,nconds-1 do e = ELLIPSE(reform(points2[0,i,h]), reform(points2[1,i,h]),MAJOR=CCall[h,i],/DATA, FILL_COLOR=reform(ct[i,*]),thick=2)
  tt=text(mean(points2[0,*,h])-200,mean(points2[1,*,h])+350,lnameg[h],font_size=24,font_style=1,/data)
end

for i=5,51,10 do begin
  num=string(loccc[i/10],FORMAT='(F5.2)')
  e = ELLIPSE(sxd+220,i*12+syd/5.5,MAJOR=i,/DATA, FILL_COLOR=reform(ctbizar[*,0]),COLOR=reform(ctbizar[*,0]),thick=2)
  t=text(sxd+20,i*12-20+syd/5.5,/data,num,font_size=24)
end

t=text(70,syd/3.3,'Clustering',/data,font_size=24,font_style=1,orientation=90)

im.save,rdir+'group/graphs/graph-clustering.png'
im.close

;Betweenness Centrality graphs
BCall=mean(bc,dimension=1)
hg=histogram(BCall,locations=locbc,nbins=5)

;rescale to node sizes to fit image window
BCall/=max(BCall)
BCall*=50d

graphim=intarr(sxd+300,syd)+!values.f_nan
im=image(graphim,margin=0,background_color='white',dimensions=[2500,500])
for h=0,nlab-1 do begin
  for i=0,nconds-1 do for j=0,nconds-1 do if j gt i and bcmatall[h,i,j] eq 1 then pl=polyline([points2[0,i,h],points2[0,j,h]],[points2[1,i,h],points2[1,j,h]],thick=linethick[h,i,j],color='black',/data)
  for i=0,nconds-1 do e = ELLIPSE(reform(points2[0,i,h]), reform(points2[1,i,h]),MAJOR=BCall[h,i],/DATA, FILL_COLOR=reform(ct[i,*]),thick=2)
  tt=text(mean(points2[0,*,h])-200,mean(points2[1,*,h])+350,lnameg[h],font_size=24,font_style=1,/data)
end

for i=5,51,10 do begin
  num=string(locbc[i/10],FORMAT='(F5.2)')
  e = ELLIPSE(sxd+220,i*12+syd/5.5,MAJOR=i,/DATA, FILL_COLOR=reform(ctbizar[*,2]),COLOR=reform(ctbizar[*,2]),thick=2)
  t=text(sxd+20,i*12-20+syd/5.5,/data,num,font_size=24)
end

t=text(70,syd/3.3,'Centrality',/data,font_size=24,font_style=1,orientation=90)

im.save,rdir+'group/graphs/graph-Centrality.png'
im.close


;Modules
modul2=total(modul,3)
cct=colortable([[255,0,0],[255,255,0],[0,255,0],[0,255,255],[0,0,255],[255,0,255]],ncolors=7)
graphim=intarr(sxd,syd)+!values.f_nan
im=image(graphim,margin=0,background_color='white',dimensions=[2250,500])
for h=0,nlab-1 do begin
  cmod=where(modul2[h,*] ne 0)
  cct=colortable([[255,0,0],[255,255,0],[0,255,0],[0,255,255],[0,0,255],[255,0,255]],ncolors=n_elements(cmod))
  for i=0,nconds-1 do for j=0,nconds-1 do if j gt i and bcmatall[h,i,j] eq 1 then pl=polyline([points2[0,i,h],points2[0,j,h]],[points2[1,i,h],points2[1,j,h]],thick=linethick[h,i,j],color='black',/data)
  for i=0,nconds-1 do begin
  tmp=where(modul[h,*,i] eq 1)
  e = ELLIPSE(reform(points2[0,i,h]), reform(points2[1,i,h]),MAJOR=30,/DATA, FILL_COLOR=reform(cct[where(tmp[0] eq cmod),*]),thick=2)
  endfor
  tt=text(mean(points2[0,*,h])-200,mean(points2[1,*,h])+350,lnameg[h],font_size=24,font_style=1,/data)
end

t=text(70,syd/3.3,'Modules',/data,font_size=24,font_style=1,orientation=90)

im.save,rdir+'group/graphs/graph-Modules.png'
im.close

end

