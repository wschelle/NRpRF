function prftimeorder,x,p

  syntx=n_params()
  if syntx lt 2 then begin
    print,'Syntax not right; parms are:'
    print,'1: independent variables'
    print,'2: parameters: [center,sigma,constant,amplitude,[positions]]'
    goto, EINDE
  endif
  
  sx=size(x)
  x2=x*0
  for i=0,sx[1]-1 do x2[i,*]=x[i,*]*exp(-(p[i+4]-p[0])^2./(2.*p[1]^2.))
  x2=total(x2,1)
  x2=p[2]+p[3]*x2
  return, x2
  
  EINDE: 
end