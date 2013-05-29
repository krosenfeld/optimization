@utils_sampler.pro
pro gw10_sampler,ensemble=ensemble,info=info,lnprob_ensemble=lnprob_ensemble

; import user settings
@gaussian_inputs.pro    ; N-dimensional gaussian example
;@rosenbrock_inputs.pro ; Rosenbrock function example

a	     = 2.   ; stretch parameter (see Foreman-Mackay+12)

; setup amount of feedback 
silent   = 0            ; 1: NO standard i/o
chatter  = 0            ; 1: LOTS of standard i/o
if (chatter eq 1) then talk = 0
if (chatter eq 0) then talk = 1
if (silent  eq 1) then talk = 2

; set up storage
if (writeme) then begin
  format = '('+string(ndim,format='(I0)')+'(e12.5,2x,e12.5))'
  cfile = dir+root+'.dat'
  status = file_exist(cfile)
  if (status eq 1) then begin
    print,'WARNING: '+cfile+' exists'
    stop,'.c to overwrite and continue'
  endif
  spawn,'rm -f '+cfile
endif

; start timer
starttime = systime(/seconds) 

; initialize walkers
ensemble        = dblarr(ndim,nwalkers,nsamples)   ; save samples
lnprob_ensemble = dblarr(nwalkers,nsamples)        ; save lnprob
ensemble[*,*,0] = rebin(mpars,ndim,nwalkers) + rebin(spars,ndim,nwalkers)*randomn(seed,ndim,nwalkers)
nacpt           = 0L

print, '================INITIALIZE==================='
; calculate log probabilities for starting ensemble
lnprob_S = dblarr(nwalkers)
for ik=0,nwalkers-1 do begin
   Y = ensemble[*,ik,0]
   lnprob_S[ik] = call_function(lnprobfn,Y,_extra=extra)
   lnprob_ensemble[ik,0] = lnprob_S[ik]
endfor

; report
print,'using '+string(nwalkers,format='(I0)')+$
      ' walkers to sample '+lnprobfn

; draw samples
f_sample = 0
if (talk le 1) then print, '==================='+'DRAW SAMPLES'+'======================'
for i=1L,nsamples-1 do begin
 if (talk eq 0) then print, '================= '+string(i,format='(I6)')+'===================='
 ; update walkers
 ensemble[*,*,i] = ensemble[*,*,i-1]

 for ik=0L,nwalkers-1 do begin
   if (talk eq 0) then print,'WALKER #'+string(ik,format='(I3)')

   ; our walker 
   Xk = ensemble[*,ik,i]

   ; randomly choose other walker
   ij = floor(randomu(seed)*nwalkers)
   while (ij eq ik) do ij = floor(randomu(seed)*nwalkers)
   Xj = ensemble[*,ij,i]

   ; propose step
   zz = ((a - 1.)*randomu(seed) + 1.)^2/a
   Y = Xj + zz*(Xk-Xj)

   ; acceptance probability
   lnprob  = call_function(lnprobfn,Y,_extra=extra)

   ; calculate log difference
   lnpdiff = (ndim-1)*alog(zz) + lnprob - lnprob_S[ik]
   u       = alog(randomu(seed))

   if (u lt lnpdiff) then begin
     ; accept
     ensemble[*,ik,i]      = Y 
     lnprob_ensemble[ik,i] = lnprob
	 lnprob_S[ik]          = lnprob
	 nacpt = nacpt + 1
   endif else $
     lnprob_ensemble[ik,i] = lnprob_S[ik] ; else reject
  endfor

  ; check i/o duty cycle
  if (writeme) then begin
    if (i mod ioduty eq 0) then begin
      l_sample = f_sample + ioduty
      openw,1,cfile,/append
      for ij=f_sample,l_sample-1 do begin
        for ik=0,nwalkers-1 do begin
          printf,1,ensemble[*,ik,ij],format=format
        endfor
      endfor
      close,1
      f_sample = l_sample
    endif
  endif

 ; report acceptance ratio
 if (talk le 1 and (i mod cduty) eq 0) then $
        print,'iteration '+string(i,format='(I0)')+', accept ratio = ' + string(nacpt/(1d*i*nwalkers),format='(F5.3)')

 ; check convergence
 if (i mod cduty eq 0) then begin
  cflag = 1
  _iburn = 0
  ; check burn in
  for idim=0,info.ndim-1 do begin
    s = reform(mean(ensemble[idim,*,0:i-1],dim=2),i)
    autoburn,s,info,_iburn=foo,ithin=1,minburn=info.mburn,status=status
    if status ne 0 then idim = info.ndim-1
    if (foo gt _iburn) then _iburn = foo
  endfor
  ; check autocorrelation time
  if (status eq 0) then begin
    for idim=0,info.ndim-1 do begin
      s = reform(mean(ensemble[idim,*,0:i-1],dim=2),i)
      status = acor(s,info,i,ithin=1,iburn=_iburn,$
                tau=tau,sigma=sigma,_C=C)
      if (status ne 0 or i*nwalkers lt acormult*tau) then cflag = 0
    endfor
  endif else cflag = 0
  if (cflag eq 1) then begin
    itot = i
    i = nsamples
  endif
 endif
endfor

if (cflag eq 0) then itot = nsamples

; save last samples 
if (writeme) then begin
  l_sample = itot
  openw,1,cfile,/append
  for ij=f_sample,l_sample-1 do begin
    for ik=0,nwalkers-1 do begin
      printf,1,ensemble[*,ik,ij],lnprob_ensemble[ik,ij],format=format
    endfor
  endfor
  close,1
endif

; trim result
ensemble=ensemble[*,*,0:itot-1]
lnprob_ensemble=lnprob_ensemble[*,0:itot-1]
info.nsamples = itot

; wrap up
if (talk le 1) then begin
  print, '==================='+'REPORT'+'======================'
  ; save result (if desired)
  if (saveme) then begin
    print,'saving...'
    save,info,ensemble,lnprob_ensemble,filename=lnprobfn+'.res'
  endif
  ; stop timer
  stoptime = systime(/seconds)
  if (cflag) then $
       print,'sampler converged in '+string(itot,format='(I0)')+' steps' $
  else print,'WARNING: sampler not converged.'
  print,'sampler took: '+string((stoptime - starttime)/60.,format='(F8.2)')+' minutes'
  print,'accept ratio:', float(nacpt)/(itot*nwalkers)
endif
end
