; ***********************************************************
; autoburn: 
;   automatically calculate burn in 
; 
;  - 5/13 : created by KR
;
pro autoburn,x,info,_iburn=_iburn,ithin=ithin,pauseme=pauseme,minburn=minburn,status=status
  if (not keyword_set(ithin)) then ithin = 1    ; no thinning
  if (not keyword_set(pauseme)) then pauseme=0  ; no pausing
  if (not keyword_set(minburn)) then minburn=0  ; no minimum burn

  ; check array size
  status = 0
  _x = reform(x,n_elements(x))
  nx = n_elements(_x)
  if (minburn gt nx) then status=1

  ; calculate burn in
  if (status eq 0) then begin
    cx = total(x,/cum)/(1+findgen(nx))
    ; mean array (from end)
    mrx = total(reverse(_x),/cum,/dou)/(1+findgen(nx))
    ; std array (from end)
    srx = sqrt(total((reverse(_x) - mrx)^2,/cum,/dou)/(1+findgen(nx)))

    flag    = 0
    _iburn  = minburn
    while flag eq 0 do begin
      mx = mrx[nx-_iburn] 
      sx = srx[nx-_iburn]
      foo = abs(_x[0:_iburn] - mx)
      ii  = where(foo gt sx)
      if ((1d*n_elements(ii)/_iburn gt 0.7) and (_iburn lt nx-1)) then begin
        _iburn++
      endif else flag = 1
    endwhile
    ; make sure we haven't burned away the series
    if (_iburn eq nx-1) then status=1
  endif
  if (pauseme) then stop
end
; ***********************************************************
; pl_means:
;   plot the moving mean of the ensemble
; 
;  - 5/13 : created by KR
;
pro pl_means,ensemble,info,iburn=iburn,ithin=ithin,pmulti=pmulti,zoomit=zoomit
  if (not keyword_set(iburn)) then iburn = -1       ; use automatic burn
  if (not keyword_set(ithin)) then ithin = 1        ; no thinning
  if (not keyword_set(pmulti)) then pmulti = [2,2]  ; 2x2 figure matrix
  if (not keyword_set(zoomit)) then zoomit = 0      ; no zoom

  pl_mkct

;
  iithin = iburn+ithin*lindgen(floor((info.nsamples-iburn)/ithin))
  nplots = pmulti[0]*pmulti[1]
  sz = size(ensemble)
  print,'plotting moving averages:'
  for idim=0,info.ndim-1 do begin
    if ((idim) mod nplots eq 0) then !p.multi = [0,pmulti[0],pmulti[1]]
;    print,'showing '+info.pname[idim]+$
;          '; mean:'+string(mu,format='(E10.2)')
    if (iburn lt 0) then begin
      _iburn = 0
      iithin = _iburn+ithin*lindgen(floor((info.nsamples-_iburn)/ithin))
      autoburn,mean(ensemble[idim,*,iithin],dim=2),info,_iburn=_iburn,minburn=info.mburn,status=status
      if status ne 0 then stop
      iithin = _iburn+ithin*lindgen(floor((info.nsamples-_iburn)/ithin))
      xx     = _iburn+findgen(n_elements(iithin))
      ; report
;      print,'automatic burn: '+string(_iburn,format='(I0)')
;      stop
    endif else $
      xx = iburn+findgen(n_elements(iithin))

    print,'showing '+info.pname[idim]
    mu = mean(ensemble[idim,*,iithin],dim=2)
    foo = ithin*lindgen(floor((info.nsamples)/ithin))
    if (zoomit) then $
      plot,xx,mu,/nodata,$
        xtit='N',ytit=info.pname[idim],/xsty,xra=[0,max(xx)],/ysty $
    else $
      plot,foo,mean(ensemble[idim,*,foo],dim=2),/nodata,$
        xtit='N',ytit=info.pname[idim],/xsty,xra=[0,max(xx)],/ysty 
    ; overplot cumulative mean
    oplot,xx,total(mu,/cum)/(1+findgen(n_elements(mu))),col=3
    ; overplot std
    oplot,!X.crange,(mean(mu) + stddev(mu)) * [1,1]
    oplot,!X.crange,(mean(mu) - stddev(mu)) * [1,1]
    ; show burned samples
    oplot,findgen(n_elements(foo)),mean(ensemble[idim,*,foo],dim=2),col=120
    ; show trimmed series
    oplot,xx,mu

    if ((idim+1) mod nplots eq 0) then stop
  endfor
  print,'done!'
  !p.multi=0
end
;
; ***********************************************************
; pl_histo2D:
;   plot 2D histograms
;
;  - 5/13 : created by KR
; 
pro pl_histo2D,ensemble,info,nbins=nbins,iburn=iburn,ithin=ithin,CI=CI
  if (not keyword_set(iburn))  then iburn = -1  ; automatic burn
  if (not keyword_set(ithin))  then ithin = 1   ; no thinning
  if (not keyword_set(nbins))  then nbins = 50  ; # of bins
  if (not keyword_set(CI))     then CI = [1]    ; sigma for CI

  nplot = info.ndim - 1
  foo = floor(sqrt(nplot))
  pmulti = [ceil(1d*nplot/foo),foo]

  pl_mkct

  ; calcualte x burn in 
  print,'calculating burn in...'
  _iburn = 0
  if (iburn lt 0) then begin
    for idim=0,info.ndim-1 do begin
       s = reform(mean(ensemble[idim,*,*],dim=2),info.nsamples)
       autoburn,s,i,minburn=info.mburn,ithin=ithin,_iburn=foo
       if (foo gt _iburn) then _iburn = foo
    endfor
  endif else _iburn = iburn
  iithin = _iburn+ithin*lindgen(floor((info.nsamples-_iburn)/ithin))
  print,'burning:',_iburn

  ; display histograms
  for idim=0,info.ndim-1 do begin
    !p.multi = [0,pmulti[0],pmulti[1]]
    ; calculate x info 
    x = reform(ensemble[idim,*,iithin],info.nwalkers*n_elements(iithin))
    ivar = 0
    min_x  = min(x)
    max_x  = max(x) 
    res_x  = nbins
    n = n_elements(x)
    bin_x = (max_x - min_x)/(res_x)
    binedge_x = bin_x*findgen(res_x) + min_x
    axis_x = 0.5*bin_x + binedge_x
    ; go through other variables
    for ix=0,pmulti[1] do begin
      for iy=0,pmulti[0] do begin
        if (ivar ne idim) then begin
          y = reform(ensemble[ivar,*,iithin],info.nwalkers*n_elements(iithin))
 
          min_y  = min(y)
          max_y  = max(y) 
          res_y  = nbins 
          bin_y = (max_y - min_y)/(res_y)
          binedge_y = bin_y*findgen(res_y) + min_y
          axis_y = 0.5*bin_y + binedge_y
          i_offplot = where(x gt max_x or y gt max_y or x lt min_x or y lt min_y, $
                      n_offplot, complement=i_onplot)
          h = hist_2d(value_locate(binedge_x,x[i_onplot]),$
                      value_locate(binedge_y, y[i_onplot]), $
                             min1=0,max1=res_x,min2=0,max2=res_y)
          hh = histogram(h,reverse_indices=ri,min=0,locations=locations)

          confidence_levels = erf(CI/sqrt(2))
          levels = [value_locate(total((hh*locations),/cumulative), $
              (1.0-confidence_levels)*n-n_offplot)]
          i = where(levels ge 0)
          i = i[uniq(levels[i])]

         contour,h[1:res_x-1, 1:res_y-1],axis_x[1:res_x-1],axis_y[1:res_y-1], $
            levels=levels[i],c_thick=thick,/close
        endif
        if ivar eq info.ndim-1 then begin
          iy = pmulti[0]
          ix = pmulti[1]
        endif
        ivar++
      endfor
    endfor
    print,info.pname[idim]
    stop
  endfor
  !p.multi = 0
end
;
; ***********************************************************
; pl_histo1D:
;   plot 1D histograms and report confidence intervals
;
;  - 5/13 : created by KR
; 
pro pl_histo1D,ensemble,info,nbins=nbins,iburn=iburn,ithin=ithin,pmulti=pmulti,$
    std=std
  if (not keyword_set(iburn))  then iburn = 0
  if (not keyword_set(ithin))  then ithin = 1
  if (not keyword_set(nbins))  then nbins = 50
  if (not keyword_set(std))    then std = 1.
  if (not keyword_set(pmulti)) then pmulti = 1
  if (n_elements(pmulti) eq 1) then begin
    foo = floor(sqrt(info.ndim))
    pmulti = [ceil(1d*info.ndim/foo),foo]
  endif
;
  nplots = pmulti[0]*pmulti[1]
  tol    = erf(std/sqrt(2.))   ; tolerance for CI
  sz = size(ensemble)

  ; calculate burn in
  print,'calculating burn in...'
  _iburn = 0
  if (iburn lt 0) then begin
    for idim=0,info.ndim-1 do begin
        s = reform(mean(ensemble[idim,*,*],dim=2),sz[3])
        autoburn,s,i,minburn=info.mburn,ithin=ithin,_iburn=foo
        if (foo gt _iburn) then _iburn = foo
    endfor
  endif else _iburn = iburn
  print,'burning:',_iburn

  ;calculate histgram
  print,'mean, [CI]'
  for idim=0,info.ndim-1 do begin
    iithin = _iburn+ithin*lindgen(floor((info.nsamples-_iburn)/ithin))
    s = ensemble[idim,*,iithin]
    ; calculate histogram
    hist1d = histogram(s,nbins=nbins,loc=binloc)
    bmax   = max(s,min=bmin)
    mbins  = binloc + 0.5*(bmax-bmin)/(nbins-1)
    lbins  = binloc
    rbins  = binloc + (bmax-bmin)/(nbins-1)
;
    ; calculate CI
    htot = total(hist1d,/dou)
    mu = mean(s)
;    foo  = max(hist1d,imax)         ; use peak
    imax = value_locate(binloc,mu)  ; use mean
    moo = hist1d
    moo[imax] *= 0.5
    lfoo = (total(reverse(moo[0:imax]),/cum,/dou))/htot
    ufoo = (total(        moo[imax:-1],/cum,/dou))/htot
    iCI  = [imax-max(where(lfoo lt 0.5*tol)),imax+max(where(ufoo lt 0.5*tol))]
    CI   = [lbins[imax-max(where(lfoo lt 0.5*tol))],rbins[imax+max(where(ufoo lt 0.5*tol))]]
;
    if ((idim) mod nplots eq 0) then !p.multi = [0,pmulti[0],pmulti[1]]
    ;print,info.pname[idim],CI[0],mu,CI[1]
    print,info.pname[idim],mu,'  [',CI[0],',',CI[1],']'
    plot,mbins,hist1d/htot,$
        ytit='N',xtit=info.pname[idim],/xsty,psym=10
    for i=0,1 do $
        oplot,CI[i]*[1,1],!Y.crange,linesty=2
    if ((idim+1) mod nplots eq 0) then stop
  endfor
  !p.multi = 0
  print,'done!'
end
;
; ***********************************************************
; acor:
;   calculate autocorrelation time following J. Goodman's "acor" program.
; 
;   reference: http://www.math.nyu.edu/faculty/goodman/software/acor/ 
;
;  - 5/13 : created by KR
;
function acor,ensemble,info,N,ithin=ithin,iburn=iburn,maxlag=maxlag,winmult=winmult,$
    tau=tau,sigma=sigma,_C=_C
  forward_function acor

  if (not keyword_set(iburn))   then iburn   = 0
  if (not keyword_set(ithin))   then ithin   = 1
  if (not keyword_set(maxlag))  then maxlag  = 10
  if (not keyword_set(winmult)) then winmult = 2    ; is 5 in acor

  iithin = iburn+ithin*lindgen(floor((N-iburn)/ithin))
  sz     = size(ensemble)
  L      = n_elements(iithin)
  C      = dblarr(maxlag+1)
  imax   = L-maxlag-1
  mu     = mean(ensemble[iithin])

  minfac = 2
  ;if (10 gt L) then return,1   ; NOT CONVERGED: autocorrelation time is too long 
  if (minfac*maxlag gt L) then return,1   ; NOT CONVERGED: autocorrelation time is too long 
                                          ;                relative to variance

  ; calcualte autocorrelation
  for it=0L,maxlag do $
    C[it] = total(    (      (ensemble[iithin])       [0:imax]-mu)*$
                      ((shift(ensemble[iithin],-1*it))[0:imax]-mu),/dou) 
  C /= (imax+1)

  _C = C    ; save autocorrelation

  ; calculate diffusion coefficient (sum of autocovariances)
  D = C[0]
  for i=1L,maxlag do D += 2*C[i]
  ;sigma = sqrt(D/L)
  sigma = sqrt(abs(D/L))
  if (finite(sigma,/nan)) then return,2 ; NOT CONVERGED: D < 0

  ; calculate integrated correlation time
  tau = D/C[0] 

  if (tau*winmult lt maxlag) then begin
    return,0
  endif else begin
    Lh  = floor(0.5*L)
    ;if (Lh le 2*maxlag) then return,3  ; NOT CONVERGED: acor is not converging
    if (Lh le maxlag) then return,3  ; NOT CONVERGED: acor is not converging
    X   = ensemble[iithin] + $
           shift((ensemble[iithin]),-1) - (2*mu)
    X = X[2L*lindgen(Lh)]

    status = acor(X,info,Lh,ithin=1,iburn=0,maxlag=maxlag,$
                    sigma=sigma,tau=tau,winmult=winmult)
    D     = 0.25*sigma*sigma*L
    tau   = D/C[0]
    sigma = sqrt(D/L)
  endelse
  return,status
end
;
; ***********************************************************
; simple acor test:
;   creates time series with nozero mean and autocovariance
;   function C(t) = C(0)*a(t) w/ a = 0.6.  
;   Follows J. Goodman's acorTest program.
;
; see: http://www.math.nyu.edu/faculty/goodman/software/acor/
;
; I get: 
;  p0: tau =        18.762740 & sigma =     0.0014326758       0
;
;  - 5/13 : created by KR
;
pro acorTest
  a    = 0.9d
  x    = 0d
  n    = 4000000L
  seed = 5L
  info = {nsamples:n,ndim:1,pname:'p0'}
  xs = dblarr(n)
  ; generate series 
  print,'generating series'
  for i=0L,n-1 do begin
    x = a*x + randomu(seed)
    xs[i] = x
  endfor
  ; calculate acor
  print,'calculating autocorrelation time'
  pl_acor,reform(xs,1,1,n),info,maxlag=25,pmulti=[1,1]
  print,'finished acorTest'
end
;
; ***********************************************************
; pl_acor:
;   wrapper to calculate autocorrelation time 
; 
;  - 5/13 : created by KR
;
pro pl_acor,ensemble,info,iburn=iburn,ithin=ithin,pmulti=pmulti,$
    maxlag=maxlag,mkfig=mkfig,winmult=winmult,report=report
  ; defaults
  if (not keyword_set(iburn))  then iburn  = 0
  if (not keyword_set(ithin))  then ithin  = 1
  if (not keyword_set(maxlag)) then maxlag = 10
  if (not keyword_set(winmult)) then winmult = 2    ; is 5 in acor
  if (not keyword_set(mkfig))  then mkfig  = 0  ; make figures
  if (not keyword_set(pmulti)) then pmulti = [2,2] 
  if (not keyword_set(report)) then report = 1
  if (n_elements(pmulti) eq 1) then begin
    foo = ceil(sqrt(info.ndim))
    pmulti = [foo,foo]
  endif

  sz = size(ensemble)
  nplots = pmulti[0]*pmulti[1]
  cnt = 0
  if (mkfig) then print,'plotting autocorrelation:' 
  for idim=0,info.ndim-1 do begin
    ; calculate
    N = info.nsamples
    s = reform(mean(ensemble[idim,*,*],dim=2),sz[3])
    if (iburn lt 0) then $
         autoburn,s,info,minburn=info.mburn,ithin=ithin,_iburn=_iburn $
    else _iburn = iburn
    status = acor(s,info,N,ithin=ithin,iburn=_iburn,maxlag=maxlag,$
                tau=tau,sigma=sigma,_C=C,winmult=winmult)
    ; report
    if (report) then print,info.pname[idim]+': tau = ',tau,' & sigma = ',sigma,status
    ; check status
    if (status ne 0) then cnt ++
    if (mkfig) then begin
      if ((idim mod nplots) eq 0) then !p.multi = [0,pmulti[0],pmulti[1]]
      ; plot
      plot,C,xtit='T',ytit='acor ('+info.pname[idim]+')'
      if ((idim+1) mod nplots eq 0) then stop
      ;if finite(tau,/nan) then begin
      if status ne 0 then begin
        xx = mean(ensemble[idim,*,iburn:-1],dim=2)
        openw,1,'out/debug.dat'
        for i=0,n_elements(xx)-1 do $
            printf,1,xx[i]
        close,1
        stop
      endif
    endif
  endfor
  ; report if not converged
  if (cnt gt 0) then $
        print,'WARNING: '+string(cnt,format='(I0)')+'/'+string(info.ndim,format='(I0)')+$
        ' parameters are not converged'
  !p.multi=0

  ; plotting autocovariance
  print,'finished pl_acor'
end

pro acorinit,e=e,i=i
  readcol,'ao',e
  i = {nsamples:n_elements(e),ndim:1,pname:'p0'}
  pl_acor,reform(e,1,1,n_elements(e)),i,maxlag=25,pmulti=[1,1]
end 
; ***********************************************************
pro pl_mkct
; set up colors (M. Cushing's color table)
  if !d.name eq 'X' then device, decomposed=0, pseudo=8
  red  =[  0,255,255,  0,  0,255,255,  0,255,127,127 ]
  green=[  0,255,  0,255,  0,255,  0,255,127,255,127 ]
  blue =[  0,255,  0,  0,255,  0,255,255,127,127,255 ]
  tvlct, red, green, blue
  loadct, 0, bottom=n_elements(red), /silent
end
