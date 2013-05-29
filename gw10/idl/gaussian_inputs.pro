; sampler settings
lnprobfn = 'gaussian'   ; name of model function
nwalkers = 150L         ; # of walkers
nsamples = 25000L       ; max # of samples/walker
acormult = 12L          ; min # of autocorrelation times
ndim	 = 50L          ; # of parameters
seed	 = 15L          ; random seed
ioduty   = 500L         ; cycles to save samples/report
cduty    = 2000L        ; cycles to check for convergence 
mburn    = 1000L        ; burn mininum
npars    = 'p'+string(findgen(ndim),format='(I0)')
mpars    = randomu(seed,ndim) 
spars    = 0.9 + dblarr(ndim) 

; simulation settings 
writeme   = 0   ; save to txt file
saveme    = 0   ; save to idl file  @ end
dir       = 'out/'
root      = 'debug'

; createan mu and sigma for multi-dimensional gaussian
means = randomu(seed,ndim)

; w/ random covariances
cov   = 0.5 - randomu(seed,ndim,ndim)
for i=0,ndim-1 do $
  for k=0,i-1 do $
    cov[i,k] = cov[k,i]

scov = cov
cov = cov##cov
icov = invert(cov,status,/double)
if (status ne 0) then stop,'ERROR in creating covariance matrix'

; package extras for function
extra = {mu:means,icov:icov}
info  = {pname:npars,ndim:ndim,nwalkers:nwalkers,nsamples:nsamples,$
         mburn:mburn}
