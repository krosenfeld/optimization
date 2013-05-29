; sampler settings
lnprobfn = 'rosenbrock'   ; name of model function
nwalkers = 150L         ; # of walkers
nsamples = 20000L       ; max # of samples/walker
acormult = 12L          ; # of autocorrelation times
ndim	 = 2L           ; # of parameters
seed	 = 15L          ; random seed
ioduty   = 100L         ; cycles to save samples/report
cduty    = 500L         ; cycles to check for convergence 
mburn    = 200L        ; burn mininum
npars    = 'p'+string(findgen(ndim),format='(I0)')
mpars    = randomu(seed,ndim) 
spars    = 0.9 + dblarr(ndim) 

; simulation settings 
writeme   = 0   ; save to txt file
saveme    = 0   ; save to idl file  @ end
dir       = 'out/'
root      = 'debug'

; package extras for function
extra = {}
info  = {pname:npars,ndim:ndim,nwalkers:nwalkers,nsamples:nsamples,$
         mburn:mburn}

;===================REPORT======================
;sampler converged in 9500 steps
;sampler took:     0.25 minutes
;accept ratio:     0.260491

