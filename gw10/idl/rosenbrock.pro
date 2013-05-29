function rosenbrock,x,_extra=_extra
  return,-1./20*(100.*(x[1] - x[0]^2)^2+(1 - x[0])^2)
end
