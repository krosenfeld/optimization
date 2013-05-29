function gaussian,x,mu=mu,icov=icov
  diff = x-mu
  return, -0.5*diff##(icov##diff)
end
