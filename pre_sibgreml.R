dges = grm %*% zdir
iges = zsb %*% grm %*% t(zsb)
idges = zdir %*% grm %*% t(zsb) + 
  zsb %*% grm %*% t(zdir)
S = zsb + zdir
