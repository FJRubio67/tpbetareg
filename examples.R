mu0 = 0.5
phi0 = 10
gamma0 = 0.75

mu1 = 0.875
phi1 = 10


tempf0 <- Vectorize(function(x) dbetamv(x, mu1, phi1, log = FALSE) )
tempf <- Vectorize(function(x) dtpbetamv(x, mu0, phi0, gamma0, log = FALSE) )

curve(tempf0,0,1,n=1000, lwd = 2, ylim = c(0,10))
curve(tempf,0,1,n=1000, lwd = 2, lty = 2, add =TRUE)



inte = Vectorize(function(x) tempf(x)*x)

integrate(inte,0,1)
