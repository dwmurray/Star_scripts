import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as pl 

def shu_derivs( x, y) :
	v = y[0]
	alpha = y[1] 
	dvdx = (alpha*(x-v) - 2e0/x)*(x-v)/((x-v)**2. - 1.)
	dalphadx = alpha*(alpha-2./x*(x-v))*(x-v)/((x-v)**2. - 1.)
	return [dvdx, dalphadx]

def shu_solution( x0, x1, v0, alpha0 ) : 
	integrator = si.ode( shu_derivs).set_integrator( "dopri5")
	integrator.set_initial_value( [v0,alpha0], x0)
	dx = 0.001*(x1-x0)
	vs = []
	alphas = [] 
	xs = []
	while( integrator.successful() and integrator.t > x1) : 
		result = integrator.integrate(integrator.t+dx)
		xs.append(integrator.t)
		vs.append( result[0])
		alphas.append( result[1])

	xs = np.array(xs)
	vs = np.array(vs)
	alphas = np.array(alphas)
	return xs, vs, alphas

pl.clf()
for A in np.arange( 2.001, 20., 1.) :
	x0 = 100.
	alpha0 = A/x0**2 - A*(A-2.)/(2.*x0**4) 
	v0 = -(A-2.)/x0 - ( 1. - A/6.)*(A-2.)/x0**3
	cs = 2.65e4
	x, v, alpha = shu_solution( x0, 0.01, v0, alpha0) 
	mdot = -cs**3./6.67e-8*x*x* alpha* v/2e33*3.15e7
	if A == 3.001:
		pl.loglog( x, mdot, label="A={0}".format(A))
pl.xlim( 0, 100.0) 
pl.legend(loc="best") 
pl.savefig("test.pdf") 


