"""Calculate dN/dm and N(>m) for a user defined powerlaw slope alpha;
takes two arguments, alpha and number_of_stars
start with m = m(t) = m_0 (t/t_0)^\alpha, so t = t_0 (m/m_0)^{1/\alpha}.
Then P(m)dm = P(t)dt = dt from which P(m) = dt/dm = (1/\alpha) (t_0/m_0)
(m/m_0)^{(1/\alpha) - 1}
"""

import numpy as n
import pylab as p
#import Norms_macros as sm
import sys

alpha = float(sys.argv[1])
number_of_stars = float(sys.argv[2])
while True:
    try:
        bins = float(sys.argv[3])
        break
    except ValueError:
        print "Need three arguments, alpha, number_of_stars, and bins"


#Powerlaw from mass=minimum_mass to very large masses
random_numbers = n.random.rand(number_of_stars) #uniform between 0 and 1
mass_max = 0.5  # minimum stellar mass
mass_random = mass_max * ( random_numbers)**(1./alpha)
#

log_mass = n.log10(mass_random)

Key = n.argsort(mass_random)
log_mass_sort = log_mass[Key]
log_mass_sort_reverse = log_mass_sort[::-1]
#number_mass = mass_sort.size - n.arange(0,mass_sort.size)

y_axis = n.histogram(log_mass,bins)
Number_of_stars_per_bin = y_axis[0]
x_bins = y_axis[1]
x_mid = (x_bins[1:] + x_bins[:-1])/2.

#Fit a straight line to the histogram above the midpoint of the curve

x_up = x_mid[bins/2:]
y_up = Number_of_stars_per_bin[bins/2:]

x = x_up
y = n.log10(y_up)
A = n.vstack([x, n.ones(len(x))]).T
m, c = n.linalg.lstsq(A,y)[0]
print('slope, intercept=',m, c)

y_fit = 10.**(m*x +c)


p.subplot(211)
p.hist(log_mass, bins, histtype='step', color='black', log=True)
p.plot(x_mid, Number_of_stars_per_bin, 'bo')
p.plot(x_up, y_up, 'ro')
p.plot(x_up, y_fit, 'r')
p.xlabel(r'$\log_{10}(m)$', size=15)
p.ylabel(r'$\log N_{*}$', size=20)

p.subplot(212)
x2 = n.arange(log_mass_sort_reverse.size) +1
log_x2 = n.log10(x2)
y2 = log_mass_sort_reverse
A2 = n.vstack([log_x2, n.ones(len(x2))]).T
m2, c2 = n.linalg.lstsq(A2,y2)[0]
print('Cumulative slope=',1./m2)

y_fit2 = m2*log_x2 +c2
p.plot(log_mass_sort_reverse, n.log10(x2),  color='black')
p.plot(y_fit2, log_x2, 'b')
p.xlabel(r'$\log_{10}(m)$', size=15)
p.ylabel(r'$\log N(>m)$', size=20)

p.savefig('dn_dm_histogram.png')
p.show()

