import _gfrd
import math

import numpy

#print gfrd.distanceSq( numpy.array( (1.,2.,3.) ), numpy.array( (4.,5.,6.,) ) )


t = 1e-4
D = 1e-12
Sigma = 1e-8
kf = 1e-18

r = 1e-7
r0 = 5e-8
a = 2e-8



def test_alpha_survival_n():

    gf = _gfrd.FirstPassagePairGreensFunction( D, kf, Sigma )
    
    lower = 1e-8
    upper = 1e10
    for i in range(10):
        alpha = gf.alpha_survival_n( a, i, lower, upper )
        #lower = alpha
        print 'alpha, f_alpha_survival: ', i, alpha, gf.f_alpha_survival( alpha, a )


#print gf.drawTime( 1e-10, 1e-8, 1.0 )
#for i in range(10000):
#    gf.drawTime( 0.2, 1e-8, 1.0 )
#    gf.drawTime( 0.5, 1e-8, 1.0 )
#    gf.drawTime( 0.8, 1e-8, 1.0 )


#for i in range(1000):
#    print gf.drawR( 0.9, r0, t )


test_alpha_survival_n()





