import gfrd
import math

r = 1e-7
r0 = 5e-8

t = 1e-4
D = 1e-12
Sigma = 1e-8
kf = 1e-18


gf = gfrd.PlainPairGreensFunction( D, kf, Sigma )

#print gf.drawTime( 1e-10, 1e-8, 1.0 )
#for i in range(10000):
#    gf.drawTime( 0.2, 1e-8, 1.0 )
#    gf.drawTime( 0.5, 1e-8, 1.0 )
#    gf.drawTime( 0.8, 1e-8, 1.0 )


#for i in range(1000):
#    print gf.drawR( 0.9, r0, t )

    
for i in range(100):
    print gf.drawTheta( 0.5, r, r0, t )
