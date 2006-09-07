import gfrd

r = 5e-8
r0 = 3e-8

t = 1e-5
D = 1e-12
Sigma = 1e-8
kf = 1e-18


gf = gfrd.PlainPairGreensFunction( D, kf, Sigma )

#for i in range(1000):
#print gf.drawR( 0.9, r0, t )

    
for i in range(1000):
    print gf.drawTheta( 0.5, r, r0, t )
