

import gfrd
import math

import numpy

t = 1e-4
D = 1e-12
a = 2e-9


gf = gfrd.FirstPassageGreensFunction( D, a )

print gf.p_survival( 1e-6 )
for i in range(1000000):
    gf.p_survival( 1e-6 )


