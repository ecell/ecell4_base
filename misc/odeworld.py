from ecell4 import *
from ecell4.ode import ODEWorld
edge_lengths = Real3(1, 2, 3)
w2 = ODEWorld(edge_lengths)
m = NetworkModel()
w2.bind_to(m)
