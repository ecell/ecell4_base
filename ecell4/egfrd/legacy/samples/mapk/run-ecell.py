



#loadModel('model3a.eml')
#loadModel('model3.eml')
loadModel('model4.eml')

lkpp = createLoggerStub('Variable:/:Kpp:Value')
lkpp.create()

#lk = createLoggerStub('Variable:/:K:Value')
#lk.create()

run(1200)

from ecell.ECDDataFile import *

#ECDDataFile(lkpp.getData()).save('Kpp_ODE_0.ecd')
ECDDataFile(lkpp.getData()).save('Kpp2.ecd')
#ECDDataFile(lk.getData()).save('K.ecd')
