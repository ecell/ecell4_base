


loadModel(MODEL_FILE)

Kpp = createEntityStub('Variable:/:Kpp') 
K = createEntityStub('Variable:/:K') 
KK = createEntityStub('Variable:/:KK') 
P = createEntityStub('Variable:/:P') 

Kpp['Value'] = N_KPP
K['Value'] = N_K
KK['Value'] = N_KK
P['Value'] = N_P

sigma = 5e-9
D_tot = 1e-12 * float(D) * 2
#D_ref = 2e-12

kD = 4 * 3.14 * sigma * D_tot
k1 = 4.514e-20
k2 = 1.359
k4 = 9.20e-20
k5 = 1.73266

Keq12 = k1 / k2
Keq45 = k4 / k5

N_A = 6e23

k1_net = (1./((1./k1) + (1./kD)))*1000*N_A
k2_net = (1./((1./k2) + (Keq12/kD)))

k4_net = (1./((1./k4) + (1./kD)))*1000*N_A
k5_net = (1./((1./k5) + (Keq45/kD)))

message(k1_net)
message(k2_net)
message(k4_net)
message(k5_net)

R1 = createEntityStub('Process:/:R1') 
R2 = createEntityStub('Process:/:R2') 
R4 = createEntityStub('Process:/:R4') 
R5 = createEntityStub('Process:/:R5') 
R7 = createEntityStub('Process:/:R7') 
R8 = createEntityStub('Process:/:R8') 
R10 = createEntityStub('Process:/:R10') 
R11 = createEntityStub('Process:/:R11') 

R1['k'] = k1_net
R2['k'] = k2_net
R4['k'] = k4_net
R5['k'] = k5_net
R7['k'] = k1_net
R8['k'] = k2_net
R10['k'] = k4_net
R11['k'] = k5_net


if MODEL_FILE == 'model4-0.eml':
    pass
else:
    try:
        R13 = createEntityStub('Process:/:R13') 
        R14 = createEntityStub('Process:/:R14') 
        R13['k'] = KI
        R14['k'] = KI
        #print KI
    except:
        # processive model doesn't have R13, 14
        pass


lkpp = createLoggerStub('Variable:/:Kpp:Value')
lkpp.create()

run(DURATION)

from ecell.ECDDataFile import *
message(OUTFILE)
ECDDataFile(lkpp.getData()).save(OUTFILE)
message(lkpp.getData())
