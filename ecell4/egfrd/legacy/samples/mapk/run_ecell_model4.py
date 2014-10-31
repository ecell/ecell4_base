loadModel(MODEL_FILE)

Kpp = createEntityStub('Variable:/:Kpp') 
K = createEntityStub('Variable:/:K') 
KK = createEntityStub('Variable:/:KK') 
P = createEntityStub('Variable:/:P') 

Kpp['Value'] = N_KPP
K['Value'] = N_K
KK['Value'] = N_KK
P['Value'] = N_P

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

run(1000000)
print Kpp['Value'] / (N_KPP + N_K)

