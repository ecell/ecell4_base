
loadModel( 'model3.eml' )


lkpp = createLoggerStub( 'Variable:/:Kpp:Value' )
lkpp.create()

lk = createLoggerStub( 'Variable:/:K:Value' )
lk.create()

run( 500 )


from ecell.ECDDataFile import *

ECDDataFile( lkpp.getData() ).save( 'Kpp.ecd' )
ECDDataFile( lk.getData() ).save( 'K.ecd' )
