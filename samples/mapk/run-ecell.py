
loadModel( 'model2.eml' )


l = createLoggerStub( 'Variable:/:Kpp:Value' )
l.create()

run( 1000 )


from ecell.ECDDataFile import *

ECDDataFile( l.getData() ).save( 'Kpp2.ecd' )
