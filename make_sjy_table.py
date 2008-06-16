import numpy
import scipy.special as special



def maxz( n ):

    z = n * n + n + 1 * 2e4

    if z >= 1000:
        z = max( 1000, n * n )

    z *= 1.01

    return z


def jnyn( n, resolution ):

    delta = numpy.pi / resolution
    zTable = numpy.mgrid[delta:maxz(n):delta]

    jTable = numpy.zeros( ( len( zTable ), n+1 ) )
    yTable = numpy.zeros( ( len( zTable ), n+1 ) )
    for i, z in enumerate( zTable ):
        jTable[i], yTable[i], _, _ = special.sph_jnyn( n, z )

    jTable = jTable.transpose()
    yTable = yTable.transpose()
    return jTable, yTable

def make_table( func, n, z0, tol ):

    max_z = maxz( n )
    z = z0

    dz = numpy.pi / 100

    j, jp = func( n, z )
    jp_prev = jp[n]

    zTable = numpy.array( [z,] )
    yTable = numpy.array( [j,] )

    z += dz


    while z < max_z:
        j, jp = func( n, z )
        abs_jpp_norm = abs( jp[n] - jp_prev ) * dz

        if abs_jpp_norm > tol:
            dz *= .5
            continue 

        zTable = numpy.append( zTable, z )
        yTable = numpy.append( yTable, j[n] )

        if abs_jpp_norm < tol / 2:
            dz *= 2

        z += dz
        jp_prev = jp[n]

        #print z, table[-1]

    return zTable, yTable


def writeHeadArray( file, name, N ):

    file.write( 'static const double* %s[%d + 1] =\n{\n' % ( name, N ) ) 
    #file.write( 'boost::array<const double**, %d + 1>%s =\n{\n' % ( N, name ) ) 

    for n in range( N+1 ):
        file.write( '    &sj_table%d[0][0],\n' % n )

    file.write( '};' )


def writeArray( file, name, table1, table2 ):

    #    head_template = '''
    #static const double %s[2][%d + 1] =
    #{\n'''

    head_template = '''
    static const double %s[2][%d + 1] =
    {\n'''

    array_template = '''{\n%s\n}'''
    number_template = '''        %.18g,\n'''
    foot_template = '''};\n'''

    file.write( head_template % ( name, len(table1) ) )

    file.write( '    {\n' )
    for n in table1:
        file.write( number_template % n )
    file.write( '    },\n' )

    file.write( '    {\n' )
    for n in table2:
        file.write( number_template % n )
    file.write( '    }\n' )

    file.write( foot_template )


if __name__ == '__main__':

    import sys

    filename = sys.argv[1]

    file = open( filename, 'w' )

    maxn = 1
    resolution = 100
    tolerance = 1e-4

    #jTable, yTable = jnyn( maxn, resolution )

    for n in range( maxn+1 ):
        print n

        zTable, jTable = make_table( special.sph_jn, n, n+1, tolerance )
        writeArray( file, 'sj_table%d' % n, zTable, jTable )
        #writeArray( file, 'name', jTable )

        zTable, yTable = make_table( special.sph_yn, n, n+1, tolerance )
        writeArray( file, 'sy_table%d' % n, zTable, yTable )
        #writeArray( file, 'name', yTable )

    writeHeadArray( file, 'sj_table', maxn )

    file.write( '\n' )

    file.close()

