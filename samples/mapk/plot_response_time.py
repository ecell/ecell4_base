#!/usr/bin/env python


# ti=1e-6
# python plot_response_time.py 09-3/data/mapk3_1e-15_0.03125_fixed_1e-6_normal_ALL_tc.dat 09-3/data/mapk3_1e-15_0.0625_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_1_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-6_normal_ALL_tc.dat 


# ti=1e-2
# python plot_response_time.py 09-3/data/mapk3_1e-15_0.03125_fixed_1e-2_normal_ALL_tc.dat 09-3/data/mapk3_1e-15_0.0625_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_0.25_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_1_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-2_normal_ALL_tc.dat 


#ODE_file_m6 = '/home/shafi/wrk/epdp/samples/mapk/Kpp_ODE_1e-6.ecd'
#ODE_file_m2 = '/home/shafi/wrk/epdp/samples/mapk/Kpp_ODE_1e-2.ecd'


import sys
import os
import glob

import numpy
import scipy.io

from matplotlib.pylab import *

def load_header( filename ):
    file = open( filename )
    header = []
    for line in file.readlines():
        if line[0:2] == '#@':
            hline = line[2:].lstrip()
            header.append( hline )

    return header

def resample( x, y, newx ):

    indices = numpy.searchsorted( x, newx )

    indices = indices.clip( 0, len(y) - 1 )
    #print indices, len(y)

    return y.take( indices )
    

def add_columns( data, ycolumns ):

    y = numpy.array([ data[:,col] for col in ycolumns ]) 

    y = y.sum(0)

    return y


def load_data( filename ):
    ycolumns = [1,]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    header = load_header( filename )
    #print header
    for l in header:
        exec( l )

    #data = numpy.loadtxt( filename )
    data = load( filename )
    x = data[:,0]
    y = add_columns( data, ycolumns )

    return x, y


def plot_file( filename, lp='-' ):

    x, y = load_data( filename )

    #plot_theory( N_K, N_P, Keq, x[-1] )
    plot( x, y, lp )

    #psd( y )
    #ylim( 1, 5e4 )


from scipy.optimize import leastsq

def fitter( a, b, x ):
    return a * (1-numpy.exp(-b*x))

def residuals(p, y, x):
    a,b=p
    return y - fitter( a, b, x )

def t_m(a,b):
    return log(2) / b

p0 = [50,10]

def response_time( filelist, end ):

    start = 0.
    interval = (end-start) / 1000.
    rx = numpy.mgrid[start:end:interval]

    data = []

    assert filelist

    tm=[]
    for filename in filelist:
        #print 'file ', filename
        x, y = load_data( filename )
        #print x,y
        ry = resample( x, y, rx )
        data.append( ry )

        res = leastsq(residuals, p0, args=(y, x),full_output=1)
        print res[0], numpy.diag(res[1])

        mry = numpy.array( data ).mean( 0 )

        ly = fitter(res[0][0],res[0][1],rx)
        #plot( rx, ly )

        tm.append(t_m(res[0][0], res[0][1]))


    #plot( rx, mry, label=l )
    tm = numpy.array(tm)

    return tm



def rt_pattern( pattern, end, x ):
    globpattern = pattern.replace('ALL','*')
    
    l = os.path.basename( os.path.splitext( pattern )[0] )
    print 'pattern ', l

    filelist = glob.glob( globpattern )

    tm = response_time( filelist, end )

    return tm.mean(), tm.std()/math.sqrt(len(tm))



if __name__ == '__main__':


    import glob
    import os

    xmax = 120

    dir = '09-4/data'

    model = 'mapk3'
    V_str = '1e-15'


    lines=[]

    #for ti_str in ['0','1e-6','1e-4','1e-2']:
    for ti_str in ['1e-2']:#['1e-6','1e-2']:

        x = []
        y = []

        for D_str in ['0.25', '1', '2','4']:##['0.03125','0.0625','0.125', '0.25','0.5','1','2','4']:#
            
            globpattern = \
                '_'.join( ( model, V_str, D_str, 'fixed', ti_str, 
                            'normal', '*' ) ) + '_tc.dat'

            filelist = glob.glob( dir + os.sep + globpattern )

            if not filelist:
                continue

            print globpattern
            print filelist
            
            filelist=filelist[-50:]


            ti = float(ti_str)
            D = float(D_str)

            tm = response_time( filelist, xmax )
            mean, std_err = tm.mean(), tm.std()/math.sqrt(len(tm))

            errorbar( D, mean, yerr=std_err, fmt='k+' )

            x.append( D )
            y.append( mean )

        line = plot( x, y )
        lines.append(line)

    ls='k--'
    for ti_str in [ '1e-6', '1e-2' ]:
        x = []
        y = []
        globpattern = 'Kpp_ODE_*_%s.ecd' % ti_str
        filelist = glob.glob( globpattern )

        for file in filelist:
            
            ODE_file = file
            D_str = file.split('_')[2]
            otm = response_time( [ODE_file], xmax )[0]
            x.append( float( D_str ) )
            y.append( otm )
            print 'otm', otm

        x,y = numpy.array(x), numpy.array(y)
        args = x.argsort()
        x = x.take(args)
        y = y.take(args)

        print x,y
        line = plot( x, y, ls )
        ls = 'k-'
        lines.append(line)


#     otm = response_time( [ODE_file_m2], xmax )[0]
#     line = plot( [1e-18,10], [otm,otm], 'k-' )
#     lines.append(line)
#     print otm

    #plot_file(ODE_file, 'k-' )



    #xticks( size=20 )
    #yticks( size=20 )

    #xlabel( r'$t_{\rm rel} {\rm [s]}$', size=22 )
    xlabel( r'$D {\rm [\mu m^2 / s]}$', size=22 )
    ylabel( r'Response time [s]', size=22 )

    xscale( 'log' )

    xlim( 0.02, 5 )
    ylim( 0, 15 )


    xticks([0.1,1,10],['0.1','1','10'], size=20 )
    yticks([0,5,10,15],['0','5','10','15'], size=20 )

#     leg =legend( lines, (r'$t_{\rm rel} = 1 {\rm \mu s}$',
#                          r'$t_{\rm rel} = 10 {\rm m s}$',
#                          r'$t_{\rm rel} = 1 {\rm \mu s} {\rm (ODE)}$',
#                          r'$t_{\rm rel} = 10 {\rm m s} {\rm (ODE)}$',
#                   ),
#                   loc=1,
#                   shadow=True,
#                   pad=0.05
#                   )
#     for l in leg.get_lines():
#         l.set_linewidth(1.5)  # the legend line width


#title( figtitle )

#savefig( 'figs/' + figtitle + '.png', dpi=80 )

    show()
