#! /usr/bin/env python
# encoding: utf-8

from waflib.Tools import waf_unit_test
from waflib import Logs


top = '.'
out = 'build'

subdirs = [
    'tests'
    ]

hppfiles = [
    'BDWorld.hpp', 'BDSimulator.hpp', 'BDPropagator.hpp', 'BDSimulatorState.hpp',
    'functions3d.hpp'
    ]

cppfiles = [
    'BDSimulator.cpp', 'BDPropagator.cpp', 'functions3d.cpp'
    ]

def options(opt):
    opt.load('compiler_cxx waf_unit_test')

def configure(conf):
    conf.load('compiler_cxx waf_unit_test')
    conf.check_cxx(lib='gsl')
    conf.check_cxx(lib='gslcblas')
    conf.check_cxx(lib='m')

    conf.check_cxx(lib='ecell4-core')

    conf.recurse(subdirs)

def summary(bld):
    '''borrowed from waf demos/unit_test/wscript
    '''
    lst = getattr(bld, 'utest_results', [])
    if lst:
        total = len(lst)
        tfail = len([x for x in lst if x[1]])

    val = 100 * (total - tfail) / (1.0 * total)
    Logs.pprint('CYAN', 'test report %3.0f%% success' % val)

    Logs.pprint('CYAN', '  tests that fail %d/%d' % (tfail, total))
    for (f, code, out, err) in lst:
        if code:
            Logs.pprint('CYAN', '    %s' % f)
            Logs.pprint('RED', 'status: %r' % code)
            if out: Logs.pprint('RED', 'out: %r' % out)
            if err: Logs.pprint('RED', 'err: %r' % err)

def build(bld):
    bld.install_files(
        '${PREFIX}/ecell/bd', hppfiles)

    bld.shlib(
        source = cppfiles,
        includes = ['.'],
        lib = ['gsl', 'gslcblas', 'm', 'ecell4-core'],
        target = 'ecell4-bd')

    bld.recurse(subdirs)

    # bld.add_post_fun(waf_unit_test.summary)
    bld.add_post_fun(summary)
    bld.options.all_tests = True
