#! /usr/bin/env python
# encoding: utf-8

from waflib.Tools import waf_unit_test


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

def build(bld):
    bld.install_files(
        '${PREFIX}/ecell/bd', hppfiles)

    bld.shlib(
        source = cppfiles,
        includes = ['.'],
        lib = ['gsl', 'gslcblas', 'm', 'ecell4-core'],
        target = 'ecell4-bd')

    bld.recurse(subdirs)

    bld.add_post_fun(waf_unit_test.summary)
    bld.options.all_tests = True
