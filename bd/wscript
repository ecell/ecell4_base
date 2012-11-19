#! /usr/bin/env python
# encoding: utf-8

from waflib.Tools import waf_unit_test


top = '.'
out = 'build'

subdirs = [
    'tests'
    ]

hppfiles = [
    'BDWorld.hpp'
    ]

cppfiles = []

def options(opt):
    opt.load('compiler_cxx waf_unit_test')

def configure(conf):
    conf.load('compiler_cxx waf_unit_test')
    conf.check_cxx(lib='gsl')
    conf.check_cxx(lib='gslcblas')
    conf.check_cxx(lib='m')

    conf.check_cxx(lib='ecell4-core')

    conf.check_cfg(package='cppunit', args='--cflags --libs', mandatory=True)
    if 'dl' not in conf.env.LIB_CPPUNIT:
        l = conf.check(lib='dl', uselib_store='CPPUNIT')

    conf.recurse(subdirs)

def build(bld):
    bld.install_files(
        '${PREFIX}/ecell/bd', hppfiles)

    # bld.shlib(
    #     source = cppfiles,
    #     includes = ['.'],
    #     lib = ['gsl', 'gslcblas', 'm'],
    #     target = 'ecell4-core')

    bld.recurse(subdirs)

    bld.add_post_fun(waf_unit_test.summary)
    bld.options.all_tests = True
