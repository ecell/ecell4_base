#! /usr/bin/env python
# encoding: utf-8


top = '.'
out = 'build'

subdirs = [
    'core',
    # 'ecell4-bd',
    # 'ecell4-ode',
    # 'ecell4-gillespie',
    # 'ecell4-vis'
    ]

def options(opt):
    opt.recurse(subdirs)

def configure(conf):
    conf.recurse(subdirs)

def build(bld):
    bld.recurse(subdirs)
