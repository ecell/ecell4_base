#! /usr/bin/env python
# encoding: utf-8


top = '.'
out = 'build'

subdirs = [
    'core',
    # 'bd',
    # 'ode',
    # 'gillespie',
    # 'vis'
    ]

def options(opt):
    opt.recurse(subdirs)

def configure(conf):
    conf.recurse(subdirs)

def build(bld):
    bld.recurse(subdirs)
