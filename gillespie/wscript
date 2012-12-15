# Build script for Gillespie Solver (draft)
# vim: syntax=python 
top = '.'
out = 'build'

# Header files which this module requires.
header_list = ['vector', 'map', 'numeric']

hppfiles = [
	'GillespieSolver.hpp', 'GillespieWorld.hpp', 'GillespieSimulator.hpp'
	]

cppfiles = [
	'GillespieSolver.cpp', 'GillespieWorld.cpp', 'GillespieSimulator.cpp', #'serialize.cpp', 
	]

def options(opt):
	opt.add_option('--enable_debug', action='store_true', default=False, help='debug')
	opt.load('compiler_cxx')

def configure(conf):
	conf.load('compiler_cxx')

	conf.check_cfg(package='gsl', uselib_store='gsl', atleat_version='1.13', args='--cflags --libs')
	conf.check_cfg(package='pficommon', uselib_store='pficommon', atleat_version='1.0.0', args='--cflags --libs')

	# Checking the existence of header files.
	for header in header_list:
		conf.check(header_name = header, features = 'c cprogram')

	# Save option flags.
	conf.env.enable_debug =  conf.options.enable_debug

	conf.env.append_unique(
		'CXXFLAGS', 
		['-Wall', '-g']
		)
	

def build(bld):
	# always build libgillespie.so or .dylib(mac)
	#bld.shlib(
	#	source = ['./GillespieSolver.cpp', './GillespieWorld.cpp', './serialize.cpp'],
	#	includes = ['.'],
	#	uselib = ['gsl', 'pficommon'],
	#	target = 'gillespie'
	#)

	bld.install_files(
		'${PREFIX}/ecell/bd', hppfiles)

	bld.shlib(
		source = cppfiles,
		includes = ['.', '..'],
		libpath = ['../build/core'],
		lib = ['gsl', 'gslcblas', 'm', 'ecell4-core'],
		target = 'ecell4-gillespie')


