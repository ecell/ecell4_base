
import os
import re
#import logging
import numpy

import logging
import h5py # added by sakurai@advancesoft.jp
from egfrd import Single, Pair, Multi # added by sakurai@advancesoft.jp

__all__ = [
    'FixedIntervalInterrupter',
    'Logger',
    'HDF5Logger',
    ]

INF = numpy.inf

log = logging.getLogger('ecell')

PARTICLES_SCHEMA = \
    [
        ('id', 'u8', ),
        ('species_id', 'u8', ),
        ('position', 'f8', (3, ))
        ]

SHELLS_SCHEMA = \
    [
        ('id', 'u8', ),
        ('radius', 'f8'),
        ('position', 'f8', (3, )),
        ]

SPECIES_SCHEMA = \
    [
        ('id', 'u8', ),
        ('name', 'S32', ),
        ('radius', 'f8', ),
        ('D', 'f8'), # diffusion coefficient
        ]

SHELL_PARTICLE_ASSOCIATION_SCHEMA = \
    [
        ('shell_id', 'u8'),
        ('particle_id', 'u8'),
        ]

DOMAINS_SCHEMA = \
    [
        ('id', 'u8', ),
        ('kind', 'u4', ),
        ]

DOMAIN_SHELL_ASSOCIATION_SCHEMA = \
    [
        ('shell_id', 'u8', ),
        ('domain_id', 'u8', ),
        ]

class FixedIntervalInterrupter(object):
    def __init__(self, sim, interval, callback):
        self.sim = sim
        self.interval = interval
        self.callback = callback
        self.last_time = 0.
        self.next_time = interval

    def step(self):
        self.sim.step()
        if self.next_time <= self.sim.t + self.sim.dt:
            self.sim.stop(self.next_time)
            self.callback(self.sim, self.next_time)
            self.last_time = self.sim.t
            self.next_time += self.interval


class HDF5Logger(object):
    def __init__(self, logname, directory='data', split=False):
        self.logname = logname
        self.directory = directory
        self.split = split
        self.file_counter = 0
        self.hdf5_file = None

    def new_hdf5_file(self, sim):
        if self.split:
            if self.hdf5_file is not None:
                self.hdf5_file.close()

            if not os.path.exists(self.directory):
                os.mkdir(self.directory)

            hdf5_filename = '%s_%04d.hdf5' % (self.logname, self.file_counter)
            self.file_counter += 1
        else:
            if self.hdf5_file is not None:
                return
            hdf5_filename = '%s.hdf5' % self.logname
        hdf5_path = os.path.join(self.directory, hdf5_filename)
        # HDF5 file must be removed before log_particles
        if os.path.exists(hdf5_path):
            os.remove(hdf5_path)
        self.hdf5_file = h5py.File(hdf5_path)

        self.create_data_group(sim)
        self.write_species(sim)

    def create_data_group(self, sim):
        data_group = self.hdf5_file.create_group('data')
        data_group.attrs['world_size'] = sim.world.world_size
        return data_group

    def write_species(self, sim):
        num_species = len(sim.world.species)

        species_dset = self.hdf5_file.create_dataset('species', (num_species, ), SPECIES_SCHEMA)
        count = 0
        for species in sim.world.species:
            species_dset[count] = (species.id.serial,
                                   str(species.id),
                                   species.radius,
                                   species.D)
            count += 1

    def write_particles(self, sim):
        "This function was created by sakurai@advancesoft.jp"
        data_group = self.hdf5_file['data']

        time_group = data_group.require_group(unicode(sim.t))
        time_group.attrs['t'] = sim.t

        # Create particles dataset on the time group

        num_particles = 0
        for species in sim.world.species:
            pid_list = sim.world.get_particle_ids(species.id)
            num_particles += len(pid_list)

        x = numpy.zeros((num_particles, ),
                        dtype = numpy.dtype(PARTICLES_SCHEMA))

        count = 0
        for species in sim.world.species:
            pid_set = sim.world.get_particle_ids(species.id)
            for pid in pid_set:
                pid, particle = sim.world.get_particle(pid)
                x['id'][count] = pid.serial
                x['species_id'][count] = species.id.serial
                x['position'][count] = particle.position
                count += 1

        dummy = time_group.create_dataset('particles', data = x)

    def write_domains(self, sim):
        if self.hdf5_file is None:
            return

        data_group = self.hdf5_file['data']        

        # Require time group
        time_group = data_group.require_group(unicode(sim.t))
        time_group.attrs['t'] = sim.t

        # Create shell dataset on the time group

        num_shells = 0
        for domain in sim.domains.itervalues():
            num_shells += len(domain.shell_list)

        x = numpy.zeros((num_shells, ), dtype = numpy.dtype(SHELLS_SCHEMA))

        count = 0
        for did, domain in sim.domains.iteritems():
            shell_list = domain.shell_list
            for shell_id, shell in shell_list:
                x['id'][count] = shell_id.serial
                x['radius'][count] = shell.shape.radius
                x['position'][count] = shell.shape.position
                count += 1

        if len(x) == 0:
            return

        dummy = time_group.create_dataset('shells', data = x)

        # Create shell particle association dataset on the time group

        num_assocs = 0
        for domain in sim.domains.itervalues():
            if isinstance(domain, Single):
                num_assocs += len(domain.shell_list)
            elif isinstance(domain, Pair):
                num_assocs += 2 * len(domain.shell_list)
            elif isinstance(domain, Multi):
                assert getattr(domain, 'pid_shell_id_map', None), 'Cannot access pid_shell_id_map'
                num_assocs += len(domain.pid_shell_id_map)

        dtype_obj = numpy.dtype(SHELL_PARTICLE_ASSOCIATION_SCHEMA)
        x = numpy.zeros((num_assocs, ), dtype = dtype_obj)

        count = 0
        for did, domain in sim.domains.iteritems():

            if(isinstance(domain, Single) or
               isinstance(domain, Pair)):

                pid_particle_pair_list = []
                if isinstance(domain, Single):
                    pid_particle_pair_list = [domain.pid_particle_pair]
                elif isinstance(domain, Pair):
                    pid_particle_pair_list = [domain.single1.pid_particle_pair,
                                              domain.single2.pid_particle_pair]

                for pid, particle in pid_particle_pair_list:
                    for shell_id, shell in domain.shell_list:
                        x['shell_id'][count] = shell_id.serial
                        x['particle_id'][count] = pid.serial
                        count += 1

            else: # for Multi
                assert getattr(domain, 'pid_shell_id_map', None), 'Cannot access pid_shell_id_map'
                for pid, shell_id in domain.pid_shell_id_map.iteritems():
                    x['shell_id'][count] = shell_id.serial
                    x['particle_id'][count] = pid.serial
                    count += 1

        dummy = time_group.create_dataset('shell_particle_association', data = x)

        # Create domain_shell_association dataset on the time group
        dtype_obj = numpy.dtype(DOMAIN_SHELL_ASSOCIATION_SCHEMA)
        x = numpy.zeros((num_shells, ), dtype = dtype_obj)

        count = 0
        for did, domain in sim.domains.iteritems():
            shell_list = domain.shell_list
            for shell_id, shell in shell_list:
                x['shell_id'][count] = shell_id.serial
                x['domain_id'][count] = did.serial
                count += 1

        dummy = time_group.create_dataset('domain_shell_association', data = x)

        # Create domain dataset on the time group
        num_domains = len(sim.domains)

        dtype_obj = numpy.dtype(DOMAINS_SCHEMA)
        x = numpy.zeros((num_domains, ), dtype = dtype_obj)

        count = 0
        for did, domain in sim.domains.iteritems():
            x['id'][count] = did.serial
            if isinstance(domain, Single):
                x['kind'][count] = 1
            elif isinstance(domain, Pair):
                x['kind'][count] = 2
            else: # must be Multi
                x['kind'][count] = 3
            count += 1

        dummy = time_group.create_dataset('domains', data = x)

    def log(self, sim, time):
        self.new_hdf5_file(sim)
        self.write_particles(sim)

    def start(self, sim):
        self.new_hdf5_file(sim)
        self.write_particles(sim)


class Logger(object):
    def __init__(self, logname='log', directory='data', comment=''):
        self.logname = logname
        self.file_counter = 0
        self.directory = directory
        self.comment = comment
        self.timecourse_file = None

    def prepare_timecourse_file(self, sim):
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        timecourse_filename = '%s_tc.dat' % self.logname
        self.timecourse_file = open(
            os.path.join(self.directory, timecourse_filename), 'w')
        self.write_timecourse_comment(self.comment)

        species_name_list = '\'' + \
            "\', \'".join(sim.world.model.get_species_type_by_id(i.id)['name']
                          for i in sim.world.species) + '\''
        columns = '[\'t\', ' + species_name_list + ']'
        self.write_timecourse_comment('@ columns= ' + columns)

    def write_timecourse_comment(self, s):
        self.timecourse_file.write('#' + s + '\n')

    def write_timecourse(self, sim):
        data = []
        self.timecourse_file.write('%g' % sim.t + '\t')
        self.timecourse_file.write('\t'.join(
            str(len(sim.world.get_particle_ids(i.id))) \
            for i in sim.world.species) + '\n')
        self.timecourse_file.flush()

    # this method will be deprecated.
    def write_particles(self, sim):
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        filename = '%s_%04d.dat' % (self.logname, self.file_counter)
        file = open(os.path.join(self.directory, filename), 'w')

        file.write('#@ name = \'%s\'\n' % str(self.logname))
        file.write('#@ count = %d\n' % int(self.file_counter))
        file.write('#@ t = %s\n' % '%g' % sim.t)
        file.write('#@ world_size = %f\n' % float(sim.world.world_size))
        file.write('#--------\n')

        for species in sim.world.species:
            pid_list = sim.world.get_particle_ids(species.id)
            for pid in pid_list:
                pid, particle = sim.world.get_particle(pid)
                st = sim.world.model.get_species_type_by_id(species.id)
                file.write('%s\t%20.14g %20.14g %20.14g %.15g\n' %
                           (st['name'], particle.position[0], 
                            particle.position[1], particle.position[2],
                            species.radius))

            file.write('#\n')

        file.close()

        self.file_counter += 1

    def log(self, sim, time):
        self.write_timecourse(sim)

    def start(self, sim):
        self.prepare_timecourse_file(sim)
        self.write_timecourse(sim)
