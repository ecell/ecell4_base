
import os
import re
#import logging
import numpy

import logging
import h5py # added by sakurai@advancesoft.jp
from egfrd import Single, Pair, Multi # added by sakurai@advancesoft.jp

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


class Logger(object):
    def __init__(self, sim, logname = 'log', directory = 'data',
                  comment = '', hdf5=False):
        self.sim = sim
        self.logname = logname
        self.file_counter = 0
        self.directory = directory
        try:
            os.mkdir(directory)
        except:
            pass

        self.particle_out_interval = INF

        self.last_time = 0.0
        self.next_time = INF

        self.particle_out_pattern = re.compile('')
        self.prepare_timecourse_file(comment)
        self.write_timecourse()

        # Added by sakurai@advancesoft.jp
        if hdf5:
            HDF5_filename = self.logname + '.hdf5'
            HDF5_path = self.directory + os.sep + HDF5_filename
            # HDF5 file must be removed before log_particles
            if os.path.exists(HDF5_path):
                os.remove(HDF5_path)
            self.HDF5_file = h5py.File(HDF5_path)
        else:
            self.HDF5_file = None

    def __del__(self):
        pass
        #if self.HDF5_file is not None:
        #    self.HDF5_file.close()

    def set_interval(self, interval):
        self.interval = interval

    def set_particle_out_pattern(self, pattern):
        self.particle_out_pattern = re.compile(pattern)

    def get_particle_out_pattern(self):
        return self.particle_out_pattern.pattern

    def set_particle_out_interval(self, interval):
        self.particle_out_interval = interval
        self.last_time = self.sim.t
        self.next_time = self.last_time + self.particle_out_interval

    def prepare_timecourse_file(self, comment):
        self.timecourse_filename = self.logname + '_tc' + '.dat'
        self.timecourse_file = open(self.directory + os.sep + \
                                    self.timecourse_filename, 'w')
        self.write_timecourse_comment(comment)

        species_name_list = '\'' + \
            "\', \'".join(str(i) for i in self.sim.world.species) + '\''
        columns = '[\'t\', ' + species_name_list + ']'
        self.write_timecourse_comment('@ columns= ' + columns)

    def write_timecourse_comment(self, s):
        self.timecourse_file.write('#' + s + '\n')

    def write_timecourse(self):
        data = []
        self.timecourse_file.write('%g' % self.sim.t + '\t')
        self.timecourse_file.write('\t'.join(
            str(len(self.sim.get_particle_pool(i.id))) \
            for i in self.sim.get_species()) + '\n')
        self.timecourse_file.flush()

    def write_particles(self):
        filename = self.logname + '_' + \
            str(self.file_counter).zfill(4) + '.dat'

        file = open(self.directory + os.sep + filename, 'w')

        file.write('#@ name = \'%s\'\n' % str(self.logname))
        file.write('#@ count = %d\n' % int(self.file_counter))
        file.write('#@ t = %s\n' % '%g' % self.sim.t)
        file.write('#@ world_size = %f\n' % float(self.sim.world.world_size))
        file.write('#--------\n')

        for species in self.sim.world.species:
            pid_list = self.sim.particle_pool[species.id]
            for pid in pid_list:
                pid, particle = self.sim.world.get_particle(pid)
                species = self.sim.world.get_species(species.id)
                file.write('%s\t%20.14g %20.14g %20.14g %.15g\n' %
                            (species.id, particle.position[0], particle.position[1], particle.position[2], species.radius))

            file.write('#\n')

        file.close()

        self.file_counter += 1

    def create_data_group(self):
        if self.HDF5_file is None:
            return
        data_group = self.HDF5_file.create_group('data')
        data_group.attrs['world_size'] = self.sim.world.world_size
        return data_group

    def write_species_byHDF5(self):
        if self.HDF5_file is None:
            return

        # Create species dataset on top level of HDF5 hierarchy

        num_species = len(self.sim.world.species)

        species_dset = self.HDF5_file.create_dataset('species', (num_species, ), SPECIES_SCHEMA)
        count = 0
        for species in self.sim.world.species:
            species_dset[count] = (species.id.serial,
                                   str(species.id),
                                   species.radius,
                                   species.D)
            count += 1

    def write_particles_byHDF5(self):
        "This function was created by sakurai@advancesoft.jp"
        if self.HDF5_file is None:
            return

        data_group = self.HDF5_file['data']

        group_name = unicode(self.next_time)
        if group_name in data_group:
            return
        time_group = data_group.create_group(group_name)
        time_group.attrs['t'] = self.next_time

        # Create particles dataset on the time group

        num_particles = 0
        for species in self.sim.world.species:
            pid_list = self.sim.particle_pool[species.id]
            num_particles += len(pid_list)

        x = numpy.zeros((num_particles, ),
                        dtype = numpy.dtype(PARTICLES_SCHEMA))

        count = 0
        for sid, pid_set in self.sim.particle_pool.iteritems():
            for pid in pid_set:
                pid, particle = self.sim.world.get_particle(pid)
                species = self.sim.world.get_species(sid)
                x['id'][count] = pid.serial
                x['species_id'][count] = sid.serial
                x['position'][count] = particle.position
                count += 1

        dummy = time_group.create_dataset('particles', data = x)

    def write_domains_byHDF5(self):
        if self.HDF5_file is None:
            return

        data_group = self.HDF5_file['data']        
        data_group.attrs['world_size'] = self.sim.world.world_size

        # Require time group
        time_group = data_group.require_group(unicode(self.next_time))
        time_group.attrs['t'] = self.next_time

        # Create shell dataset on the time group

        num_shells = 0
        for domain in self.sim.domains.itervalues():
            num_shells += len(domain.shell_list)

        x = numpy.zeros((num_shells, ), dtype = numpy.dtype(SHELLS_SCHEMA))

        count = 0
        for did, domain in self.sim.domains.iteritems():
            shell_list = domain.shell_list
            for shell_id, shell in shell_list:
                x['id'][count] = shell_id.serial
                x['radius'][count] = shell.shape.radius
                x['position'][count] = shell.shape.position
                count += 1

        dummy = time_group.create_dataset('shells', data = x)

        # Create shell particle association dataset on the time group

        num_assocs = 0
        for domain in self.sim.domains.itervalues():
            if isinstance(domain, Single):
                num_assocs += len(domain.shell_list)
            elif isinstance(domain, Pair):
                num_assocs += 2 * len(domain.shell_list)
            elif isinstance(domain, Multi):
                assert getattr(domain, 'pid_shell_id_map', None), 'Cannot access pid_shell_id_map'
                num_assocs += len(domain.pid_shell_id_map)

        shell_particle_association_schema = \
            [
                ('shell_id', 'u8'),
                ('particle_id', 'u8'),
            ]

        dtype_obj = numpy.dtype(shell_particle_association_schema)
        x = numpy.zeros((num_assocs, ), dtype = dtype_obj)

        count = 0
        for did, domain in self.sim.domains.iteritems():

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

        domain_shell_association_schema = \
            [
                ('shell_id', 'u8', ),
                ('domain_id', 'u8', ),
            ]

        dtype_obj = numpy.dtype(domain_shell_association_schema)
        x = numpy.zeros((num_shells, ), dtype = dtype_obj)

        count = 0
        for did, domain in self.sim.domains.iteritems():
            shell_list = domain.shell_list
            for shell_id, shell in shell_list:
                x['shell_id'][count] = shell_id.serial
                x['domain_id'][count] = did.serial
                count += 1

        dummy = time_group.create_dataset('domain_shell_association', data = x)

        # Create domain dataset on the time group

        domains_schema = \
            [
                ('id', 'u8', ),
                ('kind', 'u4', ),
            ]

        num_domains = len(self.sim.domains.keys())

        dtype_obj = numpy.dtype(domains_schema)
        x = numpy.zeros((num_domains, ), dtype = dtype_obj)

        count = 0
        for did, domain in self.sim.domains.iteritems():
            x['id'][count] = did.serial
            if isinstance(domain, Single):
                x['kind'][count] = 1
            elif isinstance(domain, Pair):
                x['kind'][count] = 2
            else: # must be Multi
                x['kind'][count] = 3
            count += 1

        dummy = time_group.create_dataset('domains', data = x)

    def log(self):
        self.log_time_course()
        self.log_particles()

    def log_time_course(self):
        if self.sim.last_reaction:
            self.write_timecourse()

    def log_particles(self):
        sim = self.sim
        if self.next_time <= sim.t + sim.dt:
            self.write_domains_byHDF5()
            sim.stop(self.next_time)
            self.write_particles()
            self.write_particles_byHDF5()
            self.next_time += self.particle_out_interval

    def start(self):
        self.create_data_group()
        self.write_species_byHDF5()
        self.write_particles()
        self.write_particles_byHDF5()
