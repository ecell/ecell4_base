from ecell4.core import *
from ecell4.ode import *

import h5py

def hdf5load_ode(h5filename):
    with h5py.File(h5filename, "r") as f:
        t_str = f.keys()[0]
        h5grp_path_common = "/%s/CompartmentSpace/" % t_str
        h5grp_path_num_molecules = h5grp_path_common + "num_molecules"
        h5grp_path_species = h5grp_path_common + "species"
        h5grp_path_space_attr = h5grp_path_common + "space_attributes"

        sid_species_dict = dict()
        species_num_molecules_dict = dict()
        for (d_sid, d_serial) in f[h5grp_path_species]:
            sid_species_dict[d_sid] = d_serial
        for (d_sid, d_num_molecules) in f[h5grp_path_num_molecules]:
            species_num_molecules_dict[ sid_species_dict[d_sid] ] = d_num_molecules
        t = float(t_str)
        volume = f[h5grp_path_space_attr][0]["volume"]    # XXX
        rnd = 0         # XXX
        w = ODEWorld(volume)
        w.set_t(t)
        for (sp, n) in species_num_molecules_dict.items():
            #print "{} -- {}".format(sp, n)
            w.add_molecules(Species(sp), n)
        #print volume
        return w
