from mpi4py import MPI
import skysearch_lib as ss_lib
import sys
import os
import skysurvey
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

record_table = ss_lib.record_table()


# Load halo grid into mem
grid_dir_path = os.path.join(skysurvey.grid_dir, '601')
grids = [fh for fh in os.listdir(grid_dir_path) if fh.endswith('npy')]


for g_fh in grids:

    #print name, rank, size

    grid_fh = os.path.join(grid_dir_path, g_fh)

    grd = np.load(grid_fh)

    grid = ss_lib.fix_rslice(grd)

    #print grid_fh, 'loaded'

    for r_in, r_out, deg_0, deg_1 in ss_lib.regions_mpi(rank=rank, size=size):

        r = str(r_out - r_in)
        #print 'index', r_in, r_out, deg_0, deg_1
        grid_idx = np.nonzero(
            np.logical_and(
                grid[:, :, 4] >= r_in,
                grid[:, :, 4] < r_out))


        mu = grid[:, :, 0][grid_idx].mean()
        xbox_min = mu * 0.0005
        ##print 'mu', mu

        one_before = False
        run_length = 0

        min_xbox = []
        max_xbox = []
        mean_xbox = []

        alims = np.nonzero(
            np.logical_and(
                np.logical_and(
                    grid[:, :, 5] >= deg_0,
                    grid[:, :, 5] < deg_1),
                np.logical_and(
                    grid[:, :, 4] >= r_in,
                    grid[:, :, 4] < r_out)))

        xbox = ((grid[:, :, 0][alims] - mu) / mu)

        if xbox.mean() > xbox_min:

            if not one_before:

                run_length = 0

            run_length += 1
            one_before = True
            min_xbox.append(xbox.min())
            mean_xbox.append(xbox.mean())
            max_xbox.append(xbox.max())

        else:

            if one_before:

                row = [r,
                       r_strt,
                       r_stop,
                       grids[grid_number].split('_')[0],
                       xbox_min,
                       round(np.asarray(min_xbox).min(), 4),
                       round(np.asarray(mean_xbox).mean(), 4),
                       round(np.asarray(max_xbox).max(), 4),
                       run_length,
                       0.0,
                       grid[:, :, 0][alims].sum(),
                       len(alims[0]),
                       0.0,
                       0.0,
                       0.0]

                record_table.add_row(row)

            one_before = False
            run_length = 0


        m0 = '\n' + r + ' Kpc   '
        m1 = str(deg_0) + ' deg  '
        m2 = '  '
        m3 = str(one_before) + ' '
        m4 = str(run_length) + ' '
        #sys.stdout.write(m0 + m1 + m2 + m3 + m4)
        #sys.stdout.flush()

  record_table.pprint(max_width=200, center='^')
