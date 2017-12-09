
# coding: utf-8

# In[ ]:

from new_cfg import new_rc
new_rc()


# In[5]:

clear_tables()


# In[ ]:

all_plots()


# In[2]:

from skysearch_lib import np
from skysearch_lib import os
from skysearch_lib import Table
from skysearch_lib import vstack
from skysearch_lib import MPI_TABLE_DIR
from skysearch_lib import TABLE_DIR
from skysearch_lib import PLOT_DIR
from skysearch_lib import T_FRMT
from skysearch_lib import H5_PTH
from skysearch_lib import record_table
from skysearch_lib import get_data_plot
from skysearch_lib import plot_full_halo
from skysearch_lib import save_plot
from skysearch_lib import update_plot
from skysearch_lib import grid_list
from skysearch_lib import load_grid
from skysearch_lib import mu_idx
PLOT_DIR = os.path.join(PLOT_DIR, 'summary_plots')


# In[3]:

from matplotlib import pyplot as plt
from time import sleep


# In[4]:

def collect_tables(_return=True, _save=False, _sort=True):
    master_table = None
    for fh_ in os.listdir(MPI_TABLE_DIR):
        table_fh_ = os.path.join(MPI_TABLE_DIR, fh_)
        new_table = Table.read(
            table_fh_, 
            format=T_FRMT, 
            path=H5_PTH)
        if not master_table:
            master_table = new_table
        else:
            master_table = vstack(
                [master_table, new_table], 
                join_type='exact', 
                metadata_conflicts='silent')
    if _sort:
        master_table.sort(['halo', 'radius'])
    if _save:
        master_table.pprint(max_lines=np.random.randint(500,high=50000), max_width=-1)
        master_table.write(
            os.path.join(TABLE_DIR, 'groupfinder') + 'main_table', 
            format=T_FRMT, 
            path=H5_PTH,
            overwrite=True)
    if _return:
        return master_table

def clear_tables():
    for fh in os.listdir(MPI_TABLE_DIR):
        os.remove(os.path.join(MPI_TABLE_DIR, fh))

def all_plots():
    plot_name = 'all halos - summary plot'
    fig = plt.figure(figsize=(15, 30))
    fig.suptitle(plot_name)
    for plot_panel, column_name in enumerate(table.keys()):
        ax = fig.add_subplot(8, 3, plot_panel + 1)
        ax.set_title(column_name)
        ax.scatter(table['radius'], table[column_name], s=.51, c=table['halo'])
    fig.savefig(os.path.join(PLOT_DIR, plot_name))


# In[ ]:

try:
    fntsz = 27
    blank = '                        '
    while True:
        for fh in os.listdir(PLOT_DIR):
            os.remove(os.path.join(PLOT_DIR, fh))
            os.sys.stdout.write('\rremoving: ' + fh + blank)
            os.sys.stdout.flush()
        os.sys.stdout.write('\rloading table' + blank)
        os.sys.stdout.flush()
        table = collect_tables()
        os.sys.stdout.write('\rplotting' + blank)
        os.sys.stdout.flush()
        for plot_panel, column_name in enumerate(table.keys()):
            fig = plt.figure(figsize=(20, 10))
            ax = fig.add_subplot(111)
            ax.set_title(column_name, fontsize=fntsz)
            os.sys.stdout.write('\r' + column_name + blank)
            os.sys.stdout.flush()
            ax.scatter(table['radius'], table[column_name], s=10, c=table['halo'])#c='k')#
            ax.set_xlabel('Radius Kpc', fontsize=int(fntsz * 0.75))
            ax.set_ylabel(column_name, fontsize=fntsz)
            ax.axes.grid()
            ax.tick_params(labelsize=22)
            try:
                fig.savefig(os.path.join(PLOT_DIR, str(plot_panel) + '_' + column_name))
            except SystemError as e:
                print('didnt save last plot')
            plt.close()
        for i in range(60):
            os.sys.stdout.write('\rtime till plot: ' + str(60-i) + ' seconds' + blank)
            sleep(1)
            os.sys.stdout.flush()
except KeyboardInterrupt as e:
    os.sys.stdout.write('\r' + blank + blank)
    os.sys.stdout.flush()
    print(e)


# In[ ]:

plt.show()


# In[ ]:

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='polar')
ax.bar(table[''], table[''])


# In[8]:

table = collect_tables()
table.keys()


# In[9]:

for num in np.unique(table['halo']):
    h_num = str(num)
    if len(h_num) < 2:
        halo = 'halo0' + h_num
    else:
        halo = 'halo' + h_num
    #os.sys.stdout.write('\r' + halo + ': making full plot            ')
    #os.sys.stdout.flush()
    #plot_full_halo(halo, _d_cut=1, _out_dir=PLOT_DIR)
    
    os.sys.stdout.write('\r' + halo + ': making grouping plot            ')
    os.sys.stdout.flush()
    plot = get_data_plot(halo)
    tot = 1e0 * len(table[table['halo']==num])
    for i, row in enumerate(table[table['halo']==num]):
        update_plot(plot, [row['deg0'], row['r1'] - row['r0'], row['extent'], row['r0']], halo, row['domsat_id'])
        os.sys.stdout.write('\r' + halo + ': updating ' + str(round((1e2 * i) / tot, 2)) + ' %          ')
        os.sys.stdout.flush()
    save_plot(plot, halo, _out_dir=PLOT_DIR)
    os.sys.stdout.write('\rdone                               ')
    os.sys.stdout.flush()
    
    


# In[ ]:

grids = grid_list()


# In[ ]:

plt.cm.spectral_r


# In[ ]:

# ('name', log=np.log10(), cmap=plt.cm.plasma, plot=True)
plot_config = [
    ('n_stars', True, plt.cm.bone_r, True, (0, 4)),
    ('AB_mag', False, plt.cm.Paired, True, (-1, 10)),
    ('sat_age', False, plt.cm.spectral_r, True, (3, 12)),
    ('sat_number', False, plt.cm.Paired, True, (0, 1115)),
    ('r_proj', False, None, False, (0, 0)),
    ('phi', False, None, False, (0, 0))
]


# In[ ]:

last_sat = 0
for _grid in grids:
    halo = _grid[0]
    grid = load_grid(_grid[1])
    
    
    for i in range(grid.shape[2] - 2):
        name, log, _cmap, plot, minmax = plot_config[i]
        _vmin, _vmax = minmax
        
        if not plot:
            continue
            
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        
        _std = grid[:, :, i].std()
        if name in ['AB_mag']:
            _vmin = np.power(10, grid[:, :, i]).mean() - (3.0 * _std)
            _vmax = np.power(10, grid[:, :, i]).mean() + (3.0 * _std)
            
        if name in ['sat_number']:
            _vmin = np.unique(grid[:, :, i])[1]
            _vmax = np.unique(grid[:, :, i])[-1]
        
        ax.set_title(halo + ' ' + name + '\n' + str(round(_vmin, 2)) + ' / ' + str(round(_vmax, 2)))
        
        if log:
            ax.pcolormesh(np.log10(grid[:, :, i]), cmap=_cmap, vmin=_vmin, vmax=_vmax)
            
        elif name == 'AB_mag':
            ax.pcolormesh(np.power(10, grid[:, :, i]), cmap=_cmap, vmin=_vmin, vmax=_vmax)
            
        else:
            ax.pcolormesh(grid[:, :, i], cmap=_cmap, vmin=_vmin, vmax=_vmax)
        
        #cm = ax.contour(grid[:, :, -2], labels=['25', '75', '125'], levels=[25, 75 ,125], c='k')
        #cm.clabel([25, 75, 125], format='f1%', inline=True, color='k')
        ax.set_xlim([200, 400])
        ax.set_ylim([200, 400])
        
        fig.savefig(os.path.join(PLOT_DIR, halo + '_' + name))
        
        plt.close()
    


# In[ ]:

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='polar')
grd = load_grid(grids[3][1])
idx = np.nonzero(mu_idx(grd, 50, 150)[1])[0]
ax.scatter(grd[:, :, 5], grd[:,:,4], s=.015, c=np.log10(grd[:,:,0]), cmap=plt.cm.plasma, vmin=0, vmax=3.2)
ax.
plt.show()


# In[ ]:



