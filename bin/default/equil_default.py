import os
import json

WORK_HOME = os.getenv("work")

options = {}

options['ff'] = {}
options['ff']['use_modified_CMAP'] = True
options['ff']['toppar'] = [ \
             '%s/ff/c36m.CMAPmod.hv/c36m_ecmap.prm'%WORK_HOME, \
             '%s/ff/c36m.CMAPmod.hv/c36m_ecmap.rtf'%WORK_HOME, \
             '%s/ff/c36m.CMAPmod.hv/toppar_water_ions.str'%WORK_HOME, \
             '%s/ff/c36m.CMAPmod.hv/toppar_ions_won.str'%WORK_HOME, \
             ]
options['ff']['cgenff'] = [ \
             '%s/ff/cgenff.hv/top_all36_cgenff.rtf'%WORK_HOME, \
             '%s/ff/cgenff.hv/par_all36_cgenff.prm'%WORK_HOME, \
             ]

options['openmm'] = {}
options['openmm']['platform'] = 'CUDA'

options['md'] = {}
options['md']['dyntstep'] = 0.002
options['md']['dyntemp'] = 298.15
options['md']['lang'] = True
options['md']['langfbeta'] = 0.01
options['md']['ion_conc'] = 0.
options['md']['solvate'] = 9.0
options['md']['heat'] = [20000, 50., 25.]   # n_steps, init.T, incr.T
options['md']['equil'] = [500000]           # total_steps (1 ns)

with open("equil.json", 'wt') as fout:
    fout.write(json.dumps(options, indent=2))

