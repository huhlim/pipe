import os
import json

WORK_HOME = os.getenv("work")

options = {}

options['ff'] = {}
options['ff']['toppar'] = [ \
             '%s/ff/c36m/par_all36m_prot.prm'%WORK_HOME,\
             '%s/ff/c36m/top_all36_prot.rtf'%WORK_HOME,\
             '%s/ff/c36m/toppar_water_ions.str'%WORK_HOME,\
             '%s/ff/c36m/toppar_ions_won.str'%WORK_HOME,\
             ]
options['ff']['cgenff'] = [ \
             '%s/ff/cgenff/top_all36_cgenff.rtf'%WORK_HOME, \
             '%s/ff/cgenff/par_all36_cgenff.prm'%WORK_HOME, \
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
#options['md']['equil'] = [500000]           # total_steps (1 ns)
options['md']['equil'] = [250000]           # total_steps (1 ns)
#options['md']['prod'] = [2, int(5000./options['md']['dyntstep']), 25000]
options['md']['prod'] = [2, int(1000./options['md']['dyntstep']), 5000]
options['md']['force_const'] = 0.0

with open("qa.json", 'wt') as fout:
    fout.write(json.dumps(options, indent=2))

