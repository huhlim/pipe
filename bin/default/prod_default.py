import os
import json

WORK_HOME = os.getenv("work")

options = {}

options["ff"] = {}
options["ff"]["use_modified_CMAP"] = True
options["ff"]["toppar"] = [
    "%s/ff/c36m.CMAPmod.hv/c36m_ecmap.prm" % WORK_HOME,
    "%s/ff/c36m.CMAPmod.hv/c36m_ecmap.rtf" % WORK_HOME,
    "%s/ff/c36m.CMAPmod.hv/toppar_water_ions.str" % WORK_HOME,
    "%s/ff/c36m.CMAPmod.hv/toppar_ions_won.str" % WORK_HOME,
]
options["ff"]["cgenff"] = [
    "%s/ff/cgenff.hv/top_all36_cgenff.rtf" % WORK_HOME,
    "%s/ff/cgenff.hv/par_all36_cgenff.prm" % WORK_HOME,
]

options["openmm"] = {}
options["openmm"]["platform"] = "CUDA"

options["md"] = {}

options["restart"] = True
options["md"]["dyntstep"] = 0.004
options["md"]["dyntemp"] = 360.0
options["md"]["lang"] = True
options["md"]["langfbeta"] = 0.01

SIMUL_FRAME_INTERVAL = 50.0  # 50 ps
SIMUL_TRAJ_LENGTH = 10000.0  # 10 ns
options["md"]["dynoutfrq"] = int(SIMUL_FRAME_INTERVAL / options["md"]["dyntstep"])
options["md"]["dynsteps"] = int(SIMUL_TRAJ_LENGTH / options["md"]["dyntstep"])
options["md"]["time_limit"] = 3600.0 * 36
options["md"]["iter"] = 10

options["restraint"] = {}
options["restraint"]["mode"] = "dual.anneal"
options["restraint"]["distance"] = [0.05, 2.0]
options["restraint"]["Cartesian"] = [0.025, 4.0]

with open("prod.json", "wt") as fout:
    fout.write(json.dumps(options, indent=2))
