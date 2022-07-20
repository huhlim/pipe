import os
import json

WORK_HOME = os.getenv("work")

options = {}

options["ff"] = {}
options["ff"]["toppar"] = [
    "%s/ff/c36m/par_all36m_prot.prm" % WORK_HOME,
    "%s/ff/c36m/top_all36_prot.rtf" % WORK_HOME,
    "%s/ff/c36m/toppar_water_ions.str" % WORK_HOME,
    "%s/ff/c36m/toppar_ions_won.str" % WORK_HOME,
]
options["ff"]["cgenff"] = [
    "%s/ff/cgenff/top_all36_cgenff.rtf" % WORK_HOME,
    "%s/ff/cgenff/par_all36_cgenff.prm" % WORK_HOME,
]

options["openmm"] = {}
options["openmm"]["platform"] = "CUDA"

options["md"] = {}

options["md"]["dyntstep"] = 0.002
options["md"]["dyntemp"] = 298.15
options["md"]["lang"] = True
options["md"]["langfbeta"] = 0.01
options["md"]["ion_conc"] = 0.0
options["md"]["solvate"] = 9.0
options["md"]["heat"] = [
    [50, 100, 150, 200, 250, options["md"]["dyntemp"]],
    [1000, 1000, 1000, 1000, 1000, 5000],
]

with open("average.json", "wt") as fout:
    fout.write(json.dumps(options, indent=2))
