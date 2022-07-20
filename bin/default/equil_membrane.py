import os
import glob
import json

WORK_HOME = os.getenv("work")

options = {}

options["ff"] = {}
options["ff"]["use_modified_CMAP"] = False
options["ff"]["toppar"] = [
    "%s/ff/c36m.charmm_gui/top_all36_prot.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_all36m_prot.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/top_all36_na.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_all36_na.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/top_all36_carb.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_all36_carb.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/top_all36_lipid.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_all36_lipid.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/top_all36_cgenff.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_all36_cgenff.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/top_interface.rtf" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/par_interface.prm" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_nano_lig.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_nanolig_patch.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_synthetic_polymer.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_synthetic_polymer_patch.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_polymer_solvent.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_water_ions.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_dum_noble_gases.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_ions_won.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_arg0.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_c36m_d_aminoacids.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_fluoro_alkanes.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_heme.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_na_combined.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_retinol.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_prot_modify_res.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_na_nad_ppi.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_na_rna_modified.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_bacterial.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_cardiolipin.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_cholesterol.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_dag.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_inositol.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_lps.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_miscellaneous.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_model.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_prot.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_sphingo.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_yeast.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_hmmm.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_detergent.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_lipid_ether.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_carb_glycolipid.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_carb_glycopeptide.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_carb_imlab.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_label_spin.str" % WORK_HOME,
    "%s/ff/c36m.charmm_gui/toppar_all36_label_fluorophore.str" % WORK_HOME,
]
options["ff"]["cgenff"] = [
    "%s/ff/cgenff/top_all36_cgenff.rtf" % WORK_HOME,
    "%s/ff/cgenff/par_all36_cgenff.prm" % WORK_HOME,
]

options["openmm"] = {}
options["openmm"]["platform"] = "CUDA"

options["md"] = {}
options["md"]["dyntstep"] = 0.002
options["md"]["dyntemp"] = 360.0
options["md"]["lang"] = True
options["md"]["langfbeta"] = 0.01
options["md"]["ion_conc"] = 0.0
options["md"]["solvate"] = 9.0
options["md"]["heat"] = [10000, 50.0, 25.0]  # n_steps, init.T, incr.T
options["md"]["equil"] = [500000]  # total_steps (1 ns)

with open("equil_membrane.json", "wt") as fout:
    fout.write(json.dumps(options, indent=2))
