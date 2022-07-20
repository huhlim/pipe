#!/usr/bin/env python

import os
import sys
import path
import argparse
import subprocess as sp
from importlib import import_module

from libcommon import *
from libmain import *


def init_refine_meta(arg):
    work_home = path.Dir("%s/%s" % (arg.work_dir, arg.title))
    json_job = work_home.fn("job.json")
    if json_job.status():
        job = Job.from_json(json_job)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
        sp_refine_job = Job.from_json(path.Path(job.sp_refine_job_fn))
        refine_job_s = [Job.from_json(path.Path(fn)) for fn in job.refine_job_fn_s]
        return job, sp_refine_job, refine_job_s
    #
    assert arg.input_pdb is not None
    #
    job = Job(arg.work_dir, arg.title, build=True)
    job.run_type = "refine_meta"
    #
    job.init_home = job.work_home.subdir("init", build=True)
    job.verbose = arg.verbose
    job.keep_tmp = arg.keep
    #
    out = job.init_home.fn("init.pdb")
    if not out.status():
        cmd = ["convpdb.pl", "-out", "generic", arg.input_pdb.short()]
        output = system(cmd, stdout=True, verbose=job.verbose)
        with out.open("wt") as fout:
            fout.write(output)
    job.init_pdb = [out]
    #
    job.sp_refine_job_fn = path.Path(arg.sp_refine_job_fn)
    job.refine_job_fn_s = [path.Path(fn) for fn in arg.refine_job_fn_s]
    sp_refine_job = Job.from_json(path.Path(job.sp_refine_job_fn))
    refine_job_s = [Job.from_json(path.Path(fn)) for fn in job.refine_job_fn_s]
    #
    match_topology(job, sp_refine_job, refine_job_s)
    #
    job.to_json()
    job.append_to_joblist()
    return job, sp_refine_job, refine_job_s


def match_topology(job, sp_refine_job, refine_job_s):
    mdtraj = import_module("mdtraj")
    np = import_module("numpy")
    #
    top_s = []
    top_s.append(mdtraj.load(sp_refine_job.top_fn.short()).top)
    top_s.extend([mdtraj.load(job.top_fn.short()).top for job in refine_job_s])

    def get_atoms(top):
        atom_s = []
        for atom in top.atoms:
            atom_s.append((atom.residue.chain.index, atom.residue.resSeq, atom.name))
        return atom_s

    atom_s = [get_atoms(top) for top in top_s]

    common_atom_s = set(atom_s[0])
    common_atom_s = common_atom_s.intersection(*[set(atom) for atom in atom_s[1:]])
    #
    def get_index_sp(common_s, atom_s):
        select = []
        for i, atom in enumerate(atom_s):
            if atom in common_s:
                select.append(i)
        select = np.array(select, dtype=int)
        return select

    def get_index(common_s, atom_s):
        select = []
        k = -1
        index = np.zeros(len(common_s), dtype=int)
        for i, atom in enumerate(atom_s):
            if atom not in common_s:
                continue
            select.append(i)
            k += 1
            index[common_s.index(atom)] = k
        select = np.array(select, dtype=int)
        return select, index

    #
    #
    job.top_index_s = []
    top_index_fn = job.init_home.fn("top_index.0.npz")
    job.top_index_s.append(top_index_fn)
    sp_refine_top_select = get_index_sp(common_atom_s, atom_s[0])
    sp_refine_top_index = np.arange(sp_refine_top_select.shape[0]).astype(int)
    np.savez(top_index_fn.short(), select=sp_refine_top_select, index=sp_refine_top_index)
    #
    common_top = [atom_s[0][i] for i in sp_refine_top_select]  # -> will be final common top
    #
    for i, atom in enumerate(atom_s[1:]):
        select, index = get_index(common_top, atom)
        top_index_fn = job.init_home.fn("top_index.%d.npz" % (i + 1))
        job.top_index_s.append(top_index_fn)
        np.savez(top_index_fn.short(), select=select, index=index)


def main():
    arg = argparse.ArgumentParser(prog="PREFMD")
    arg.add_argument(dest="title", help="Job title")
    arg.add_argument("-i", "--input", dest="input_pdb", help="input PDB file")
    arg.add_argument("--sp", dest="sp_refine_job_fn")
    arg.add_argument("--meta", dest="refine_job_fn_s", nargs="*")
    arg.add_argument(
        "-d", "--dir", dest="work_dir", default="./", help="working directory (default=./)"
    )
    arg.add_argument(
        "--keep",
        dest="keep",
        action="store_true",
        default=False,
        help="set temporary file mode (default=False)",
    )
    arg.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="set verbose mode (default=False)",
    )
    arg.add_argument(
        "-w",
        "--wait",
        dest="wait_after_run",
        action="store_true",
        default=False,
        help="set running type (default=False)",
    )

    if len(sys.argv) == 1:
        return arg.print_help()
    #
    arg = arg.parse_args()
    if arg.input_pdb is not None:
        arg.input_pdb = path.Path(arg.input_pdb)
    #
    # init
    if arg.title.endswith("job.json"):
        input_json = path.Path(arg.title)
        job = Job.from_json(input_json)
        job.verbose = arg.verbose
        job.keep_tmp = arg.keep
        job.append_to_joblist()
        #
        sp_refine_job = Job.from_json(path.Path(job.sp_refine_job_fn))
        refine_job_s = [Job.from_json(path.Path(fn)) for fn in job.refine_job_fn_s]
    else:
        job, sp_refine_job, refine_job_s = init_refine_meta(arg)

    n_meta = 1 + len(refine_job_s)

    # define topology
    locPREFMD_out = get_outputs(sp_refine_job, "locPREFMD")[:1]
    import_module("define_topology").prep(job, locPREFMD_out[0][0])

    # prod
    import_module("prod_meta").prep(job, 0, sp_refine_job)
    for i, refine_job in enumerate(refine_job_s):
        import_module("prod_meta").prep(job, i + 1, refine_job)
    if not run(job, arg.wait_after_run):
        return
    prod_out = get_outputs(job, "prod_meta")
    #
    # score
    import_module("score").prep_meta(job, job.get_task("prod_meta"))
    if not run(job, arg.wait_after_run):
        return

    # average
    import_module("average_meta").prep(
        job,
        "%s.meta" % job.title,
        [i for i in range(n_meta)],
        path.Path("%s/average.json" % DEFAULT_HOME),
        rule="casp12",
    )
    import_module("average_meta").prep(
        job,
        "%s.cluster" % job.title,
        [i for i in range(n_meta)],
        path.Path("%s/average.json" % DEFAULT_HOME),
        rule="cluster",
    )
    if not run(job, arg.wait_after_run):
        return
    average_out = []
    for output in get_outputs(job, "average_meta", expand="pdb_s"):
        average_out.extend(output[0])
    average_out = average_out[:N_MODEL]

    # model
    job.work_home.chdir()
    model_home = job.work_home.subdir("model", build=True)
    prep_s = []
    for i, out in enumerate(average_out):
        prep_fn = model_home.fn("prep_%d.pdb" % (i + 1))
        if not prep_fn.status():
            system(["cp", out.short(), prep_fn.short()], verbose=job.verbose)
        prep_s.append(prep_fn)
    #
    import_module("scwrl").prep(job, prep_s)
    if not run(job, arg.wait_after_run):
        return
    scwrl_out = get_outputs(job, "scwrl")
    #
    import_module("locPREFMD").prep(job, [out[0] for out in scwrl_out])
    if not run(job, arg.wait_after_run):
        return
    locPREFMD_out = get_outputs(job, "locPREFMD")[-N_MODEL:]
    #
    model_s = []
    for i, out in enumerate(locPREFMD_out):
        model_fn = model_home.fn("model_%d.pdb" % (i + 1))
        if not model_fn.status():
            system(["cp", out[0].short(), model_fn.short()], verbose=job.verbose)
        model_s.append(model_fn)

    # qa
    import_module("qa").prep(job, model_s, path.Path("%s/qa.json" % DEFAULT_HOME))
    if not run(job, arg.wait_after_run):
        return
    qa_out = get_outputs(job, "qa")

    # final
    final_home = job.work_home.subdir("final", build=True)
    for i, out in enumerate(qa_out):
        pdb_fn = final_home.fn("model_%d.pdb" % (i + 1))
        if not pdb_fn.status():
            system(["cp", out[0].short(), pdb_fn.short()], verbose=job.verbose)
    #
    job.remove_from_joblist()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit()
