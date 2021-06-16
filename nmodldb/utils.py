import errno
import logging
import multiprocessing
import os
import shutil
import yaml

from config import *
from model import *


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def process_name():
    return multiprocessing.current_process().name


def write_yaml_report(models, report):
    with open(report, "w") as outfile:
        output = {}
        output["models"] = models
        yaml.dump(output, outfile, default_flow_style=False)


def is_subdir_in_dir(subdir, rootdir):
    # make both absolute
    rootdir = os.path.join(os.path.realpath(rootdir), "")
    subdir = os.path.realpath(subdir)

    # return true, if the common prefix of both is equal to rootdir
    # e.g. /a/b/c/d.rst and rootdir is /a/b, the common prefix is /a/b
    return os.path.commonprefix([subdir, rootdir]) == rootdir


def copy_mod_files(model, source, dest):
    def is_synapse(filename):
        with open(filename, "r") as f:
            return True if "POINT_PROCESS" in f.read() else False

    def verbatim_count(filename):
        with open(filename, "r") as f:
            return f.read().count("ENDVERBATIM")

    def line_count(filename):
        return sum(1 for line in open(filename))

    def modfile_name(filename):
        fname = os.path.basename(filename)
        return os.path.splitext(fname)[0]

    for root, _, filenames in os.walk(source):
        for filename in filenames:
            if filename.endswith((".mod", ".inc")):
                srcfile = os.path.join(root, filename)
                destfile = os.path.join(dest, filename)
                shutil.copy(srcfile, destfile)
                model.nmod += 1
                nlines = line_count(srcfile)
                nverbatim = verbatim_count(srcfile)

                if is_synapse(srcfile):
                    model.nsynapses += 1
                    name = modfile_name(srcfile)
                    synapse = SynapseInfo(name)
                    synapse.nlines = nlines
                    synapse.nverbatim = nverbatim
                    model.synapses.append(synapse)
                else:
                    model.nchannels += 1
                    name = modfile_name(srcfile)
                    channel = ChannelInfo(name)
                    channel.nlines = nlines
                    channel.nverbatim = nverbatim
                    model.channels.append(channel)

                model.nverbatim += nverbatim

    # for stats report to be consistent across runs,
    # sort channels/synapses list
    model.channels.sort(key=lambda x: x.name)
    model.synapses.sort(key=lambda x: x.name)


def remove_dir(subdirectory, directory=None):
    if os.path.isdir(subdirectory):
        tid = process_name()
        if not directory or is_subdir_in_dir(subdirectory, directory):
            logging.debug("%s : Removing directory %s", tid, subdirectory)
            shutil.rmtree(subdirectory)
        else:
            logging.error("%s : Logical error? Trying to remove %s", tid, subdirectory)


def get_bbp_models_info():
    def flatten_model(models, name, submodels):
        modfiles = models[name]["mod"]
        for modelname in submodels:
            model = models[modelname]
            includes = model["includes"]
            files = flatten_model(models, modelname, includes)
            modfiles.extend(files)
            model["includes"] = []
        return modfiles

    file = open(BBP_MODEL_INFO_FILE, "r")
    models = yaml.load(file)
    result = {}

    for name, model in models.items():
        if "includes" in model:
            submodels = model["includes"]
            flatten_model(models, name, submodels)
            if model["model"]:
                result[name] = model

    return result
