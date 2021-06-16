import argparse
import logging
import os
import shutil
import signal
import sys
import time
import traceback
import urllib.request
import wget

from bs4 import BeautifulSoup
from functools import partial
from git import Repo
from multiprocessing import Pool, Value
from zipfile import ZipFile, BadZipfile

from config import *
from colorer import *
import model
from utils import *

max_models_download_count = None


def counter_init(args):
    """store the counter for later use"""
    global max_models_download_count
    max_models_download_count = args


def execute_with_workers(function, data, max_downloads, nstreams):
    """Execute tasks in parallel using nstreams processes.
    Handle Ctl-C keyboard interrupt if user want to terminate."""
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)

    max_models_download_count = Value("i", 0)
    pool = Pool(
        nstreams, initializer=counter_init, initargs=(max_models_download_count,)
    )

    try:
        result = pool.map_async(function, data, chunksize=1)
        total_task = len(data)
        while True:
            remaining_task = result._number_left
            logger.info(
                "{0}/{1} processed, downloaded {2} !".format(
                    remaining_task, total_task, max_models_download_count.value
                )
            )
            if remaining_task == 0:
                logger.info("Done!")
                break
            if max_models_download_count.value > max_downloads:
                logger.exception(
                    "Maximum models download count {0} reached!".format(max_downloads)
                )
                pool.terminate()
                sys.exit(-1)
            time.sleep(2)
    except KeyboardInterrupt:
        logger.exception("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
    else:
        logger.debug("ThreadPool normal termination")
        pool.close()

    pool.join()
    return result.get()


def download_zip(download, url, scratch_dir, filename):
    """download zip tarball using wget. To avoid progress
    bars from different threads getting mixed, simply set
    callback outputing none."""
    if download:

        def progress_bar(current, total, width=80):
            return ""

        logger.debug("%s downloading %s", process_name(), url)
        wget.download(url, filename, progress_bar)
        global max_models_download_count
        with max_models_download_count.get_lock():
            max_models_download_count.value += 1

    if os.stat(filename).st_size == 0:
        os.remove(filename)
        message = "Empty archive, download failed"
        status = 1
    else:
        mkdir_p(scratch_dir)
        try:
            with ZipFile(filename, "r") as zf:
                zf.extractall(scratch_dir)
        except BadZipfile as e:
            logger.warn("Zip Error : {}".format(e))
            message = str(e)
            status = 1
        else:
            message = "Archive successfully extracted"
            status = 0

    return status, message


def download_git(download, url, branch_name, data_dir, scratch_dir):
    if download:
        logger.debug("%s cloning %s", process_name(), branch_name)
        Repo.clone_from(url, data_dir, branch=branch_name, depth=1)
        global max_models_download_count
        with max_models_download_count.get_lock():
            max_models_download_count.value += 1

    def ignore_files(dir, files):
        return [f for f in files if os.path.islink(os.path.join(dir, f))]

    mkdir_p(data_dir)
    remove_dir(scratch_dir, NMODLDB_CACHE_DIR)
    shutil.copytree(data_dir, scratch_dir, ignore=ignore_files)

    message = "Repository uccessfully cloned"
    status = 0
    return status, message


def download_model(use_cache, keep_download, mod):
    """Download model from ModelDB or BBP repository.
    Currently we have considered only ModelDB and BBP
    models which are stored as zip or git repository.
    """
    tid = process_name()
    mod_dir = "%s/%s" % (mod.output_dir, "mod")
    data_dir = "%s/%s" % (mod.cache_dir, "data")
    scratch_dir = "%s/%s" % (mod.cache_dir, "scratch")

    result = {"status": 0, "message": "", "url": mod.url}

    try:
        if os.path.exists(data_dir) and use_cache:
            download = False
        else:
            download = True
            remove_dir(data_dir, NMODLDB_CACHE_DIR)

        mkdir_p(data_dir)
        mkdir_p(mod_dir)

        if mod.is_zip:
            filename = mod.output_filename
            status, message = download_zip(download, mod.url, scratch_dir, filename)
        elif mod.is_git:
            model_dir = "%s/%s" % (data_dir, mod.output_filename)
            status, message = download_git(
                download, mod.url, mod.git_branch, model_dir, scratch_dir
            )

        copy_mod_files(mod, scratch_dir, mod_dir)
        remove_dir(scratch_dir, NMODLDB_CACHE_DIR)
        logger.debug("%s MOD files copied to %s", tid, mod_dir)

    except BaseException:
        error = traceback.format_exception(*sys.exc_info())
        message = "Error %s while cloning model %s" % (error, mod.name)
        status = 1

    if not keep_download:
        remove_dir(data_dir, NMODLDB_CACHE_DIR)

    result["status"] = status
    result["message"] = message
    mod.result = result
    return mod


def download_bluebrain_models(max_downloads, nstreams, reuserepo, keep_download):
    """Read model information from YAML file and download the models
    from git repository"""
    logger.info("Downloading BlueBrain models")
    models = []
    modelsinfo = get_bbp_models_info()

    for tag, modelinfo in modelsinfo.iteritems():
        url = modelinfo["url"]
        branch = modelinfo["branch"]
        mod = model.ModelInfo(url)
        mod.name = tag
        mod.is_git = True
        mod.output_dir = "%s/bluebrain/%s" % (NMODLDB_ROOT_DIR, tag)
        mod.cache_dir = "%s/bluebrain/%s" % (NMODLDB_CACHE_DIR, tag)
        mod.output_filename = "repository"
        mod.git_branch = branch
        mod.neuron_model = True
        mod.group = "BlueBrain"
        mod.simulator = "NEURON"
        models.append(mod)

    func = partial(download_model, reuserepo, keep_download)
    models = execute_with_workers(func, models, max_downloads, nstreams)

    return models


def download_modeldb_models(
    nmodels, max_downloads, nstreams, syncmodeldb, reuserepo, keep_download
):
    """Download zip files from ModelDB"""
    if syncmodeldb or (not os.path.isfile(MDB_YAML_DATA_FILE)):
        logger.info("Syncing models information from ModelDB")
        sync_modeldb_models()
    else:
        logger.info("Using local ModelDB models cache")

    models = yaml.load(open(MDB_YAML_DATA_FILE), Loader=yaml.Loader)
    func = partial(download_model, reuserepo, keep_download)

    logger.info("Downloading ModelDB models")
    models = execute_with_workers(func, models, max_downloads, nstreams)

    return models


def extract_modeldb_models(response, nmodels=None):
    """Extract ModelDB models information from HTML file

    ModelDB doesn't have REST api to query/download models.
    From HTML response, we iterate and extract model id which
    can be used to downloa zip file of model"""
    soup = BeautifulSoup(response, "html.parser")
    records = soup.select("tr.searchTable")

    def get_element(x, i):
        if len(x) > i:
            return x[i].get_text()
        return None

    def get_citation(msg):
        citation = msg[msg.find("(") + 1 : msg.find(")")]
        return citation

    def is_neuron_model(sim):
        return True if "neuron" in sim.lower() else False

    def download_url(model_id):
        return "%s?o=%s&a=23&mime=application/zip" % (MDB_DOWNLOAD_URL, model_id)

    models = []
    nmodels = len(records) if nmodels is None else min(nmodels, len(records))

    for record in records[:nmodels]:
        columns = record.find_all("td")
        modeldb_id = get_element(columns, 1)
        description = get_element(columns, 2)
        simulator = get_element(columns, 3)
        citation = get_citation(description)
        url = download_url(modeldb_id)

        mod = model.ModelInfo(url)
        mod.name = modeldb_id
        mod.is_zip = True
        mod.output_dir = "%s/modeldb/%s" % (NMODLDB_ROOT_DIR, modeldb_id)
        mod.cache_dir = "%s/modeldb/%s" % (NMODLDB_CACHE_DIR, modeldb_id)
        mod.output_filename = "%s/data/%s.zip" % (mod.cache_dir, modeldb_id)
        mod.description = description
        mod.citation = citation
        mod.neuron_model = is_neuron_model(simulator)
        mod.group = "ModelDB"
        mod.simulator = simulator

        models.append(mod)

    return models


def sync_modeldb_models():
    """Query ModelDB for all models andstore response as HTML file"""
    url = MDB_SEARCH_URL
    url += "?ddl_class=19"
    url += "&fd_oid=on"
    url += "&fd_name=on"
    url += "&fc_a25="
    url += "&fd_a114=on"
    url += "&fd_a299=on"
    url += "&fulltextSearch="
    url += "&btn_Search=Search"

    response = urllib.request.urlopen(url).read()
    models = extract_modeldb_models(response)

    with open(MDB_YAML_DATA_FILE, "w") as outfile:
        yaml.dump(models, outfile, default_flow_style=False)


def generate_database(options):
    if options.clean_db:
        logger.warn("Removing NMODL database from %s", NMODLDB_ROOT_DIR)
        remove_dir(NMODLDB_ROOT_DIR)

    if options.clean_cache:
        logger.warn("Removing NMODL cache from %s", NMODLDB_CACHE_DIR)
        remove_dir(NMODLDB_CACHE_DIR)

    mkdir_p(MDB_ROOT_DIR)
    mkdir_p(NMODLDB_CACHE_DIR)

    models = []
    if options.download_modeldb:
        result = download_modeldb_models(
            options.max_update_count,
            options.max_download_count,
            options.download_streams,
            options.modeldb_sync,
            options.modeldb_cache,
            options.keep_download,
        )
        result = sorted(result, key=lambda k: int(k.name))
        models.extend(result)

    if options.download_bbp:
        result = download_bluebrain_models(
            options.max_download_count,
            options.download_streams,
            options.bbp_cache,
            options.keep_download,
        )
        result = sorted(result, key=lambda k: k.name)
        models.extend(result)

    if not options.skip_stats and (options.download_modeldb or options.download_bbp):
        failed = sum([1 for mod in models if mod.result["status"] == 1])
        success = len(models) - failed

        result = {}
        result["downloaded"] = success
        result["failed"] = failed
        result["models"] = models

        write_yaml_report(result, NMODLDB_STATS_FILE)
        logger.warn(
            "DBGEN STATS : SUCCESS %d, FAILED %d, LOG FILE %s",
            success,
            failed,
            NMODLDB_STATS_FILE,
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="NMODL database generator from ModelDB"
    )
    parser.add_argument(
        "--download-modeldb",
        dest="download_modeldb",
        action="store_true",
        default=False,
        help="Download models from ModelDB (False).",
    )
    parser.add_argument(
        "--download-bbp",
        dest="download_bbp",
        action="store_true",
        default=False,
        help="Download models from BBP Neurodamus (False).",
    )
    parser.add_argument(
        "--max-models-update-count",
        dest="max_update_count",
        type=int,
        default=0,
        help="Maximum number of models to update (0).",
    )
    parser.add_argument(
        "--use-modeldb-cache",
        dest="modeldb_cache",
        action="store_true",
        default=False,
        help="Use existing zip files from ModelDB if available (False).",
    )
    parser.add_argument(
        "--use-bbp-cache",
        dest="bbp_cache",
        action="store_true",
        default=False,
        help="Use exisiting BBP models if available (False).",
    )
    parser.add_argument(
        "--sync-modeldb",
        dest="modeldb_sync",
        action="store_true",
        default=False,
        help="Synchronize models information from ModelDB (False).",
    )
    parser.add_argument(
        "--keep-downloads",
        dest="keep_download",
        action="store_true",
        default=False,
        help="Keep zip/git downloads from ModelDB/BBP (False).",
    )
    parser.add_argument(
        "--clean-db",
        dest="clean_db",
        action="store_true",
        default=False,
        help="Remove exisiting NMODL database (False).",
    )
    parser.add_argument(
        "--clean-cache",
        dest="clean_cache",
        action="store_true",
        default=False,
        help="Remove exisiting cached files (False).",
    )
    parser.add_argument(
        "--download-streams",
        dest="download_streams",
        type=int,
        default=8,
        help="Number of parallel download streams to use (8).",
    )
    parser.add_argument(
        "--max-models-download-count",
        dest="max_download_count",
        type=int,
        default=25,
        help="Maximum number of models to download from ModelDB or BBP (25).",
    )
    parser.add_argument(
        "--skip-stats",
        dest="skip_stats",
        action="store_true",
        default=False,
        help="Do not write models/stats.yaml file (False).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    generate_database(parse_args())
