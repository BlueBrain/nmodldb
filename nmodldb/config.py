"""Global configuration variables used in the projectz"""

import os

MDB_SEARCH_URL = "https://senselab.med.yale.edu/ModelDB/searchFullText.cshtml"
MDB_DOWNLOAD_URL = "https://senselab.med.yale.edu/ModelDB/eavBinDown.cshtml"

ROOT_DIR = os.path.abspath(__file__ + "/../../")
ROOT_DIR = os.path.relpath(ROOT_DIR, os.getcwd())

NMODLDB_ROOT_DIR = "%s/models/db" % ROOT_DIR
NMODLDB_CACHE_DIR = "%s/cache" % ROOT_DIR

BBP_ROOT_DIR = "%s/bluebrain" % NMODLDB_ROOT_DIR
MDB_ROOT_DIR = "%s/modeldb" % NMODLDB_ROOT_DIR

MDB_YAML_DATA_FILE = "%s/models.yaml" % MDB_ROOT_DIR
NMODLDB_STATS_FILE = "%s/models/stats.yaml" % ROOT_DIR
BBP_MODEL_INFO_FILE = "%s/models/bbp.yaml" % ROOT_DIR
