# NMODL Database and Generator

This repository stores all NMODL models (>6300 mod files) parsed by [NMODL Framework](https://github.com/BlueBrain/nmodl).

For more information, see this issue : [BlueBrain/nmodl/issues/192](https://github.com/BlueBrain/nmodl/issues/192)

## Installing dependencies

We use tox for dependency management and script execution. That means, the only thing you need to
install yourself is tox:

```
pip install tox
```

## syncing / downloading DBs


Downloading modeldb is done with:

```
tox -e download-modeldb
```

Downloading BBP's models is done with:

```
tox -e download-bbpdb
```

> There is a catch, though. You need to have access to BBPs internal repositories and you need a
> special YAML config file. If you are not BBP HPC Team, you will probably not have the appropriate
> access.


### First-time download

The default for maximum number of file downloads is set to 25, this is to avoid large numbers of
requests hammering the modeldb server. For the initial download you can use :

```
tox -e init-modeldb
```

to pass to the `dbgen.py` script the extra argument to fetch up to 2000 models.

### Other operations

To clean/remove existing databases:
```
tox -e clean-db
```

To clean/remove existing cache:
```
tox -e clean-cache
```

