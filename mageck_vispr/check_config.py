__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import sys


class ConfigError(Exception):
    def __init__(self, msg, key=None, entry=None):
        entry = " (key={}, entry={})".format(key, entry) if entry is not None and key is not None else ""
        super().__init__("Error in configuration file{}: {}".format(entry, msg))


def is_file(key, entry):
    if not os.path.exists(str(entry)):
        raise ConfigError("File does not exist.", key, entry)


def is_str(key, entry):
    if not isinstance(entry, str):
        raise ConfigError("Expecting a string.", key, entry)


def is_bool(key, entry):
    if not isinstance(entry, bool):
        raise ConfigError("Expecting a string.", key, entry)


def is_int(key, entry):
    if not isinstance(entry, int):
        raise ConfigError("Expecting an integer.", key, entry)

def is_str_or_int(key, entry):
    if not isinstance(entry, int) and not isinstance(entry,str):
        raise ConfigError("Expecting an integer or a string.", key, entry)


def is_sample(key, entry):
    if isinstance(entry, str):
        is_file(key, entry)
    elif isinstance(entry, list):
        for f in entry:
            is_file(key, f)
    else:
        raise ConfigError("Expecting a list of files or a single file.", key, entry)


def is_samples(key, entry):
    if not isinstance(entry, dict):
        raise ConfigError("Expecting an assignment of samples to FASTQ files.", key, entry)
    for sample, entry in entry.items():
        is_sample(sample, entry)


def is_experiment(key, entry, msg="Expecting treatment and control samples or a design matrix."):
    if key.startswith("myexperiment"):
        print("Warning: You use the experiment name '{}' in your config file which is "
              "intended as a placeholder. This won't affect functionality, but "
              "you are encourage to consider giving your experiment a more "
              "descriptive name.".format(key), file=sys.stderr)
    if not isinstance(entry, dict):
        raise ConfigError(msg, key, entry)
    if "designmatrix" not in entry:
        if "treatment" not in entry or "control" not in entry:
            raise ConfigError(msg, key, entry)
        for condition in "treatment control".split():
            if not isinstance(entry[condition], list):
                raise ConfigError("Expecting a list of {} samples.".format(condition), key, entry)
            for sample in entry[condition]:
                is_str(key, sample)
    else:
        is_file(key, entry["designmatrix"])


def is_experiments(key, entry):
    if not isinstance(entry, dict):
        raise ConfigError("Expecting an assignment of experiments to treatment and control samples or a design matrix.", key, entry)
    for experiment, entry in entry.items():
        is_experiment(experiment, entry)


config_structure = {
    "library": (is_file, False),
    "species": (is_str, True),
    "assembly": (is_str, True),
    "targets": {
        "genes": (is_bool, True),
        "controls": (is_file, False)
    },
    "sgrnas": {
        "trim-5": (is_str, True),
        "len": (is_str_or_int, True),
        "annotation": (is_file, False),
        "adapter": (is_str, False),
    },
    "samples": (is_samples, False),
    "correct_cnv": (is_bool, True),
    "cnv_norm": (is_file, False),
    "experiments": (is_experiments, True)
}


def _check_config(subconfig, substructure):
    for key, entry in substructure.items():
        if isinstance(entry, dict):
            test = None
            required = True
        else:
            test, required = entry
        if key not in subconfig:
            if required:
                raise ConfigError("Missing {} entry.".format(key))
        else:
            if test:
                test(key, subconfig[key])
            else:
                _check_config(subconfig[key], substructure[key])


def check_config(config):
    _check_config(config, config_structure)
