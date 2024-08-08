__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import os
import glob
import shutil
import re    #[cuiyb]++
from collections import defaultdict

import yaml

from mageck_vispr.check_config import check_config, ConfigError


COMBAT_SCRIPT_PATH = os.path.join(os.path.dirname(__file__), "combat.R")


def get_norm_method(config):
    if "norm_method" in config:
        return config["norm_method"]
    elif "batchmatrix" in config:
        return "none"
    else:
        return "median"


def get_counts(config, normalized=False):
    if "batchmatrix" in config:
        return "results/count/all.count.batchcorrected.txt"
    else:
        suffix = "_normalized" if normalized else ""
        return config.get("counts", "results/count/all.count{}.txt".format(suffix))

def get_sample_name(config):
    """
    Get sample name from count tables
    """
    if "day0label" in config:
        day0sample=config["day0label"]
    else:
        day0sample=None
    day0samplelist=day0sample.split(',')
    if "samples" in config:
        return [fx for fx in config["samples"].keys() if fx not in day0samplelist]
    else:
        # count table provided
        count_table_file=get_counts(config)
        c_f=open(count_table_file)
        fields=c_f.readline().strip().split()
        c_f.close()
        return [fx for fx in fields[2:] if fx not in day0samplelist]



def annotation_available(config):
    #if 'annotate-sgrna' not in config['sgrnas']:
    #    return False
    if config['sgrnas']['annotate-sgrna'] == False:
        return False
    return "library" in config and config["assembly"] in ["mm9", "mm10", "hg19", "hg38"]  


def design_available(config):
    """
    Returns true only when it's an MLE experiment and a real design matrix (not /dev/null) is provided
    """
    if re.search('designmatrix', str(config["experiments"])) and not re.search('\/dev\/null', str(config["experiments"])):
        if "day0label" in config:
            raise ConfigError("Users can either specify a day0 label (with an empty design matrix), or a design matrix file, but can't specify the both.")
        return True
    else:  # no design matrix assigned
        if "day0label" not in config and not (re.search('treatment', str(config["experiments"])) and re.search('control', str(config["experiments"]))):
            raise ConfigError("No design matrix, day0 label, or treatment and control assigned.")
        return False





def efficiency_estimation_available(config, experiment):
    return "designmatrix" in config["experiments"][experiment]


def postprocess_config(config):
    check_config(config)
    config["replicates"] = {}
    if "samples" in config:
        if not "library" in config:
            raise ConfigError("A library file has to be specified in the configuration.")
        if "counts" in config:    # [cuiyb]++ check the available of "samples" and "counts"
            raise ConfigError("\"samples\" and \"counts\" cannot be set as valid at the same time.")
        for sample in config["samples"]:
            replicates = config["samples"][sample]
            if not isinstance(replicates, list):
                replicates = [replicates]
            config["samples"][sample] = []
            for i, replicate in enumerate(replicates):
                name = "{}_{}".format(sample, i)
                config["replicates"][name] = replicate
                config["samples"][sample].append(name)
    # [wubing] Enable mageck compatible to pair-end sequencing data
    config["paired_rep"] = {}
    if "paired" in config:
        for sample in config["paired"]:
            paired_rep = config["paired"][sample]
            if not isinstance(paired_rep, list):
                paired_rep = [paired_rep]
            config["paired"][sample] = []
            for i, replicate in enumerate(paired_rep):
                name = "{}_{}".format(sample, i)
                config["paired_rep"][name] = replicate
                config["paired"][sample].append(name)
    # add default options if the options are missing
    # this allows MAGeCK-VISPR to work on earlier versions of config.yaml
    if 'annotate-sgrna' not in config['sgrnas']:
        config['sgrnas']['annotate-sgrna']=False
    if 'annotate-sgrna-efficiency' not in config['sgrnas']:
        config['sgrnas']['annotate-sgrna-efficiency']=False
    if 'correct_cnv' not in config:
        config['correct_cnv']=False
    if 'cnv_norm' not in config:
        config['cnv_norm']='/dev/null'



def get_fastq(replicate, config):
    if "adapter" in config["sgrnas"]:
        return "results/trimmed_reads/{}.fastq".format(replicate)
    return config["replicates"][replicate]

def need_run_rra_in_mle(wildcards, config):
    # return True
    if "day0label" not in config: # not specifying day0label
        return False
    if "designmatrix" not in config["experiments"][wildcards.experiment]: # specify day0label but it's an RRA experiment
        return False
    return True

def rra_treatment_string(wildcards, config):
    if "day0label" in config:
        return ""
    if "treatment" not in config["experiments"][wildcards.experiment]:
        return ""
    return "--treatment-id "+",".join(config["experiments"][wildcards.experiment]["treatment"])


def rra_control_string(wildcards, config):
    if "day0label" in config:
        return ""
    if "control" not in config["experiments"][wildcards.experiment]:
        return ""
    return "--control-id "+",".join(config["experiments"][wildcards.experiment]["control"])

def need_annotate_bed_with_lfc(config):
    if annotation_available(config) and not design_available(config): # only activates when annotation is enabled and (either day0 is provided or an RRA experiment). Note: MLE will not generate individual sgRNA information.
        return True
    else:
        return False




def vispr_config(input, output, wildcards, config):
    relpath = lambda path: os.path.relpath(path, "results")
    copy = lambda path: shutil.copy(path, "results")
    vispr_config = {
        "experiment": wildcards.experiment,
        "species": config["species"],
        "assembly": config["assembly"],
        "targets": {
            "results": relpath(input.results),
            "genes": config["targets"]["genes"]
        },
        "sgrnas": {
            "counts": relpath(input.counts)
        }
    }
    #if "fastqc" in input.keys():
    if "fastqc" in input.keys() and len(input["fastqc"])>0 :
        samples = {
            rep: sample
            for sample, replicates in config["samples"].items()
            for rep in replicates
        }
        qc = defaultdict(list)
        for replicate, fastqc in zip(config["replicates"], input.fastqc):
            qc[samples[replicate]].extend(
                relpath(os.path.join(data, "fastqc_data.txt"))
                for data in sorted(glob.iglob("{}/*_fastqc".format(fastqc))))
        vispr_config["fastqc"] = dict(qc)
    #if "pairedfastqc" in input.keys():
    if "pairedfastqc" in input.keys() and len(input["pairedfastqc"])>0 :
        samples = {
            rep: sample
            for sample, replicates in config["paired"].items()
            for rep in replicates
        }
        qc = defaultdict(list)
        for replicate, fastqc in zip(config["paired_rep"], input.pairedfastqc):
            qc[samples[replicate]].extend(
                relpath(os.path.join(data, "fastqc_data.txt"))
                for data in sorted(glob.iglob("{}/*_fastqc".format(fastqc))))
        vispr_config["pairedfastqc"] = dict(qc)
    if "mapstats" in input.keys() and len(input["mapstats"])>0 :
        vispr_config["sgrnas"]["mapstats"] = relpath(input.mapstats)
    if "controls" in config["targets"]:
        copy(config["targets"]["controls"])
        vispr_config["targets"]["controls"] = os.path.basename(
            config["targets"]["controls"])
    if annotation_available(config):
        copy("annotation/sgrnas.bed")
        vispr_config["sgrnas"]["annotation"] = "sgrnas.bed"
    if efficiency_estimation_available(config, wildcards.experiment):
        vispr_config["sgrnas"]["results"] = relpath(input.sgrna_results)
    with open(output[0], "w") as f:
        yaml.dump(vispr_config, f, default_flow_style=False)
