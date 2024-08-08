__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import argparse
import logging
import sys
import os
import shutil
import jinja2

import yaml

from mageck_vispr.version import __version__
from mageck_vispr import annotation


def init_workflow(directory, reads, keep_config=False):
    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    def guess_sample_name(f):
        if f.endswith(".fastq"):
            return os.path.splitext(os.path.basename(f))[0]
        elif f.endswith(".fastq.gz"):
            return os.path.splitext(
                os.path.splitext(os.path.basename(f))[0])[0]
        else:
            logging.error("Reads must be given as .fastq or .fastq.gz files.")

    def get_resource(f):
        return os.path.join(os.path.dirname(__file__), f)

    def get_target(f):
        return os.path.join(directory, f)

    def backup(f):
        target = get_target(f)
        if os.path.exists(target):
            shutil.copy(target, target + ".old")

    def install(f, backup=True):
        source = get_resource(f)
        target = get_target(f)
        shutil.copy(source, target)

    backup("Snakefile")
    install("Snakefile")
    install("README.txt")
    if not keep_config:
        backup("config.yaml")
        with open(get_resource("config.yaml")) as f:
            config_template = jinja2.Template(f.read(),
                                              trim_blocks=True,
                                              lstrip_blocks=True)
            samples = None
            if reads is not None:
                samples = {
                    guess_sample_name(f): os.path.relpath(f, directory)
                    for f in reads
                }
            #print('samples:'+str(samples))
            render_results=config_template.render(samples=samples)
            #print(render_results)
            with open(get_target("config.yaml"), "w") as out:
                out.write(render_results)


def annotate_library(args):
    library=args.library
    #assembly=args.assembly
    #sgrna_len=args.sgrna_len
    #annotation_table=args.annotation_table

    # a = annotation.Annotator(library, annotation_table)
    a = annotation.Annotator(library )

    # if --bedvalue is specified, read the bedvalue file
    a.add_value_frame(args)

    a.annotate(args)


def main():
    # create arg parser
    parser = argparse.ArgumentParser(
        "MAGeCK-VISPR is a comprehensive quality control, analysis and "
        "visualization pipeline for CRISPR/Cas9 screens.")
    parser.add_argument("--version",
                        action="store_true",
                        help="Print version info.")
    subparsers = parser.add_subparsers(dest="subcommand")

    workflow = subparsers.add_parser(
        "init",
        help="Initialize the MAGeCK/VISPR workflow "
        "in a given directory. This will "
        "install a Snakefile, a README and a "
        "config file in this directory. "
        "Configure the config file according "
        "to your needs, and run the workflow "
        "with Snakemake "
        "(https://bitbucket.org/johanneskoester/snakemake).")
    workflow.add_argument("directory",
                          help="Path to the directory where the "
                          "workflow shall be initialized.")
    workflow.add_argument("--reads",
                          nargs="+",
                          help="Paths to FastQ files with reads that shall be "
                          "added to the config file. You can edit the sample "
                          "sample names and assignment to experiments "
                          "in the config file.")
    workflow.add_argument("--keep-config", action="store_true",
                          help="Keep existing config file.")

    annotate = subparsers.add_parser(
        "annotate-library",
        help="Annotate an sgRNA library design with information about "
        "sgRNA position and predicted efficiency. Annotation is printed in "
        "BED format.")
    annotate.add_argument(
        "library",
        help="Path to sgRNA library design file (comma separated, columns "
        "identifier, sequence, gene).")
    annotate.add_argument("--sgrna-len",
                          #type=int,
                          choices=['19', '20', 'AUTO'],
                          help="Length of sgrnas in library file.")
    annotate.add_argument("--assembly",
                          choices=["mm10", "mm9", "hg38", "hg19"],
                          help="Assembly to use.")
    annotate.add_argument(
        "--bedvalue",
        help="Instead of providing an efficiency value in the output bed file, "
            "use the values provided in a given txt file."
            "The txt file must be tab separated, with header."
            "The first column must be sgRNA ID that matches the identifier in sgRNA library design file.")
    annotate.add_argument(
        "--bedvalue-column",
        help="Provide a column name in the file in --bedvalue option as the column to fill in. "
            "For example, the 'LFC' column in sgrna_summary.txt in MAGeCK RRA represents the log fold change value. ")
    annotate.add_argument(
        "--annotation-table-folder",
        help="After specifying the sgrna length and assembly, instead of downloading directly from bitbucket, search in the folder for corresponding annotation table.")
    annotate.add_argument(
        "--annotation-table",
        help="As an alternative to specifying the sgrna length and assembly, "
        "a path to an annotation table can be provided "
        "(tab separated, no header; with columns chromosome, "
        "start, end, gene, score, strand, sequence). This can also be a URL. "
        "See https://bitbucket.org/liulab/mageck-vispr/downloads for precomputed tables.")

    logging.basicConfig(format="%(message)s",
                        level=logging.INFO,
                        stream=sys.stderr)

    args = parser.parse_args()

    if args.version:
        print(__version__)
        exit(0)
    elif args.subcommand == "init":
        init_workflow(args.directory, args.reads, keep_config=args.keep_config)
    elif args.subcommand == "annotate-library":
        if args.assembly or args.annotation_table:
            annotate_library(args)
        else:
            parser.print_help()
            print("Error: need to specify one of the following: path to an annotation table (--annotation-table); or  assembly (--assembly).")
            exit(1)
    else:
        parser.print_help()
        exit(1)
    exit(0)
