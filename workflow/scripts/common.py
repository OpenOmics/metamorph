#!/usr/bin/env python
# ~~~~~~~~~~
# Common utility functions used throughout the workflow
# ~~~~~~~~~~
from os.path import join

def get_dna_sid(config, rna_sid):
    if 'sample_map' not in config:
        return None
    if not list_bool(config.get("rna", 'false')):
        return None
    sample_map = config['sample_map']
    return sample_map[rna_sid]

def get_paired_dna(config, rna_sid):
    if 'sample_map' not in config:
        return None
    if not list_bool(config.get("rna", 'false')):
        return None
    sample_map = config['sample_map']
    dna_sid = sample_map[rna_sid]
    
    # dna directories
    top_refine_dir = join(config["project"]["workpath"], config['project']['id'], "metawrap_bin_refine")
    top_mags_dir = join(config["project"]["workpath"], config['project']['id'], "mags")
    # construct necessary linked dna artifacts
    derep_reqs = join(top_refine_dir, "dRep", "dereplicated_genomes")
    mag_index = join(top_mags_dir, dna_sid, "index")
    return [derep_reqs, mag_index]


def provided(samplelist, condition):
    """
        Determines if optional rules should run. If an empty list is provided to rule all,
        snakemake will not try to generate that set of target files. If a given condition
        is not met (i.e. False) then it will not try to run that rule.
    """
    if not condition:
        # If condition is False, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist


def references(config, reflist):
    """
        Checks if a set of required reference files were provided. Some rules depend
        on a set of required reference files that may only exist for specific reference
        genomes. An example of this would be blasklists arriba. The blacklist are manually
        curated and only exist for a few reference genomes (mm10, hg38, hg19).
        If one of the required reference files does not exist, then it will return
        an empty list.
    """
    _all = True
    for ref in reflist:
        try: tmp = config['references'][ref]
        # Check if ref exists in config
        except KeyError:
            _all = False
            break
        # Check if ref is empty key string
        if not tmp: _all = False
    return _all


def allocated(resource, rule, lookup, default="__default__"):
    """
        Pulls resource information for a given rule. If a rule does not have any information 
        for a given resource type, then it will pull from the default. Information is pulled from
        definitions in the cluster.json (which is used a job submission). This ensures that any 
        resources used at runtime mirror the resources that were allocated.
        :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
        :param rule <str>: rule to lookup its information
        :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
        :param default <str>: default information to use if rule information cannot be found
        :return allocation <str>: 
            allocation information for a given resource type for a given rule
    """
    try: 
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]
    
    return allocation


def str_bool(s):
    """
        Deprecated from stdlib in 3.10, see: https://peps.python.org/pep-0632/

        Converts a string to boolean. It is dangerous to try to
        typecast a string into a boolean value using the built-in 
        `bool()` function. This function avoids any issues that can
        arise when using `bool()`. 
        Example:
            boolean('True') returns True
            boolean('False') returns False
            boolean('asdas') raises TypeError
    """
    if not isinstance(s, str):
        s = str(s)
    val = s.lower()
    if val in ['true', '1', 'y', 'yes']:
        return True
    elif val in ['false', '0', 'n', 'no', '']:
        return False
    else:
        # Provided value could not be
        # type casted into a boolean
        raise TypeError('Fatal: cannot type cast {} into a boolean'.format(val))


def list_bool(l):
    # some lists are strings of "None" if no files input
    # type checking instead of string manip
    if isinstance(l, list):
        if len(l) > 0:
            return True
    return False
