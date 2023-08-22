#!/usr/bin/env python3
import argparse
import gzip
import shutil
import os
import logging
from datetime import datetime
from itertools import chain
from pathlib import Path


logging.basicConfig(level=logging.INFO)
sampleid = lambda x: str(x.stem).split('_')[0]
repid = lambda x: str(x.stem).split('_')[2]
tmp_dir = Path('/', 'gpfs', 'gsfs12', 'users', 'NHLBI_IDSS', 'rawdata', 'metamorph-test', 'merged', 'rna', 'tmp')


def gunzip(from_gzfile, to_regfile):
    """
        Uncompress a file from a gzip compressed file to a regular uncompressed file

        @param from_gzfile <pathlib.Path> or <str>, gzipped from file
        @param to_regfile <pathlib.Path> or <str>, regular file output
    """
    with gzip.open(from_gzfile, 'rb') as f_in:
        with open(to_regfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return Path(to_regfile).exists()


def catfiles(outputfile, *files):
    """
        Concatenate multiple files into one

        @param outputfile <pathlib.Path> or <str>, file path to ouput concatination to
        @param *files <list>, list of valid file paths to concatinate
    """
    with open(outputfile, 'wb') as wfd:
        for f in files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return Path(outputfile).exists()


def main(args):
    """
        Collapse multiplex samples across lanes into single fastq.gz files. Optimally uncompressing first and recompressing.
        This script matches partner fastq.gz's by the a common sample identifier and lane identifier. Then uses lexicographical
        order to merge each consituent file into one.

        E.g.:
            Control.rna62_L001_R1.fastq.gz and Control.rna62_L004_R1.fastq.gz
            Control.rna62_L001_R2.fastq.gz and Control.rna62_L004_R2.fastq.gz
        File name template:
            <sample id>_<lane id>_<replicate id>.fastq.gz

        @param args.outputdir <pathlib.Path> output directory
        @param args.inputdir <pathlib.Path> input directory
    """

    args.inputdir = Path(args.inputdir).resolve()
    args.outputdir = Path(args.outputdir).resolve()
    
    if not args.outputdir.exists():
        Path.mkdir(args.outputdir)

    # pairing
    all_files = sorted(args.inputdir.glob('*.fastq.gz'))
    pairs = []
    for _file in all_files:
        if _file in list(chain.from_iterable(pairs)):
            continue
        sid = sampleid(_file)
        rid = repid(_file)
        partner = None
        for _partner_file in all_files:
            if _file == _partner_file: continue
            partner_sid = sampleid(_partner_file)
            partner_rid = repid(_partner_file)
            if partner_sid == sid and partner_rid == rid:
                if partner is not None:
                    logging.warn(
                        f'More than two files with the sample identifier: {partner_sid}\n'
                        f'--- Files: {_file}, {_partner_file}, {partner}'
                    )
                partner = _partner_file
        if partner is None:
            logging.warn(f'No partner match for file: {_file}')
            continue
        pairs.append(sorted([_file, partner], key=lambda x: x.stem))


    # decompressing-recompressing
    log = '/gpfs/gsfs12/users/NHLBI_IDSS/rawdata/metamorph-test/merged/rna/.collapse_log'
    for pair in pairs:
        final_path = Path(args.outputdir, sampleid(pair[0]) + '_' + repid(pair[0]) + '.gz')
        if final_path.exists():
            continue
        
        p0tmp = str(Path(tmp_dir, pair[0].name))
        p1tmp = str(Path(tmp_dir, pair[1].name))
        p_concat = str(Path(tmp_dir, sampleid(pair[0])))
        p_concat_gz = str(Path(tmp_dir, sampleid(pair[0]) + '_' + repid(pair[0]) + '.gz'))

        log_entry = r">> " + str(pair[0]) + r"\n" + \
                    r">> " + str(pair[1]) + r"\n" + \
                    r"\t - merged to: " + str(final_path) + r"\n\n"

        cmd = \
            f"gunzip -c {pair[0]} > {p0tmp} ; " + \
            f"gunzip -c {pair[1]} > {p1tmp} ; " + \
            f"cat {p0tmp} {p1tmp} > {p_concat} ; " + \
            f"gzip -c {p_concat} > {p_concat_gz} ; " + \
            f"mv {p_concat_gz} {str(final_path)} ; " + \
            f"rm {p0tmp} {p1tmp} {p_concat} ; " + \
            f"echo \"{log_entry}\" >>  {log}\n"
        
        with open('merge.swarm', 'a') as f:
            f.write(cmd)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputdir", type=Path, help="directory containing all fastq.qz files to be combined ")
    parser.add_argument("outputdir", type=Path, help="directory to copy combined fastq.gzs to")
    main(parser.parse_args())
