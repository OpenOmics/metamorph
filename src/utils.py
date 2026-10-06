#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import os, sys, hashlib
import subprocess, json
from csv import DictReader


def md5sum(filename, first_block_only = False, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks/chunks defined by the blocksize parameter. This is 
    a safer option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param first_block_only <bool>:
        Calculate md5 checksum of the first block/chunk only
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        if first_block_only:
            # Calculate MD5 of first block or chunck of file.
            # This is a useful heuristic for when potentially 
            # calculating an MD5 checksum of thousand or 
            # millions of file.
            hasher.update(buf)
            return hasher.hexdigest()
        while len(buf) > 0:
            # Calculate MD5 checksum of entire file
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


def longest_common_parent_path(data):
    data = [x if x[-1] == '/' else x + '/' for x in data]
    substrs = lambda x: {x[i:i+j] for i in range(len(x)) for j in range(len(x) - i + 1)}
    s = substrs(data[0])
    for val in data[1:]:
        s.intersection_update(substrs(val))
    longest = max(s, key=len)
    if longest.count('/') < 3:
        return None
    return os.path.dirname(longest)


def permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the user is authorized to access
    a file/directory. Checks for existence, readability, writability and executability via:
    os.F_OK (tests existence), os.R_OK (tests read), os.W_OK (tests write), os.X_OK (tests exec).
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        Returns abs path if it exists and permissions are correct
    """
    if not exists(path):
        parser.error("Path '{}' does not exists! Failed to provide vaild input.".format(path))
    if not os.access(path, *args, **kwargs):
        parser.error("Path '{}' exists, but cannot read path due to permissions!".format(path))

    return os.path.abspath(path)


def standard_input(parser, path, *args, **kwargs):
    """Checks for standard input when provided or permissions using permissions().
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        If path exists and user can read from location
    """
    # Checks for standard input
    if not sys.stdin.isatty():
        # Standard input provided, set path as an
        # empty string to prevent searching of '-'
        path = ''
        return path

    # Checks for positional arguments as paths
    path = permissions(parser, path, *args, **kwargs)

    return path


def exists(testpath):
    """Checks if file exists on the local filesystem.
    @param parser <argparse.ArgumentParser() object>:
        argparse parser object
    @param testpath <str>:
        Name of file/directory to check
    @return does_exist <boolean>:
        True when file/directory exists, False when file/directory does not exist
    """
    does_exist = True
    if not os.path.exists(testpath):
        does_exist = False # File or directory does not exist on the filesystem

    return does_exist


def ln(files, outdir):
    """Creates symlinks for files to an output directory.
    @param files list[<str>]:
        List of filenames
    @param outdir <str>:
        Destination or output directory to create symlinks
    """
    # Create symlinks for each file in the output directory
    for file in files:
        ln = os.path.join(outdir, os.path.basename(file))
        if not exists(ln):
            os.symlink(os.path.abspath(file), ln)


def which(cmd, path=None):
    """Checks if an executable is in $PATH
    @param cmd <str>:
        Name of executable to check
    @param path <list>:
        Optional list of PATHs to check [default: $PATH]
    @return <boolean>:
        True if exe in PATH, False if not in PATH
    """
    if path is None:
        path = os.environ["PATH"].split(os.pathsep)

    for prefix in path:
        filename = os.path.join(prefix, cmd)
        executable = os.access(filename, os.X_OK)
        is_not_directory = os.path.isfile(filename)
        if executable and is_not_directory:
            return True
    return False


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def require(cmds, suggestions, path=None):
    """Enforces an executable is in $PATH
    @param cmds list[<str>]:
        List of executable names to check
    @param suggestions list[<str>]:
        Name of module to suggest loading for a given index
        in param cmd.
    @param path list[<str>]]:
        Optional list of PATHs to check [default: $PATH]
    """
    error = False
    for i in range(len(cmds)):
        available = which(cmds[i])
        if not available:
            c = Colors
            error = True
            err("""\n{}{}Fatal: {} is not in $PATH and is required during runtime!{}
            └── Possible solution: please 'module load {}' and run again!""".format(
                c.bg_red, c.white, cmds[i], c.end, suggestions[i])
            )

    if error: fatal()

    return 


def safe_copy(source, target, resources = []):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            copytree(os.path.join(source, resource), destination)


def git_commit_hash(repo_path):
    """Gets the git commit hash of the repo.
    @param repo_path <str>:
        Path to git repo
    @return githash <str>:
        Latest git commit hash
    """
    try:
        githash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT, cwd = repo_path).strip().decode('utf-8')
        # Typecast to fix python3 TypeError (Object of type bytes is not JSON serializable)
        # subprocess.check_output() returns a byte string
        githash = str(githash)
    except Exception as e:
        # Github releases are missing the .git directory,
        # meaning you cannot get a commit hash, set the 
        # commit hash to indicate its from a GH release
        githash = 'github_release'
    return githash


def join_jsons(templates):
    """Joins multiple JSON files to into one data structure
    Used to join multiple template JSON files to create a global config dictionary.
    @params templates <list[str]>:
        List of template JSON files to join together
    @return aggregated <dict>:
        Dictionary containing the contents of all the input JSON files
    """
    # Get absolute PATH to templates in git repo
    repo_path = os.path.dirname(os.path.realpath(__file__))
    aggregated = {}

    for file in templates:
        with open(os.path.join(repo_path, file), 'r') as fh:
            aggregated.update(json.load(fh))

    return aggregated


def check_cache(parser, cache, *args, **kwargs):
    """Check if provided SINGULARITY_CACHE is valid. Singularity caches cannot be
    shared across users (and must be owned by the user). Singularity strictly enforces
    0700 user permission on on the cache directory and will return a non-zero exitcode.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param cache <str>:
        Singularity cache directory
    @return cache <str>:
        If singularity cache dir is valid
    """
    if not exists(cache):
        # Cache directory does not exist on filesystem
        os.makedirs(cache)
    elif os.path.isfile(cache):
        # Cache directory exists as file, raise error
        parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
        The provided --singularity-cache already exists on the filesystem as a file.
        Please run {} again with a different --singularity-cache location.
        """.format(sys.argv[0]))
    elif os.path.isdir(cache):
        # Provide cache exists as directory
        # Check that the user owns the child cache directory
        # May revert to os.getuid() if user id is not sufficent
        if exists(os.path.join(cache, 'cache')) and os.stat(os.path.join(cache, 'cache')).st_uid != os.getuid():
                # User does NOT own the cache directory, raise error
                parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
                The provided --singularity-cache already exists on the filesystem with a different owner.
                Singularity strictly enforces that the cache directory is not shared across users.
                Please run {} again with a different --singularity-cache location.
                """.format(sys.argv[0]))

    return cache


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary 
    """
    # Iterate over all values of given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over
            # all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type then
            # yield the value
            yield value


class Colors():
    """Class encoding for ANSI escape sequeces for styling terminal text.
    Any string that is formatting with these styles must be terminated with
    the escape sequence, i.e. `Colors.end`.
    """
    # Escape sequence
    end = '\33[0m'
    # Formatting options
    bold   = '\33[1m'
    italic = '\33[3m'
    url    = '\33[4m'
    blink  = '\33[5m'
    higlighted = '\33[7m'
    # Text Colors
    black  = '\33[30m'
    red    = '\33[31m'
    green  = '\33[32m'
    yellow = '\33[33m'
    blue   = '\33[34m'
    pink  = '\33[35m'
    cyan  = '\33[96m'
    white = '\33[37m'
    # Background fill colors
    bg_black  = '\33[40m'
    bg_red    = '\33[41m'
    bg_green  = '\33[42m'
    bg_yellow = '\33[43m'
    bg_blue   = '\33[44m'
    bg_pink  = '\33[45m'
    bg_cyan  = '\33[46m'
    bg_white = '\33[47m'


def hashed(l):
    """Returns an MD5 checksum for a list of strings. The list is sorted to
    ensure deterministic results prior to generating the MD5 checksum. This 
    function can be used to generate a batch id from a list of input files.
    It is worth noting that path should be removed prior to calculating the 
    checksum/hash.
    @Input:
        l list[<str>]: List of strings to hash
    @Output:
        h <str>: MD5 checksum of the sorted list of strings
    Example:
        $ echo -e '1\n2\n3' > tmp
        $ md5sum tmp
        # c0710d6b4f15dfa88f600b0e6b624077  tmp
        hashed([1,2,3])   # returns c0710d6b4f15dfa88f600b0e6b624077
    """
    # Sort list to ensure deterministic results
    l = sorted(l)
    # Convert everything to strings
    l = [str(s) for s in l]
    # Calculate an MD5 checksum of results
    h = hashlib.md5()
    # encode method ensure cross-compatiability 
    # across python2 and python3 
    h.update("{}\n".format("\n".join(l)).encode())
    h = h.hexdigest()

    return h


def valid_input(sheet):
    from argparse import ArgumentTypeError
    """
    Valid sample sheets should contain two columns: "DNA" and "RNA"

             _________________
             |  DNA  |  RNA  |
             |---------------| 
    pair1    |  path |  path |
    pair2    |  path |  path |
    """
    # check file permissions
    sheet = os.path.abspath(sheet)
    if not os.path.exists(sheet):
        raise ArgumentTypeError(f'Sample sheet path {sheet} does not exist!')
    if not os.access(sheet, os.R_OK):
        raise ArgumentTypeError(f"Path `{sheet}` exists, but cannot read path due to permissions!")

    # check format to make sure it's correct
    if sheet.endswith('.tsv') or sheet.endswith('.txt'):
        delim = '\t'
    elif sheet.endswith('.csv'):
        delim = ','

    rdr = DictReader(open(sheet, 'r'), delimiter=delim)

    if 'DNA' not in rdr.fieldnames:
        raise ArgumentTypeError("Sample sheet does not contain `DNA` column")
    if 'RNA' not in rdr.fieldnames:
        print("-- Running in DNA only mode --")
    else:
        print("-- Running in paired DNA & RNA mode --")

    data = [row for row in rdr]
    RNA_included = False
    seen = {}
    for i, row in enumerate(data):
        row['DNA'] = os.path.abspath(row['DNA'])
        if not os.path.exists(row['DNA']):
            raise ArgumentTypeError(f"Sample sheet path `{row['DNA']}` does not exist")
        if 'RNA' in row and not row['RNA'] in ('', None, 'None'):
            RNA_included = True
            row['RNA'] = os.path.abspath(row['RNA'])
            if not os.path.exists(row['RNA']):
                raise ArgumentTypeError(f"Sample sheet path `{row['RNA']}` does not exist")

        # every file in the sheet must be unique -- catches copy-paste
        # mistakes like the same fastq reused as R1 for one sample and
        # R2 for another (or DNA for one sample and RNA for another)
        for col in ('DNA', 'RNA'):
            if col not in row or row[col] in ('', None, 'None'):
                continue
            path = row[col]
            if path in seen:
                raise ArgumentTypeError(
                    f"Sample sheet path `{path}` is used more than once: "
                    f"as {seen[path]} and again as {col} in sheet row {i + 1}. "
                    f"Each file in the sample sheet must be unique."
                )
            seen[path] = f"{col} in sheet row {i + 1}"

    # Structural pairing checks: every sample must resolve to exactly
    # one R1 and one R2 file, and no two different input files may
    # collide on the same renamed symlink name once staged. metamorph
    # stages inputs via src.run.sym_safe(), which calls src.run.rename()
    # to rewrite each file's basename to <sample>_R[12].fastq.gz before
    # symlinking it into <output>/dna or <output>/rna. Reuse that exact
    # transform here (lazy import to avoid a circular import with
    # src.run, which imports from this module) so sheet validation can
    # never drift from actual staging behavior.
    from .run import rename, FASTQ_R1_POSTFIX, FASTQ_R2_POSTFIX

    for col in ('DNA', 'RNA'):
        entries = [
            (i, row[col]) for i, row in enumerate(data)
            if col in row and row[col] not in ('', None, 'None')
        ]
        if not entries:
            continue

        staged = {}   # staged symlink name -> (row index, original path)
        samples = {}  # sample id -> {'R1': [(row, path), ...], 'R2': [...]}
        for i, path in entries:
            try:
                staged_name = rename(os.path.basename(path))
            except NameError as e:
                raise ArgumentTypeError(str(e))

            if staged_name in staged:
                other_row, other_path = staged[staged_name]
                raise ArgumentTypeError(
                    f"Sample sheet paths `{other_path}` ({col} row {other_row + 1}) "
                    f"and `{path}` ({col} row {i + 1}) both stage to the same "
                    f"renamed filename `{staged_name}` when symlinked into the "
                    f"pipeline's {col.lower()} input directory. Each input file "
                    f"must produce a unique staged filename."
                )
            staged[staged_name] = (i, path)

            if staged_name.endswith(FASTQ_R1_POSTFIX):
                mate, sample_id = 'R1', staged_name[:-len(FASTQ_R1_POSTFIX)]
            elif staged_name.endswith(FASTQ_R2_POSTFIX):
                mate, sample_id = 'R2', staged_name[:-len(FASTQ_R2_POSTFIX)]
            else:
                # Unreachable: rename() only ever returns names ending
                # in FASTQ_R1_POSTFIX or FASTQ_R2_POSTFIX (or raises).
                continue
            samples.setdefault(sample_id, {'R1': [], 'R2': []})[mate].append((i, path))

        # Skip pairing checks for columns that are entirely single-end
        # (no R2 present at all), consistent with src.run.get_nends().
        if not any(mates['R2'] for mates in samples.values()):
            continue

        for sample_id, mates in samples.items():
            r1s, r2s = mates['R1'], mates['R2']

            # Each mate list can only ever hold 0 or 1 entries here: two
            # inputs sharing a sample id and mate would necessarily also
            # share the exact same staged name, and that collision is
            # already caught above before entries are ever grouped here.
            if r1s and r2s:
                r1_row, r1_path = r1s[0]
                r2_row, r2_path = r2s[0]
                if os.path.realpath(r1_path) == os.path.realpath(r2_path):
                    raise ArgumentTypeError(
                        f"Sample `{sample_id}` ({col}): R1 path `{r1_path}` "
                        f"(row {r1_row + 1}) and R2 path `{r2_path}` "
                        f"(row {r2_row + 1}) resolve to the same file on disk. "
                        f"R1 and R2 must be different files."
                    )
                continue

            missing_mate = 'R1' if not r1s else 'R2'
            present_row, present_path = (r1s or r2s)[0]
            raise ArgumentTypeError(
                f"Sample `{sample_id}` ({col}) is missing its {missing_mate} file "
                f"(only found `{present_path}`, row {present_row + 1}). Every "
                f"paired-end sample must have exactly one R1 and one R2 file."
            )

    return data, RNA_included


def valid_trigger(trigger_given):
    from argparse import ArgumentTypeError
    snk_triggers = ('mtime', 'code', 'input', 'params', 'software-env')
    trigger_list =  trigger_given.split(',')
    invalid = [t for t in trigger_list if t not in snk_triggers]
    if invalid:
            raise ArgumentTypeError('Invalid trigger selected please only use one of: ' + ', '.join(snk_triggers))
    return trigger_given


def valid_host_genome(genome_given):
    """Validates a --host-genome value: must be an existing file (a FASTA,
    optionally gzipped, of a genome to screen out during dehosting in
    addition to the default host reference). Returns the resolved absolute
    path so downstream bind-path/Snakemake logic never has to re-resolve it.
    """
    from argparse import ArgumentTypeError
    resolved = os.path.realpath(os.path.abspath(os.path.expanduser(genome_given)))
    if not exists(resolved):
        raise ArgumentTypeError(f'--host-genome file does not exist: {genome_given}')
    if not os.path.isfile(resolved):
        raise ArgumentTypeError(f'--host-genome path is not a file: {genome_given}')
    return resolved


if __name__ == '__main__':
    # Calculate MD5 checksum of entire file 
    print('{}  {}'.format(md5sum(sys.argv[0]), sys.argv[0]))
    # Calcualte MD5 cehcksum of 512 byte chunck of file,
    # which is similar to following unix command: 
    # dd if=utils.py bs=512 count=1 2>/dev/null | md5sum 
    print('{}  {}'.format(md5sum(sys.argv[0], first_block_only = True, blocksize = 512), sys.argv[0]))
