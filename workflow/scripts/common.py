import re

def ingest_sample_files(dna, rna):
    """
        :param dna <list>: 
        :param rna <list>: 
        :return parsed_info <tuple>: return tuple for each sample consisting of it's
            condition <str>
            name <str>
            rep id <int>
            dna reads <str>
            rna reads <str> or None (if only meatgenomic analysis)

        Files processed here must adhere to the convention:
            <condition> . <name> _ <replicate> .fastq.gz
        e.g.:
            Control.dna97_R2.fastq.gz
    """
    pattern = r"(^.+?(?=\.))\.(.+?(?=\_))\_(.+?(?=\.))"

    out = []
    for this_dna, this_rna in zip(dna, rna):
        matches = re.finditer(pattern, this_dna, re.MULTILINE)
        
    return out


def str_bool(s):
    """
        Converts a string to boolean. It is dangerous to try to
        typecast a string into a boolean value using the built-in 
        `bool()` function. This function avoids any issues that can
        arise when using `bool()`. 
        Example:
            boolean('True') returns True
            boolean('False') returns False
            boolean('asdas') raises TypeError
    """
    val = s.lower()
    if val in ['true', '1', 'y', 'yes']:
        return True
    elif val in ['false', '0', 'n', 'no', '']:
        return False
    else:
        # Provided value could not be
        # type casted into a boolean
        raise TypeError('Fatal: cannot type cast {} into a boolean'.format(val))