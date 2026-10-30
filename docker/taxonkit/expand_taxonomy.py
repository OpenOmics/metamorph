#!/usr/bin/env python3
"""Reconcile hierarchical abundance CSVs, optionally fill ranks with TaxonKit,
and aggregate exclusive contributions using phyloseq. See README.md.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import sys
import numpy as np

RANKS = 'kpcofgst'
NAMES = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
RESERVED = 'UNRESOLVED_BY_WORKFLOW'


def parse_path(path):
    if path in ('UNCLASSIFIED', 'UNKNOWN'):
        return {0: 'k__' + RESERVED + '_global'}
    result = {}
    previous = -1
    for token in path.split('|'):
        if len(token) < 3 or token[1:3] != '__' or token[0] not in RANKS:
            raise ValueError(f'Unsupported taxonomy token: {token!r} in {path!r}')
        rank = RANKS.index(token[0])
        if rank <= previous:
            raise ValueError(f'Ranks duplicated or out of order: {path}')
        if RESERVED in token:
            raise ValueError(f'Reserved placeholder string in input: {path}')
        result[rank] = token
        previous = rank
    return result


def canonicalize(path, looked_up=None):
    parts = parse_path(path)
    last = max(parts)
    for i in range(last + 1):
        if i not in parts or len(parts[i]) == 3:
            name = looked_up[i] if looked_up and i < len(looked_up) else ''
            if name and name != 'Unassigned':
                parts[i] = RANKS[i] + '__' + name.replace(' ', '_')
            else:
                parts[i] = RANKS[i] + '__' + RESERVED + '_missing_rank'
    return '|'.join(parts[i] for i in range(last + 1))


def taxonkit_lookup(mapping, data_dir, scratch):
    """Map EXACT supplied clade paths to explicit NCBI taxids; no name guessing."""
    if not mapping:
        return {}
    executable = shutil.which('taxonkit')
    if not executable:
        raise ValueError('TaxonKit is required when --taxid-map is supplied')
    with open(mapping, newline='') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != ['clade_name', 'taxid']:
            raise ValueError('Taxid map must have exactly clade_name,taxid columns')
        pairs = list(reader)
    if len({r['clade_name'] for r in pairs}) != len(pairs):
        raise ValueError('Duplicate clade_name in taxid map')
    if not pairs:
        raise ValueError('Taxid map is empty')
    ids = sorted({r['taxid'] for r in pairs})
    if any(not x.isdigit() or int(x) <= 0 for x in ids):
        raise ValueError('Taxids must be positive integer NCBI taxonomy identifiers')
    infile, outfile = scratch / 'taxids.txt', scratch / 'taxonkit.tsv'
    infile.write_text('\n'.join(ids) + '\n')
    fmt = '{domain|acellular root|superkingdom};{phylum};{class};{order};{family};{genus};{species}'
    cmd = [executable, 'reformat2', '-I', '1', '-r', 'Unassigned', '-f', fmt]
    if data_dir:
        cmd += ['--data-dir', data_dir]
    cmd += [str(infile)]
    with outfile.open('w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    found = {}
    for row in csv.reader(outfile.open(), delimiter='\t'):
        if len(row) != 2 or len(row[1].split(';')) != 7:
            raise ValueError('Unexpected TaxonKit output; reformat2 support is required')
        values = row[1].split(';')
        if not any(x and x != 'Unassigned' for x in values):
            raise ValueError(f'TaxonKit could not resolve taxid {row[0]}')
        if any('|' in x for x in values):
            raise ValueError('TaxonKit returned a name containing a reserved delimiter')
        found[row[0]] = values
    if set(found) != set(ids):
        raise ValueError('TaxonKit did not return every requested taxid')
    return {r['clade_name']: found[r['taxid']] for r in pairs}


def build_tree(paths):
    nodes = set()
    for path in paths:
        tokens = path.split('|')
        nodes.update('|'.join(tokens[:i]) for i in range(1, len(tokens) + 1))
    nodes = sorted(nodes, key=lambda x: (x.count('|'), x))
    index = {p: i for i, p in enumerate(nodes)}
    depths = np.array([p.count('|') for p in nodes])
    parent = np.array([index[p.rsplit('|', 1)[0]] if '|' in p else -1 for p in nodes])
    source = np.array([index[p] for p in paths])
    return nodes, depths, parent, source


def reconcile(values, nodes, depths, parent, source, atol, rtol):
    """Return nonoverlapping node contributions, resolved totals, and diagnostics.
    Explicit parents are authoritative except small child overruns within tolerance:
    those parents are raised to their child sum and the adjustment is logged.
    """
    n, width = len(nodes), values.shape[1]
    total = np.zeros((n, width))
    childsum = np.zeros_like(total)
    residual = np.zeros_like(total)
    supplied = np.zeros(n, dtype=bool)
    total[source] = values
    supplied[source] = True
    issues = []
    for depth in range(int(depths.max()), -1, -1):
        idx = np.flatnonzero(depths == depth)
        for i in idx:
            if supplied[i]:
                delta = total[i] - childsum[i]
                tolerance = atol + rtol * np.maximum(total[i], childsum[i])
                for j in np.flatnonzero(delta < 0):
                    kind = 'ERROR_children_exceed_parent' if -delta[j] > tolerance[j] else 'rounding_parent_raised'
                    issues.append((i, int(j), kind, float(total[i, j]), float(childsum[i, j]), float(delta[j])))
                residual[i] = np.maximum(delta, 0)
                total[i] = np.maximum(total[i], childsum[i])
            else:
                total[i] = childsum[i]
        p = parent[idx]
        valid = p >= 0
        np.add.at(childsum, p[valid], total[idx[valid]])
    # Log unresolved parent mass only for supplied nodes with children.
    internal = np.zeros(n, dtype=bool)
    internal[parent[parent >= 0]] = True
    for i in np.flatnonzero(supplied & internal):
        for j in np.flatnonzero(residual[i] > atol + rtol * total[i]):
            issues.append((int(i), int(j), 'unresolved_parent_mass', float(total[i, j]), float(childsum[i, j]), float(residual[i, j])))
    return residual, total, issues


def leaf_lineage(node):
    parts = node.split('|')[:7]  # SGB/strain contributions collapse to their species.
    while len(parts) < 7:
        parts.append(RANKS[len(parts)] + '__' + RESERVED + '_below_' + RANKS[len(parts)-1])
    return tuple(parts)


def exclusive_species(residual, nodes):
    active = np.flatnonzero(np.any(residual > 0, axis=1))
    groups = {}
    for i in active:
        key = leaf_lineage(nodes[i])
        if key not in groups:
            groups[key] = residual[i].copy()
        else:
            groups[key] += residual[i]
    keys = sorted(groups)
    if not keys:
        raise ValueError('No positive abundance in sample batch')
    return keys, np.stack([groups[k] for k in keys])


def run(args):
    if (args.batch_size < 1 or not np.isfinite(args.atol) or not np.isfinite(args.rtol)
            or args.atol < 0 or args.rtol < 0):
        raise ValueError('Batch size must be positive; tolerances must be nonnegative')
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be empty (prevents mixing runs)')
    rscript = shutil.which('Rscript')
    if not rscript:
        raise ValueError('Rscript is required; install R, phyloseq, and data.table')
    worker = Path(__file__).with_name('aggregate_phyloseq.R')
    if not worker.is_file():
        raise ValueError('Keep aggregate_phyloseq.R beside this Python script')
    subprocess.run([rscript, str(worker), '--check'], check=True)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.temp_dir) as temp:
        temp = Path(temp)
        lookup = taxonkit_lookup(args.taxid_map, args.taxonkit_data, temp)
        paths, originals, samples = [], [], None
        with open(args.input, newline='', encoding='utf-8-sig') as f, (temp / 'values.bin').open('wb') as binary:
            reader = csv.reader(line for line in f if not line.startswith('#'))
            header = next(reader)
            if not header or header[0] != 'clade_name' or len(header) < 2:
                raise ValueError('Expected CSV header clade_name,sample1,sample2,...')
            samples = header[1:]
            if (len(set(samples)) != len(samples) or any(not s for s in samples)
                    or set(samples) & set(NAMES + ['clade_name'])):
                raise ValueError('Sample names must be nonempty, unique, and different from taxonomy column names')
            for line, row in enumerate(reader, 2):
                if len(row) != len(header):
                    raise ValueError(f'Wrong number of columns on data line {line}')
                original = row[0].strip()
                path = canonicalize(original, lookup.get(original))
                values = np.asarray(row[1:], dtype=np.float64)
                if not np.all(np.isfinite(values)) or np.any(values < 0):
                    raise ValueError(f'Nonfinite or negative abundance on line {line}')
                binary.write(values.tobytes())
                paths.append(path)
                originals.append(original)
        if not paths:
            raise ValueError('No abundance rows')
        if len(set(paths)) != len(paths):
            raise ValueError('Duplicate clade paths after canonicalization; resolve rather than summing parent rows')
        unused = set(lookup) - set(originals)
        if unused:
            raise ValueError(f'Taxid map contains {len(unused)} paths absent from input')
        np_values = np.memmap(temp / 'values.bin', dtype=np.float64, mode='r', shape=(len(paths), len(samples)))
        nodes, depths, parent, source = build_tree(paths)
        with (out / 'lineage_changes.tsv').open('w', newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['input_clade', 'canonical_clade'])
            w.writerows((a, b) for a, b in zip(originals, paths) if a != b)
        with (out / 'nodes.tsv').open('w', newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['clade_name', 'origin'])
            supplied_paths = set(paths)
            w.writerows((p, 'supplied' if p in supplied_paths else 'inferred_ancestor') for p in nodes)
        checks_file = (out / 'hierarchy_checks.tsv').open('w', newline='')
        totals_file = (out / 'totals.tsv').open('w', newline='')
        input_rank_file = (out / 'input_rank_totals.tsv').open('w', newline='')
        input_rank_writer = csv.writer(input_rank_file, delimiter='\t')
        input_rank_writer.writerow(['sample', 'rank', 'supplied_row_count', 'supplied_rank_total'])
        check_writer = csv.writer(checks_file, delimiter='\t')
        total_writer = csv.writer(totals_file, delimiter='\t')
        check_writer.writerow(['sample', 'clade_name', 'status', 'parent_value', 'child_sum', 'difference'])
        total_writer.writerow(['sample', 'observed_total', 'multiplier', 'output_total', 'normalization'])
        manifest = {'input': str(Path(args.input).resolve()), 'samples': len(samples), 'input_rows': len(paths),
                    'tree_nodes': len(nodes), 'normalization': args.normalize, 'atol': args.atol, 'rtol': args.rtol,
                    'taxid_map': args.taxid_map, 'taxonkit_data': args.taxonkit_data,
                    'python': sys.version, 'numpy': np.__version__, 'status': 'incomplete'}
        (out / 'run.json').write_text(json.dumps(manifest, indent=2))
        try:
            for start in range(0, len(samples), args.batch_size):
                batch_samples = samples[start:start + args.batch_size]
                raw = np.array(np_values[:, start:start + args.batch_size])
                for d, rank_name in enumerate(NAMES + ['Strain_or_SGB']):
                    select = depths[source] == d
                    for j, name in enumerate(batch_samples):
                        input_rank_writer.writerow([name, rank_name, int(select.sum()), float(raw[select, j].sum())])
                residual, totals, issues = reconcile(raw, nodes, depths, parent, source, args.atol, args.rtol)
                for i, j, kind, p, c, d in issues:
                    check_writer.writerow([batch_samples[j], nodes[i], kind, p, c, d])
                checks_file.flush()
                if any(x[2].startswith('ERROR') for x in issues):
                    raise ValueError('Children exceed parents; see hierarchy_checks.tsv. No valid completed run produced.')
                observed = totals[parent < 0].sum(axis=0)
                if np.any(observed <= 0):
                    raise ValueError('Zero-total sample: ' + ', '.join(np.array(batch_samples)[observed <= 0]))
                multiplier = 100 / observed if args.normalize == 'observed' else np.ones(len(observed))
                residual *= multiplier
                keys, abundances = exclusive_species(residual, nodes)
                for j, name in enumerate(batch_samples):
                    total_writer.writerow([name, observed[j], multiplier[j], abundances[:, j].sum(), args.normalize])
                block = temp / 'species.csv'
                with block.open('w', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(NAMES + batch_samples)
                    for key, vals in zip(keys, abundances):
                        w.writerow(list(key) + vals.tolist())
                destination = out / f'batch_{start // args.batch_size + 1:04d}'
                destination.mkdir()
                subprocess.run([rscript, str(worker), str(block), str(destination)], check=True)
                print(f'Completed samples {start + 1}-{start + len(batch_samples)} / {len(samples)}', file=sys.stderr)
            manifest['status'] = 'complete'
            (out / 'run.json').write_text(json.dumps(manifest, indent=2))
        finally:
            checks_file.close(); totals_file.close(); input_rank_file.close()
            del np_values
    print(f'Completed: {out.resolve()}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', help='Merged MetaPhlAn-style CSV; percentages or other nonnegative abundances')
    p.add_argument('--out', required=True, help='New or empty output directory')
    p.add_argument('--batch-size', type=int, default=32, help='Samples per phyloseq batch (default 32)')
    p.add_argument('--normalize', choices=['observed', 'none'], default='observed',
                   help='observed: scale each sample to 100 once; none: preserve input denominator')
    p.add_argument('--atol', type=float, default=0.001, help='Absolute parent discrepancy tolerance, input units')
    p.add_argument('--rtol', type=float, default=1e-8)
    p.add_argument('--taxid-map', help='Optional CSV: clade_name,taxid; exact input paths with verified NCBI IDs')
    p.add_argument('--taxonkit-data', help='Directory containing NCBI taxonomy dump for TaxonKit')
    p.add_argument('--temp-dir', help='Disk location for temporary memory-mapped abundance matrix')
    args = p.parse_args()
    try:
        run(args)
    except (ValueError, OSError, subprocess.CalledProcessError, StopIteration) as e:
        print(f'ERROR: {e}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
