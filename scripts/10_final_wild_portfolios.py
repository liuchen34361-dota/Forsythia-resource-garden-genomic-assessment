#!/usr/bin/env python3
"""Calculate final portfolios directly from the biallelic diversity VCF GT.

Wild-polymorphic minor alleles define the universe (ties choose ALT).
Called alleles in partial genotypes count. Selection itself uses complete
diploid GT; this independent scan measures the resulting portfolios.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from analysis_utils import (default_config_path, load_config, input_path,
    read_sample_list, read_tsv, open_text, prepare_step_output, step_output_dir,
    write_tsv)


def called_gt(value, n_alt=1):
    tokens = value.replace('|', '/').split('/')
    alleles = [int(t) for t in tokens if t != '.']
    if any(a < 0 or a > n_alt for a in alleles):
        raise ValueError(f'GT outside allele range: {value}')
    return alleles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(default_config_path()))
    cfg = load_config(parser.parse_args().config)
    out = prepare_step_output(cfg, '10_final_wild_portfolios')
    wild = read_sample_list(input_path(cfg, 'wild_samples'))
    current = read_sample_list(input_path(cfg, 'exsitu_samples'))
    core = [r['SampleID'] for r in read_tsv(step_output_dir(cfg, '05_PIHAT_components_and_core') / '06.deredundant_exsitu_core.tsv')]
    selection = step_output_dir(cfg, '07_greedy_plan_selection')
    plan_a = [r['SampleID'] for r in read_tsv(selection / '01.Greedy_Plan_A_Wild_supplements.tsv')]
    plan_b = [r['SampleID'] for r in read_tsv(selection / '02.Greedy_Plan_B_Wild_supplements.tsv')]
    if set(wild) & set(current) or not set(core) <= set(current):
        raise ValueError('Invalid Wild/Current/Core partition')
    if not set(plan_a) <= set(plan_b) <= set(wild):
        raise ValueError('Plans must be nested subsets of Wild')
    groups = {'Current': current, 'Core': core, 'Plan_A': core + plan_a, 'Plan_B': core + plan_b}
    if any(len(v) != len(set(v)) for v in groups.values()):
        raise ValueError('Duplicate portfolio sample')
    counts = Counter()
    header = None
    previous = None
    closed_chroms = set()
    with open_text(input_path(cfg, 'diversity_snp_vcf')) as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')[9:]
                if len(header) != len(set(header)):
                    raise ValueError('Duplicate VCF sample')
                if not set(wild + current) <= set(header):
                    raise ValueError('VCF lacks requested samples')
                index = {s: i for i, s in enumerate(header)}
                wi = [index[s] for s in wild]
                gi = {k: [index[s] for s in v] for k, v in groups.items()}
                continue
            if line.startswith('#'):
                continue
            if header is None:
                raise ValueError('No VCF sample header')
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 + len(header):
                raise ValueError('Malformed VCF row')
            chrom, pos = fields[0], int(fields[1])
            if previous:
                if chrom == previous[0] and pos <= previous[1]:
                    raise ValueError('Duplicate or unsorted VCF coordinate')
                if chrom != previous[0]:
                    closed_chroms.add(previous[0])
                    if chrom in closed_chroms:
                        raise ValueError('Noncontiguous chromosome')
            previous = chrom, pos
            ref, alt = fields[3].upper(), fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
                raise ValueError(f'Expected a biallelic SNP at {chrom}:{pos}')
            counts['records'] += 1
            fmt = fields[8].split(':')
            gt_idx = fmt.index('GT')
            genotypes = []
            for cell in fields[9:]:
                vals = cell.split(':')
                if gt_idx >= len(vals):
                    raise ValueError(f'Missing GT field at {chrom}:{pos}')
                genotypes.append(called_gt(vals[gt_idx]))
            ac = sum(sum(genotypes[i]) for i in wi)
            an = sum(len(genotypes[i]) for i in wi)
            if not 0 < ac < an:
                continue
            target = 1 if ac <= an - ac else 0
            counts['target_total'] += 1
            present = {k: any(target in genotypes[i] for i in indices) for k, indices in gi.items()}
            for k, yes in present.items():
                counts[k + '_missing'] += not yes
            new_loss = present['Current'] and not present['Core']
            counts['core_added_loss'] += new_loss
            for plan in ('Plan_A', 'Plan_B'):
                counts[plan + '_rescued_loss'] += new_loss and present[plan]
                counts[plan + '_recovered_original'] += not present['Current'] and present[plan]
            if counts['records'] % 1000000 == 0:
                print(f"[INFO] Portfolio scan: {counts['records']:,} records", flush=True)
    if header is None or counts['target_total'] == 0:
        raise ValueError('Empty VCF or wild target universe')
    if counts['Core_missing'] != counts['Current_missing'] + counts['core_added_loss']:
        raise AssertionError('Core loss identity failed')
    for plan in ('Plan_A', 'Plan_B'):
        expected = counts['Current_missing'] - counts[plan + '_recovered_original'] + counts['core_added_loss'] - counts[plan + '_rescued_loss']
        if counts[plan + '_missing'] != expected:
            raise AssertionError(f'{plan} loss/recovery identity failed')
    rows = []
    for name, samples in groups.items():
        n, missing = counts['target_total'], counts[name + '_missing']
        rows.append({'Portfolio': name, 'Sample_N': len(samples), 'Target_total': n,
            'Retained': n - missing, 'Missing': missing, 'Retention_rate': (n - missing) / n})
    write_tsv(out / '01.final_portfolios.tsv', rows, list(rows[0]))
    (out / '02.loss_and_rescue.json').write_text(json.dumps(dict(counts), indent=2) + '\n')
    expected = cfg.get('expected_portfolios', {})
    checks = [{'Metric': k, 'Observed': counts[k], 'Expected': v,
        'Status': 'PASS' if counts[k] == v else 'FAIL'} for k, v in expected.items()]
    write_tsv(out / '03.anchor_validation.tsv', checks, ['Metric', 'Observed', 'Expected', 'Status'])
    if any(r['Status'] != 'PASS' for r in checks):
        raise AssertionError('Portfolio anchors differ; see 03.anchor_validation.tsv')
    print(f'[OK] Direct portfolio calculation complete: {out}')


if __name__ == '__main__':
    main()
