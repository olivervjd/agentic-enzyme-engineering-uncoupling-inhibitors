"""Fresh protein-only G177 refolding and matched-WT structural comparison."""
import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import tarfile
import time

SEEDS = (211, 223, 227)
ORDER = 'GACDEFHIKLMNPQRSTVWY'
NAMES = dict(zip('ACDEFGHIKLMNPQRSTVWY', ['Alanine','Cysteine','Aspartate','Glutamate',
    'Phenylalanine','Glycine','Histidine','Isoleucine','Lysine','Leucine','Methionine',
    'Asparagine','Proline','Glutamine','Arginine','Serine','Threonine','Valine','Tryptophan','Tyrosine']))

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')
    tmp.replace(p)

def alignment_counts(raw):
    """Foldseek alnlen includes gaps; coverage is not the paired-residue count."""
    q, t = raw['qaln'], raw['taln']
    columns = int(raw['alnlen'])
    if len(q) != len(t) or len(q) != columns:
        raise ValueError('Alignment string lengths differ from alnlen')
    paired = sum(a != '-' and b != '-' for a, b in zip(q, t))
    qcount, tcount = sum(a != '-' for a in q), sum(a != '-' for a in t)
    qlen, tlen = int(raw['qlen']), int(raw['tlen'])
    if not (0 < paired <= qcount <= qlen and paired <= tcount <= tlen):
        raise ValueError('Invalid paired-residue counts')
    if qcount + tcount - paired != columns:
        raise ValueError('Alignment contains a double-gap column')
    if abs(qcount / qlen - float(raw['qcov'])) >= .001 or abs(tcount / tlen - float(raw['tcov'])) >= .001:
        raise ValueError('Alignment coverage differs from residue counts')
    return dict(alignment_columns=columns, paired_residues=paired, paired_coverage=paired/tlen)

def prepare(root, donor):
    if root.exists():
        raise FileExistsError(root)
    old = json.loads((donor / 'manifest.json').read_text())
    full = old['full_sequence']
    assert len(full) == 520 and full[176] == 'G'
    assert hashlib.sha256(full.encode()).hexdigest() == 'ac2d95791932eca90a6986cbe18b4d9ab948b08d42fb68b871b9b777d8b51bcf'
    manifest = dict(created_at=now(), target=old['target'], full_sequence=full,
        canonical_position=177, mature_position=101, mature_range=[77,520],
        donor_manifest_sha256=sha(donor/'manifest.json'), seeds=list(SEEDS), records=[],
        state='apo / protein only; no ligand or template',
        protocol=dict(boltz='2.2.1', recycling_steps=3, sampling_steps=100, diffusion_samples=1),
        msa_policy='Reuse identical homolog rows from prior ColabFold MSA, changing only the mutant query',
        reference_policy='Predicted WT of the same mature construct, method, ligand-free state and seed; not an experimental structure')
    common = None
    for aa in ORDER:
        name = 'WT' if aa == 'G' else f'G177{aa}'
        sequence = full[76:176] + aa + full[177:]
        assert len(sequence) == 444 and sequence[100] == aa
        assert sum(a != b for a,b in zip(sequence, full[76:])) == (aa != 'G')
        source = next(r for r in old['records'] if r['variant'] == name)
        msa = donor / source['msa']
        assert sha(msa) == source['msa_sha256']
        lines = msa.read_text().splitlines(keepends=True)
        headers = [i for i,l in enumerate(lines) if l.startswith('>')]
        assert ''.join(l.strip() for l in lines[headers[0]+1:headers[1]]) == sequence
        tail = ''.join(lines[headers[1]:])
        common = common or tail
        assert tail == common
        dest = root / 'msa' / f'{name}.a3m'
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(msa, dest)
        for seed in SEEDS:
            rid = f'{name}-apo-{seed}'
            request = {'version':1, 'sequences':[{'protein':{'id':'A', 'sequence':sequence, 'msa':str(dest)}}]}
            path = root / 'inputs' / f'seed-{seed}' / f'{rid}.yaml'
            save(path, request)
            manifest['records'].append(dict(id=rid, variant=name, amino_acid=aa, seed=seed,
                sequence=sequence, request=str(path.relative_to(root)), request_sha256=sha(path),
                msa=str(dest.relative_to(root)), msa_sha256=sha(dest)))
    assert len(manifest['records']) == 60
    manifest['homolog_rows_sha256'] = hashlib.sha256(common.encode()).hexdigest()
    save(root/'manifest.json', manifest)
    save(root/'status.json', dict(state='PREPARED', completed=0, expected=60, updated_at=now()))
    print('Prepared 60 fresh apo predictions', flush=True)

def output_paths(root, row):
    folder = root/'boltz'/f"seed-{row['seed']}"/f"boltz_results_seed-{row['seed']}"/'predictions'/row['id']
    return folder/f"{row['id']}_model_0.cif", folder/f"confidence_{row['id']}_model_0.json"

def run(root, cache, executable):
    m = json.loads((root/'manifest.json').read_text())
    assert importlib.metadata.version('boltz') == m['protocol']['boltz']
    if (root/'execution.json').exists():
        raise FileExistsError('Refusing duplicate run')
    for r in m['records']:
        assert sha(root/r['request']) == r['request_sha256']
        assert sha(root/r['msa']) == r['msa_sha256']
    events = []
    save(root/'execution.json', events)
    for seed in SEEDS:
        out = root/'boltz'/f'seed-{seed}'
        assert not out.exists()
        cmd = [str(executable),'predict',str(root/'inputs'/f'seed-{seed}'), '--out_dir',str(out),
            '--cache',str(cache),'--model','boltz2','--accelerator','gpu','--devices','1',
            '--seed',str(seed),'--recycling_steps','3','--sampling_steps','100',
            '--diffusion_samples','1','--num_workers','0','--no_kernels','--output_format','mmcif']
        event = dict(seed=seed, started_at=now(), command=cmd, state='RUNNING')
        events.append(event)
        save(root/'execution.json', events)
        save(root/'status.json', dict(state='RUNNING', seed=seed, updated_at=now()))
        start = time.monotonic()
        with (root/f'seed-{seed}.log').open('w') as log:
            p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        event.update(returncode=p.returncode, finished_at=now(), elapsed_seconds=time.monotonic()-start)
        count = sum(all(p.is_file() for p in output_paths(root,r)) for r in m['records'])
        event['state'] = 'COMPLETED' if p.returncode == 0 else 'FAILED'
        save(root/'execution.json', events)
        save(root/'status.json', dict(state='RUNNING' if p.returncode == 0 else 'FAILED', seed=seed,
            completed=count, expected=60, updated_at=now()))
        print(f'Seed {seed}: {count}/60 structures, {event["elapsed_seconds"]:.1f}s', flush=True)
        if p.returncode:
            raise RuntimeError('Boltz failed; inspect log')
    assert count == 60
    save(root/'status.json', dict(state='REFOLDING_COMPLETE', completed=count, expected=60, updated_at=now()))

def analyze(root, foldseek):
    from Bio.PDB import MMCIFParser, PDBIO, Superimposer
    from Bio.SeqUtils import seq1
    m = json.loads((root/'manifest.json').read_text())
    structs, artifacts, identities = {}, [], {}
    for r in m['records']:
        cif, confidence = output_paths(root, r)
        assert cif.is_file() and confidence.is_file()
        assert sha(root/r['request']) == r['request_sha256'] and sha(root/r['msa']) == r['msa_sha256']
        model = MMCIFParser(QUIET=True).get_structure(r['id'], str(cif))[0]
        assert list(model.child_dict) == ['A']
        residues = list(model['A'].get_residues())
        assert len(residues) == 444 and all('CA' in x for x in residues)
        assert ''.join(seq1(x.resname) for x in residues) == r['sequence']
        assert seq1(residues[100].resname) == r['amino_acid']
        structs[r['id']] = residues
        conf = json.loads(confidence.read_text())
        assert math.isfinite(conf['complex_plddt'])
        dest = root/'structures'/f"seed-{r['seed']}"/f"{r['id']}.pdb"
        dest.parent.mkdir(parents=True, exist_ok=True)
        io = PDBIO()
        io.set_structure(model)
        io.save(str(dest))
        identities[dest.stem] = r
        for path in (cif,confidence,dest):
            artifacts.append(dict(path=str(path.relative_to(root)), sha256=sha(path)))
    fsroot = root/'foldseek'
    fsroot.mkdir(exist_ok=True)
    version = subprocess.check_output([str(foldseek),'version'], text=True).strip()
    columns = ['query','target','qtmscore','ttmscore','lddt','qcov','tcov','alnlen','qlen','tlen','qaln','taln']
    pairs, executions = [], []
    # Each of the 60 queries is compared to all three WT predictions. Main results
    # use matched seeds; cross-seed WT comparisons quantify prediction variability.
    wt_dir = root/'references'
    wt_dir.mkdir(exist_ok=True)
    for seed in SEEDS:
        src = root/'structures'/f'seed-{seed}'/f'WT-apo-{seed}.pdb'
        shutil.copyfile(src, wt_dir/src.name)
    for seed in SEEDS:
        tsv = fsroot/f'seed-{seed}.tsv'
        cmd = [str(foldseek),'easy-search',str(root/'structures'/f'seed-{seed}'),str(wt_dir),str(tsv),
            str(fsroot/f'tmp-{seed}'),'--alignment-type','1','--exhaustive-search','1',
            '--format-output',','.join(columns),'-a','1','--threads','4','--max-seqs','1000']
        with (fsroot/f'seed-{seed}.log').open('w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        executions.append(cmd)
        for line in tsv.read_text().splitlines():
            vals = line.split('\t')
            assert len(vals) == len(columns)
            raw = dict(zip(columns,vals))
            def identify(value):
                for name in identities:
                    if value in (name, name+'.pdb', name+'_A', name+'.pdb_A'):
                        return identities[name]
                raise ValueError(f'Unexpected Foldseek identifier {value}')
            q, t = identify(raw['query']), identify(raw['target'])
            assert t['variant'] == 'WT'
            metrics = {k:float(raw[k]) for k in ['qtmscore','ttmscore','lddt','qcov','tcov']}
            assert all(math.isfinite(v) and 0 <= v <= 1 for v in metrics.values())
            assert int(raw['qlen']) == int(raw['tlen']) == 444
            counts = alignment_counts(raw)
            sup = Superimposer()
            sup.set_atoms([r['CA'] for r in structs[t['id']]], [r['CA'] for r in structs[q['id']]])
            conf = json.loads(output_paths(root,q)[1].read_text())
            pairs.append(dict(variant=q['variant'], amino_acid=q['amino_acid'], query_id=q['id'],
                target_id=t['id'], seed=q['seed'], reference_seed=t['seed'],
                matched_seed=q['seed']==t['seed'], **metrics, **counts,
                full_chain_ca_rmsd_angstrom=float(sup.rms), complex_plddt=conf['complex_plddt']))
    assert len(pairs) == 180
    assert len({(r['query_id'],r['target_id']) for r in pairs}) == 180
    matched = [r for r in pairs if r['matched_seed']]
    assert len(matched) == 60
    def stats(xs):
        return dict(n=len(xs), mean=statistics.mean(xs), sd=statistics.stdev(xs), min=min(xs), max=max(xs))
    rows = []
    keys = ['ttmscore','qtmscore','lddt','qcov','tcov','alignment_columns','paired_residues','paired_coverage','full_chain_ca_rmsd_angstrom','complex_plddt']
    for aa in ORDER:
        selected = [r for r in matched if r['amino_acid'] == aa]
        assert len(selected) == 3
        rows.append(dict(variant=selected[0]['variant'], amino_acid=aa, name=NAMES[aa],
            **{key:stats([r[key] for r in selected]) for key in keys}))
    baseline = [r for r in pairs if r['variant']=='WT' and r['seed'] < r['reference_seed']]
    assert len(baseline) == 3
    report = dict(complete=True, created_at=now(), predictions=60, foldseek_comparisons=180,
        matched_comparisons=60, foldseek_version=version, foldseek_executable_sha256=sha(foldseek),
        protocol=m, rows=rows, matched_pairs=matched, all_pairs=pairs, wt_cross_seed_pairs=baseline,
        wt_cross_seed_summary={key:stats([r[key] for r in baseline]) for key in keys},
        foldseek_commands=executions, artifacts=artifacts,
        alignment_note='alnlen includes gap columns; paired_residues excludes them. qcov/tcov are alignment span coverage.',
        limitations=['Predicted apo structures; reference is a prediction, not experimental wild type',
            'Boltz2 predictions may be insensitive to subtle single-residue mutation effects',
            'Structural similarity does not demonstrate PEP retention, glyphosate resistance or catalytic activity',
            'Three seeds measure sampling variation; SD is not a calibrated confidence interval',
            'Foldseek TM-align metrics describe its aligned residues; full-chain CA RMSD separately uses all 444 matched positions',
            'WT matched-seed comparisons are self-comparisons; cross-seed WT variation is reported separately'])
    save(root/'results.json', report)
    save(root/'validation.json', dict(status='PASSED', predictions=60, comparisons=180,
        checks=['exact 444-residue mature sequence', 'correct residue 101 / canonical 177',
            'protein-only chain A', 'all CA atoms present', 'input and MSA SHA256',
            'unique 60 predictions and 180 comparisons', 'finite bounded Foldseek scores',
            'raw alignment strings checked; gap columns distinguished from paired residues'], updated_at=now()))
    save(root/'status.json', dict(state='COMPLETED', completed=60, expected=60, updated_at=now()))
    print(json.dumps(dict(complete=True, predictions=60, comparisons=180, foldseek_version=version)), flush=True)

def pack(root):
    report = json.loads((root/'results.json').read_text())
    validation = json.loads((root/'validation.json').read_text())
    if not report['complete'] or report['predictions'] != 60 or validation['status'] != 'PASSED':
        raise ValueError('Refusing to package incomplete or unvalidated results')
    for artifact in report['artifacts']:
        path = (root/artifact['path']).resolve()
        path.relative_to(root.resolve())
        if sha(path) != artifact['sha256']:
            raise ValueError('Artifact hash changed since analysis')
    selected = [root/n for n in ['manifest.json','execution.json','results.json','validation.json','status.json']]
    for pattern in ['structures/**/*.pdb','inputs/**/*.yaml','foldseek/*.tsv','foldseek/*.log',
                    'boltz/**/predictions/**/*.cif','boltz/**/predictions/**/confidence*.json','seed-*.log']:
        selected += sorted(root.glob(pattern))
    archive = root.parent/(root.name+'-artifacts.tar.gz')
    with tarfile.open(archive,'w:gz') as out:
        for p in selected:
            out.add(p, arcname=str(p.relative_to(root)), recursive=False)
    print(json.dumps(dict(archive=str(archive), files=len(selected), sha256=sha(archive))), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['prepare','run','analyze','pack'])
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--donor',type=Path)
    p.add_argument('--cache',type=Path)
    p.add_argument('--boltz',type=Path)
    p.add_argument('--foldseek',type=Path)
    a = p.parse_args()
    a.root = a.root.resolve()
    if a.donor is not None:
        a.donor = a.donor.resolve()
    if a.mode == 'prepare' and a.donor is None:
        p.error('prepare requires --donor')
    if a.mode == 'run' and (a.cache is None or a.boltz is None):
        p.error('run requires --cache and --boltz')
    if a.mode == 'analyze' and a.foldseek is None:
        p.error('analyze requires --foldseek')
    if a.mode=='prepare': prepare(a.root,a.donor)
    elif a.mode=='run': run(a.root,a.cache,a.boltz)
    elif a.mode=='analyze': analyze(a.root,a.foldseek)
    elif a.mode=='pack': pack(a.root)
