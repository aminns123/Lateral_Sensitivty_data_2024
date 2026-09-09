"""New, seeded reversal-subset profiles; never overwrite archived measurements.

Uses the exclusion candidate that reproduces input zero. Compatibility with all
historical inputs is audited separately; see docs/KNOWN_LIMITATIONS.md.
"""
from pathlib import Path
import argparse,csv,gzip,itertools,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SUBSETS=np.array(list(itertools.combinations(range(1,10),5)))

def candidate_pool(cid):
    groups=[g for g in json.loads((ROOT/'config/staircase_groups.json').read_text())
            if g['condition_id']==cid and g['included_for_baseline_candidate']]
    positions=sorted({g['position'] for g in groups if g['experiment']=='base'} &
                     {g['position'] for g in groups if g['experiment']=='flanker'})
    ratios=[];spreads=[]
    for pos in positions:
        thresholds=[]
        for role in ['base','flanker']:
            stairs=[np.array(g['reversals']) for g in groups if g['position']==pos and g['experiment']==role]
            if any(len(r)<10 for r in stairs):raise ValueError(f'{cid}: fewer than ten reversals in an included staircase')
            thresholds.append(np.array([np.exp(np.mean(np.log(r[SUBSETS]),axis=1)) for r in stairs]))
        b,f=thresholds
        ratios.append(np.log((1/f.mean(axis=0))/(1/b.mean(axis=0))))
        pairwise=np.array([np.log((1/ff)/(1/bb)) for bb in b for ff in f])
        spreads.append(np.std(pairwise,axis=0,ddof=1))
    return np.array(positions),np.array(ratios),np.array(spreads)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--condition',required=True);p.add_argument('--seed',required=True,type=int)
    p.add_argument('--samples',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.samples<1:raise ValueError('samples must be positive')
    known={c['condition_id'] for c in json.loads((ROOT/'config/conditions.json').read_text())}
    if args.condition not in known:raise ValueError('Unknown condition')
    out=args.output.resolve()
    # New runs belong outside the immutable reference package.
    if out==ROOT or ROOT in out.parents:raise ValueError('Choose an output directory outside this reference package')
    x,r,s=candidate_pool(args.condition)
    out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(args.seed);seen=set()
    with gzip.open(out/'new_profiles.csv.gz','wt',newline='',encoding='utf-8') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['new_sample_index','point_index','probe_position_deg','log_sensitivity_ratio','spread','subset_index','reversal_indices_zero_based'])
        for sample in range(args.samples):
            while True:
                vector=tuple(int(v) for v in rng.integers(0,len(SUBSETS),size=len(x)))
                if vector not in seen:seen.add(vector);break
            for j,k in enumerate(vector):w.writerow([sample,j,x[j],r[j,k],s[j,k],k,';'.join(map(str,SUBSETS[k]))])
    (out/'run.json').write_text(json.dumps(dict(condition=args.condition,seed=args.seed,samples=args.samples,
        generator='numpy.default_rng',numpy_version=np.__version__,subsets=126,
        scope='New draws with unique joint subset vectors. No deterministic input-zero entry. No historical seed or ordering claimed.',
        exclusion_policy='baseline-matching candidate',fit_results_included=False),indent=2)+'\n')
    print(out)
if __name__=='__main__':main()
