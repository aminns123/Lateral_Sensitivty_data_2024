"""Check every saved nonzero input against the 126 paired reversal-subset candidates."""
import csv,gzip,json
import numpy as np
from resample_profiles import ROOT,candidate_pool
from reproduce_raw import write_csv

def main():
    results=[]
    for c in json.loads((ROOT/'config/conditions.json').read_text()):
        cid=c['condition_id'];x,ratios,spreads=candidate_pool(cid)
        with gzip.open(ROOT/f'data/bootstrap/{cid.lower()}_profiles.csv.gz','rt') as f:
            a=np.loadtxt(f,delimiter=',',skiprows=1)
        a=a.reshape(50000,len(x),5)
        if not np.array_equal(a[:,:,2],np.broadcast_to(x,(50000,len(x)))):raise ValueError(cid+' positions differ')
        for j,pos in enumerate(x):
            # Joint ratio AND spread must come from the SAME subset; duplicates allowed.
            observed=a[1:,j,3:5];matched=np.zeros(len(observed),bool);closest=np.full(len(observed),np.inf)
            for rr,ss in zip(ratios[j],spreads[j]):
                delta=np.max(np.abs(observed-[rr,ss]),axis=1)
                closest=np.minimum(closest,delta);matched|=delta<=1e-12
            results.append(dict(condition_id=cid,probe_position_deg=pos,checked_inputs=49999,
                compatible_inputs=int(matched.sum()),incompatible_inputs=int((~matched).sum()),
                max_nearest_joint_difference=float(closest.max()),tolerance=1e-12))
        print(cid,'compatible',sum(r['compatible_inputs'] for r in results if r['condition_id']==cid),flush=True)
    write_csv(ROOT/'validation/resampling_subset_checks.csv',results)
    summary=dict(checked_points=sum(r['checked_inputs'] for r in results),compatible_points=sum(r['compatible_inputs'] for r in results),
        incompatible_points=sum(r['incompatible_inputs'] for r in results),scope='Nonzero input indices only; baseline-matching exclusions; common reversal subset at each position.',
        historical_seed_recovered=False,historical_fit_row_link_verified=False)
    (ROOT/'validation/resampling_subset_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
if __name__=='__main__':main()
