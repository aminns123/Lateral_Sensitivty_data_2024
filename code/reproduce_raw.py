"""Recompute lateral profiles from preserved trials; audit input zero without altering it."""
from pathlib import Path
import collections,csv,gzip,io,itertools,json
import numpy as np
import raw_helpers as helper
ROOT=Path(__file__).resolve().parents[1]
def read_csv(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def matrix(rel):return np.loadtxt(io.BytesIO(gzip.decompress((ROOT/rel).read_bytes())),ndmin=2)
def threshold(reversals,indices=None):
    rev=np.array(reversals);subset=rev[-5:] if indices is None else rev[indices]
    return float(np.exp(np.nanmean(np.log(subset))))
def profile(groups,excluded):
    thresholds=collections.defaultdict(list)
    for g in groups:
        if (g['experiment'],g['file_index'],g['source_staircase_id']) in excluded:continue
        thresholds[(g['experiment'],g['position'])].append(threshold(g['reversals']))
    output=[]
    positions=sorted({pos for role,pos in thresholds if role=='base'} & {pos for role,pos in thresholds if role=='flanker'})
    for pos in positions:
        base=thresholds[('base',pos)];flank=thresholds[('flanker',pos)]
        pairs=[np.log((1/f)/(1/b)) for b in base for f in flank]
        output.append(dict(probe_position_deg=pos,distance_from_flanker_edge_deg=pos+.5,base_threshold=np.mean(base),flanker_threshold=np.mean(flank),base_sensitivity=1/np.mean(base),flanker_sensitivity=1/np.mean(flank),log_sensitivity_ratio=np.log((1/np.mean(flank))/(1/np.mean(base))),spread=np.std(pairs,ddof=1),base_staircases=len(base),flanker_staircases=len(flank)))
    return output
def main():
    conditions=json.loads((ROOT/'config/conditions.json').read_text());filelist=read_csv(ROOT/'data/raw/trial_files.csv');recorded=read_csv(ROOT/'data/raw/recorded_exclusions.csv')
    allgroups=[];stairs=[];reversals=[];checks=[];profiles=[];candidate_checks=[]
    for c in conditions:
        cid=c['condition_id'];groups=[]
        for file in [r for r in filelist if r['condition_id']==cid]:
            a=matrix(file['public_file'])
            positions=list(dict.fromkeys(abs(v) for v in a[:,0]));ids=list(dict.fromkeys(a[:,4]))
            ordinal=0
            for pos in positions:
                for originalid in ids:
                    selected=np.flatnonzero((abs(a[:,0])==pos)&(a[:,4].astype(int)==int(originalid)))
                    if len(selected)==0:continue
                    ordinal+=1;sid=file['source_id']+f'_A{ordinal:02d}'
                    used=selected[:200];contrast=helper.weberContrast(a[used,6].tolist(),c['background_digital_intensity'],1.0,np.log10,False)
                    count,indices,rev=helper.count_reversals_HighLow(contrast)
                    g=dict(condition_id=cid,experiment=file['experiment'],file_index=int(file['source_directory_index']),source_id=file['source_id'],source_staircase_id=int(originalid),staircase_id=sid,position=float(a[used[0],0]),reversals=[float(v) for v in rev],source_rows=[int(v)+1 for v in used])
                    groups.append(g)
                    stairs.append(dict(condition_id=cid,experiment=file['experiment'],staircase_id=sid,source_id=file['source_id'],source_staircase_id=originalid,probe_position_deg=g['position'],trials_recorded=len(selected),trials_used=len(used),distinct_signed_positions=len(set(a[used,0])),detected_reversals=count,final_five_geometric_threshold=threshold(rev)))
                    for j,(idx,v) in enumerate(zip(indices,rev)):
                        reversals.append(dict(condition_id=cid,staircase_id=sid,reversal_index=j,trial_within_staircase=int(idx),source_row=int(used[int(idx)])+1,normalized_contrast=v,used_in_final_five=j>=len(rev)-5))
        candidates={'no_exclusions':set()}
        recorded_set={(r['experiment'],int(float(r['recorded_file_index'])),int(float(r['recorded_staircase_id']))) for r in recorded if r['condition_id']==cid and float(r['recorded_file_index'])!=222}
        candidates['recorded_exclusion_list']=recorded_set
        if c['participant_id']=='P01':
            candidates['storage_script_active_exclusions']={('base',3,3),('base',2,2),('flanker',6,2),('flanker',21,8),('flanker',16,5),('flanker',3,0),('flanker',4,0)}
        saved=matrix(f'data/source/baseline_profiles/{cid.lower()}_profile.tsv.gz')
        spread=matrix(f'data/source/baseline_profiles/{cid.lower()}_spread.tsv.gz')
        results=[]
        for name,excluded in candidates.items():
            out=profile(groups,excluded);a=np.array([[r['probe_position_deg'],r['log_sensitivity_ratio'],r['spread']] for r in out]);same_x=a.shape[0]==saved.shape[0] and np.array_equal(a[:,0],saved[:,0])
            err=float(np.max(abs(a[:,1]-saved[:,1]))) if same_x else None
            serr=float(np.max(abs(a[:,2]-spread[:,1]))) if same_x else None
            check=dict(condition_id=cid,exclusion_candidate=name,points=len(out),same_positions=same_x,max_ratio_difference=err,max_spread_difference=serr,ratio_match_at_1e_12=same_x and err<1e-12,spread_match_at_1e_12=same_x and serr<1e-12)
            candidate_checks.append(check);results.append((check,out,excluded))
        best=min(results,key=lambda r:float('inf') if r[0]['max_ratio_difference'] is None else r[0]['max_ratio_difference'])
        check,out,excluded=best
        checks.append(check)
        for r in out:profiles.append(dict(condition_id=cid,exclusion_candidate=check['exclusion_candidate'],**r))
        for g in groups:g['included_for_baseline_candidate']=(g['experiment'],g['file_index'],g['source_staircase_id']) not in excluded
        allgroups.extend(groups)
        print(cid,check,flush=True)
    (ROOT/'validation/recomputed').mkdir(exist_ok=True)
    write_csv(ROOT/'data/processed/staircase_summary.csv',stairs)
    write_csv(ROOT/'data/processed/staircase_reversals.csv',reversals)
    write_csv(ROOT/'data/processed/lateral_profiles_recomputed.csv',profiles)
    write_csv(ROOT/'validation/raw_input_zero_checks.csv',checks)
    write_csv(ROOT/'validation/exclusion_candidate_checks.csv',candidate_checks)
    (ROOT/'config/staircase_groups.json').write_text(json.dumps(allgroups,separators=(',',':'))+'\n',encoding='utf-8')
    result={'status':'completed','conditions':len(conditions),'staircases':len(stairs),'detected_reversals':len(reversals),'profile_points':len(profiles),'conditions_ratio_match_at_1e_12':sum(r['ratio_match_at_1e_12'] for r in checks),'conditions_spread_match_at_1e_12':sum(r['spread_match_at_1e_12'] for r in checks)}
    (ROOT/'validation/raw_reproduction_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
