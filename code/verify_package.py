"""Verify preserved bytes, table conversions, counts, source helpers and all bootstrap inputs.

Default verification is read-only. --report writes a JSON report at a requested
path; --skip-manifest is for maintainers before generating a new manifest.
"""
from pathlib import Path
import argparse,ast,csv,gzip,hashlib,itertools,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
FIT_COLUMNS=['fms_ra','fms_lof','fms_ransac','phase_consistency_ra','phase_consistency_lof','phase_consistency_ransac','tolerance_fraction_ra','tolerance_fraction_lof','tolerance_fraction_ransac','f_n_ra_cpd','f_n_lof_cpd','f_n_ransac_cpd','f_n_main_cpd','fms_main','fisher_information_recorded','cr_bound_recorded','f_n_aux_recorded_cpd']
EXPECTED_CLUSTERS={'P01_L010':[11748,2254],'P01_L020':[22042,875],'P01_L025':[3068,20473],'P01_L050':[22027,33],'P01_L200':[1,15346],'P02_L010':[14930,14871],'P02_L026':[8119,3943],'P02_L049':[1888,1312],'P02_L100':[2177,3281],'P04_L010':[8536,3739],'P04_L049':[32944,811],'P04_L100':[330,24882],'P04_L200':[20605,4289]}

def read(path):
    with (gzip.open(path,'rt',encoding='utf-8',newline='') if path.suffix=='.gz' else path.open(encoding='utf-8',newline='')) as f:return list(csv.DictReader(f))
def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
def files():return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(v in {'.git','.venv','__pycache__'} for v in p.relative_to(ROOT).parts) and p.name!='manifest-sha256.csv')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--skip-manifest',action='store_true');p.add_argument('--report',type=Path);args=p.parse_args()
    manifest_count=0
    if not args.skip_manifest:
        rows=read(ROOT/'manifest-sha256.csv');actual={p.relative_to(ROOT).as_posix():p for p in files()}
        require(set(actual)=={r['path'] for r in rows},'Manifest file set differs')
        for r in rows:
            b=actual[r['path']].read_bytes();require(len(b)==int(r['bytes']) and sha(b)==r['sha256'],'Manifest mismatch: '+r['path'])
        manifest_count=len(rows)
    sources=read(ROOT/'source_manifest.csv')
    for r in sources:
        path=ROOT/r['public_file'];b=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
        require(len(b)==int(r['source_bytes']) and sha(b)==r['source_sha256'],'Source bytes: '+r['source_id'])
    trials=read(ROOT/'data/raw/trial_files.csv');total_trials=0
    for r in trials:
        with gzip.open(ROOT/r['public_file'],'rt') as f:a=np.loadtxt(f,ndmin=2)
        require(a.shape==(int(r['rows']),int(r['columns'])) and np.isfinite(a).all(),'Trial matrix '+r['source_id']);total_trials+=len(a)
    require(len(trials)==560 and total_trials==115445,'Raw inventory count')
    for file in (ROOT/'code').glob('*.py'):ast.parse(file.read_text(encoding='utf8'))
    for r in json.loads((ROOT/'config/function_provenance.json').read_text()):
        src=(ROOT/r['public_file']).read_text(encoding='utf8');tree=ast.parse(src);node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==r['function'])
        text='\n'.join(src.splitlines()[node.lineno-1:node.end_lineno]);require(sha(text.encode())==r['function_sha256'],'Author function changed: '+r['function'])
    conditions=json.loads((ROOT/'config/conditions.json').read_text());require(len(conditions)==17 and {c['participant_id'] for c in conditions}=={'P01','P02','P03','P04'},'Subject coverage')
    fitted=accepted=bootstrap_sources=bootstrap_points=0
    baseline=read(ROOT/'data/processed/lateral_profiles_recomputed.csv')
    modes=read(ROOT/'data/processed/isf_mode_summaries.csv')
    for c in conditions:
        cid=c['condition_id'];stem=cid.lower()
        with gzip.open(ROOT/f'data/source/saved_fits/{stem}.tsv.gz','rt') as f:original=[line.split() for line in f]
        converted=read(ROOT/f'data/processed/{stem}_saved_fits.csv.gz');require(len(original)==len(converted)==50000,cid+' fit rows')
        acc=0
        for index,(a,b) in enumerate(zip(original,converted)):
            require(len(a)==17 and a==[b[k] for k in FIT_COLUMNS],cid+' fit tokens')
            keep=.8<=float(a[0])<=1.;acc+=keep
            require(b['condition_id']==cid and int(b['saved_fit_row'])==index+1 and (b['accepted_fms_ra']=='True')==keep and b['fit_success_recorded']=='',cid+' fit flags')
        require(acc==c['thesis_accepted_count'],cid+' thesis acceptance count');fitted+=50000;accepted+=acc
        if cid in EXPECTED_CLUSTERS:require([int(r['cluster_count']) for r in modes if r['condition_id']==cid and r['variant']=='thesis_percentiles']==EXPECTED_CLUSTERS[cid],cid+' thesis clusters')
        br=[r for r in baseline if r['condition_id']==cid]
        for kind,col in [('profile','log_sensitivity_ratio'),('spread','spread')]:
            with gzip.open(ROOT/f'data/source/baseline_profiles/{stem}_{kind}.tsv.gz','rt') as f:a=np.loadtxt(f)
            b=np.array([[float(r['probe_position_deg']),float(r[col])] for r in br]);require(np.array_equal(a,b),cid+' input zero '+kind)
        with gzip.open(ROOT/f'data/bootstrap/{stem}_profiles.csv.gz','rt',newline='') as data, gzip.open(ROOT/f'validation/{stem}_bootstrap_sources.csv.gz','rt',newline='') as meta:
            source_iter=iter(csv.DictReader(meta));num_inputs=0
            for index,rows in itertools.groupby(csv.DictReader(data),key=lambda r:int(r['source_input_index'])):
                rows=list(rows);require(index==num_inputs,cid+' noncontiguous input index');num_inputs+=1
                require([int(r['source_point_index']) for r in rows]==list(range(len(rows))),cid+' point indices')
                require(all(np.isfinite([float(r[k]) for k in ['probe_position_deg','log_sensitivity_ratio','spread_recorded']]).all() for r in rows),cid+' nonfinite bootstrap data')
                for kind,col in [('profile','log_sensitivity_ratio'),('spread','spread_recorded')]:
                    m=next(source_iter);require(int(m['source_input_index'])==index and m['source_kind']==kind,cid+' bootstrap source order')
                    newline='\r\n' if m['newline']=='CRLF' else '\n'
                    b=(newline.join(r['probe_position_deg']+'\t'+r[col] for r in rows)+(newline if m['terminal_newline']=='True' else '')).encode()
                    require(len(rows)==int(m['source_rows']) and len(b)==int(m['source_bytes']) and sha(b)==m['source_sha256'],cid+f' bootstrap original {index} {kind}');bootstrap_sources+=1
                bootstrap_points+=len(rows)
            require(num_inputs==50000 and next(source_iter,None) is None,cid+' bootstrap input count')
        print(cid,'verified',flush=True)
    # Independently check the new primary measurement exports against their sources.
    groups={g['staircase_id']:g for g in json.loads((ROOT/'config/staircase_groups.json').read_text())}
    conditions_by_id={c['condition_id']:c for c in conditions}
    thresholds=read(ROOT/'data/processed/staircase_thresholds.csv')
    require(len(thresholds)==2469,'Threshold coverage')
    for r in thresholds:
        g=groups[r['staircase_id']]
        value=float(np.exp(np.mean(np.log(g['reversals'][-5:]))))
        require(abs(value-float(r['threshold_last_five']))<1e-14 and (r['included_for_baseline_candidate']=='True')==g['included_for_baseline_candidate'],'Threshold/inclusion export')
    trial_sources={r['source_id']:r for r in trials};cached={};trace_rows=0;detected=0;final_five=0
    recorded_reversals={(r['staircase_id'],int(r['source_row'])):r for r in read(ROOT/'data/processed/staircase_reversals.csv')}
    with gzip.open(ROOT/'data/processed/staircase_trials.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            sid=r['source_id']
            if sid not in cached:
                with gzip.open(ROOT/trial_sources[sid]['public_file'],'rt') as h:cached[sid]=np.loadtxt(h,ndmin=2)
            intensity=cached[sid][int(r['source_row'])-1,6];bkg=conditions_by_id[r['condition_id']]['background_digital_intensity']
            require(intensity==float(r['probe_digital_intensity']) and abs((intensity-bkg)/(1-bkg)-float(r['normalized_contrast']))<1e-14,'Trial export contrast')
            rv=recorded_reversals.get((r['staircase_id'],int(r['source_row'])))
            require((r['is_detected_reversal']=='True')==(rv is not None),'Trial reversal flag')
            if rv:require(r['reversal_index']==rv['reversal_index'] and r['used_in_final_five']==rv['used_in_final_five'],'Trial reversal membership')
            else:require(r['reversal_index']=='' and r['used_in_final_five']=='False','Nonreversal membership')
            g=groups[r['staircase_id']]
            require(g['source_rows'][int(r['trial_within_staircase'])]==int(r['source_row']) and (r['included_for_baseline_candidate']=='True')==g['included_for_baseline_candidate'],'Trial identity')
            trace_rows+=1;detected+=r['is_detected_reversal']=='True';final_five+=r['used_in_final_five']=='True'
    require(trace_rows==115445 and detected==25133 and final_five==2469*5,'Trace counts')
    endpoints=read(ROOT/'data/processed/isf_endpoints.csv');endpoint_count=0
    require({r['condition_id'] for r in endpoints}==set(conditions_by_id),'Endpoint condition coverage')
    for c in conditions:
        cid=c['condition_id']
        with gzip.open(ROOT/f'data/source/saved_fits/{cid.lower()}.tsv.gz','rt') as f:a=np.loadtxt(f)
        d=a[(a[:,0]>=.8)&(a[:,0]<=1),9]
        mm=[r for r in modes if r['condition_id']==cid and r['variant']=='thesis_percentiles' and r['retained_for_display']=='True']
        ee=[r for r in endpoints if r['condition_id']==cid]
        require(len(ee)==max(1,len(mm)),cid+' endpoint count')
        if not mm:require(ee[0]['isf_estimate_cpd']=='' and ee[0]['candidate_mode']=='',cid+' empty endpoint')
        for m in mm:
            e=next(r for r in ee if r['candidate_mode']==m['candidate_mode'])
            subset=d[(d>=float(m['interval_low_cpd']))&(d<=float(m['interval_high_cpd']))]
            low,mid,high=np.percentile(subset,[.5,50,99.5])
            require(np.allclose([float(e['isf_lower_cpd']),float(e['isf_estimate_cpd']),float(e['isf_upper_cpd'])],[low,mid,high],rtol=0,atol=1e-12),'Endpoint quantiles')
            require(abs(float(e['error_minus_cpd'])-(mid-low))<1e-12 and abs(float(e['error_plus_cpd'])-(high-mid))<1e-12,'Endpoint error bars');endpoint_count+=1
    require(endpoint_count==27,'Retained endpoint count')
    # Every distributed CSV header must have field documentation.
    dictionary=read(ROOT/'data_dictionary.csv');documented={(r['file'],r['field']) for r in dictionary}
    for path in files():
        if path.name.endswith(('.csv','.csv.gz')) and path.name not in {'data_dictionary.csv','manifest-sha256.csv'}:
            with (gzip.open(path,'rt',encoding='utf8') if path.suffix=='.gz' else path.open(encoding='utf8')) as f:header=next(csv.reader(f))
            require(all((path.relative_to(ROOT).as_posix(),h) in documented for h in header),'Dictionary coverage: '+path.name)
    result=dict(status='PASS',manifest_files_verified=manifest_count,preserved_source_files_verified=len(sources),raw_files=len(trials),raw_trials=total_trials,
        saved_fit_rows=fitted,accepted_fit_rows=accepted,bootstrap_original_files_reconstructed_and_verified=bootstrap_sources,bootstrap_points=bootstrap_points,
        baseline_points_exactly_reproduced=len(baseline),thesis_condition_counts_matched=17,thesis_cluster_pairs_matched=13,
        exported_trial_rows_verified=trace_rows,staircase_thresholds_verified=len(thresholds),isf_endpoint_quantiles_verified=endpoint_count,
        limitations='Integrity and table checks do not establish historical fitted-row-to-input identity. See docs/KNOWN_LIMITATIONS.md.')
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
