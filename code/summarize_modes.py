"""Summarize preserved RA-screened distributions with explicit source/thesis settings."""
from pathlib import Path
import csv,gzip,io,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.stats import gaussian_kde
from scipy.integrate import simpson
import fit_helpers as helper
from reproduce_raw import write_csv,matrix
ROOT=Path(__file__).resolve().parents[1]

def intervals(data):
    for n in [3,2]:
        try:
            peaks,valleys=helper.get_peaks_vally_distributions(data,n,bandwidthScale=1,nBins=31,kindInterp='linear',lam=.02)
            peaks=sorted(float(p) for p in peaks)
            groups=sorted(helper.find_peak_vally_trio(peaks,valleys,data),key=lambda a:a[1])
            if len(peaks) and all(len(a)==3 for a in groups):return groups,'source_histogram_spline'
        except (IndexError,ValueError):pass
    return [[float(min(data)),float(np.median(data)),float(max(data))]],'median_fallback'

def legacy_choice(mode_ids,accepted_ids,accepted_scores):
    ids=[];scores=[]
    for i in mode_ids:
        i=int(i)
        try:score=accepted_scores[accepted_ids[i]]
        except (IndexError,TypeError):
            try:score=accepted_scores[i]
            except (IndexError,TypeError):continue
        ids.append(i);scores.append(score)
    return ids[int(np.argmax(scores))] if scores else None

def main():
    rows=[];densities=[];assignments=[];requests=[];comparisons=[]
    for c in json.loads((ROOT/'config/conditions.json').read_text()):
        cid=c['condition_id'];a=matrix(f'data/source/saved_fits/{cid.lower()}.tsv.gz')
        accepted=np.flatnonzero((a[:,0]>=.8)&(a[:,0]<=1));data=a[accepted,9];scores=a[accepted,0]
        groups,algorithm=intervals(data)
        x=np.linspace(float(data.min()),float(data.max()),1000);density=gaussian_kde(data,bw_method='scott')(x)
        for xx,yy in zip(x,density):densities.append(dict(condition_id=cid,frequency_cpd=xx,density_per_cpd=yy))
        for variant,lo,hi in [('thesis_percentiles',.5,99.5),('current_driver_percentiles',2.5,97.5)]:
            candidates=[]
            for n,(low,peak,high) in enumerate(groups,1):
                mask=(data>=low)&(data<=high);cluster=data[mask];ids=accepted[mask]
                if len(cluster)==0:continue
                qlo,median,qhi=np.percentile(cluster,[lo,50,hi]);region=(x>=qlo)&(x<=qhi)
                mass=float(simpson(density[region],x=x[region])) if region.sum()>=2 else 0.
                correct=int(ids[np.argmax(a[ids,0])]);legacy=legacy_choice(ids,accepted,scores)
                candidates.append(dict(condition_id=cid,variant=variant,candidate_mode=n,peak_locator=algorithm,interval_low_cpd=low,spline_peak_cpd=peak,interval_high_cpd=high,cluster_count=len(cluster),cluster_fraction_of_accepted=len(cluster)/len(data),percentile_low=lo,percentile_high=hi,quantile_low_cpd=qlo,median_cpd=median,quantile_high_cpd=qhi,kde_mass=mass,retained_for_display=mass>=.1,representative_input_index_row_linked=correct,representative_input_index_source_code=legacy,representative_selection_agrees=correct==legacy,representative_saved_fms=a[correct,0],included_in_thesis_main_isf=c['included_in_thesis_main_isf']))
            qualifying=sorted([r for r in candidates if r['retained_for_display']],key=lambda r:r['kde_mass'],reverse=True)[:2]
            for r in candidates:
                r['retained_for_display']=r in qualifying
                r['dominant_retained_mode']=bool(qualifying) and r is qualifying[0]
                rows.append(r)
                if variant=='thesis_percentiles' and r['retained_for_display']:
                    for selection,key in [('row_linked_max_fms','representative_input_index_row_linked'),('source_plotting_index_rule','representative_input_index_source_code')]:
                        if r[key] is not None:requests.append({'condition_id':cid,'candidate_mode':r['candidate_mode'],'selection':selection,'source_input_index':int(r[key]),'spline_peak_cpd':float(r['spline_peak_cpd'])})
            if variant=='thesis_percentiles':
                for index in accepted:
                    membership=[r['candidate_mode'] for r in candidates if r['interval_low_cpd']<=a[index,9]<=r['interval_high_cpd']]
                    assignments.append([cid,int(index)+1,int(index),float(a[index,9]),';'.join(map(str,membership))])
        print(cid,'accepted',len(accepted),'clusters',[r['cluster_count'] for r in rows if r['condition_id']==cid and r['variant']=='thesis_percentiles'],flush=True)
    write_csv(ROOT/'data/processed/isf_mode_summaries.csv',rows)
    write_csv(ROOT/'data/processed/isf_kde.csv',densities)
    with gzip.open(ROOT/'data/processed/mode_assignments.csv.gz','wt',encoding='utf-8',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['condition_id','saved_fit_row','assumed_input_index','f_n_ra_cpd','candidate_modes']);w.writerows(assignments)
    (ROOT/'config/representative_requests.json').write_text(json.dumps(requests,indent=2)+'\n')
    summary={'conditions':17,'accepted_rows':len(assignments),'candidate_clusters_thesis_percentiles':sum(r['variant']=='thesis_percentiles' for r in rows),'representative_rule_disagreements':sum(r['variant']=='thesis_percentiles' and r['retained_for_display'] and not r['representative_selection_agrees'] for r in rows),'source_peak_locator':'smoothed histogram spline; KDE used separately for display and mass','original_figure_pixel_identity_claimed':False,'saved_row_to_input_pairing_verified':False}
    (ROOT/'validation/mode_summary_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
