"""Export an easy-to-follow trial/reversal/threshold path and descriptive ISF endpoints."""
import csv,gzip,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from reproduce_raw import ROOT,read_csv,write_csv,matrix
import raw_helpers as helper

def main():
    cs=json.loads((ROOT/'config/conditions.json').read_text())
    groups=json.loads((ROOT/'config/staircase_groups.json').read_text())
    sources={r['source_id']:r for r in read_csv(ROOT/'data/raw/trial_files.csv')}
    reversals=read_csv(ROOT/'data/processed/staircase_reversals.csv')
    revmap={(r['staircase_id'],int(r['source_row'])):r for r in reversals}
    summary={r['staircase_id']:r for r in read_csv(ROOT/'data/processed/staircase_summary.csv')}
    modes=read_csv(ROOT/'data/processed/isf_mode_summaries.csv')
    thresholds=[];endpoints=[];condition_rows=[];examples={};cached={};ntrials=0
    fields=['condition_id','experiment','staircase_id','source_id','source_row','trial_within_staircase','probe_position_deg','probe_digital_intensity','normalized_contrast','is_detected_reversal','reversal_index','used_in_final_five','included_for_baseline_candidate']
    with gzip.open(ROOT/'data/processed/staircase_trials.csv.gz','wt',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader()
        for c in cs:
            cid=c['condition_id'];cg=[g for g in groups if g['condition_id']==cid]
            for g in cg:
                sid=g['source_id']
                if sid not in cached:cached[sid]=matrix(sources[sid]['public_file'])
                a=cached[sid][np.array(g['source_rows'])-1]
                contrast=helper.weberContrast(a[:,6].tolist(),c['background_digital_intensity'],1.0,np.log10,False)
                for i,(rowno,ii,cc) in enumerate(zip(g['source_rows'],a[:,6],contrast)):
                    rv=revmap.get((g['staircase_id'],rowno))
                    w.writerow(dict(condition_id=cid,experiment=g['experiment'],staircase_id=g['staircase_id'],source_id=sid,source_row=rowno,trial_within_staircase=i,probe_position_deg=g['position'],probe_digital_intensity=ii,normalized_contrast=cc,is_detected_reversal=rv is not None,reversal_index=rv['reversal_index'] if rv else '',used_in_final_five=rv['used_in_final_five'] if rv else False,included_for_baseline_candidate=g['included_for_baseline_candidate']));ntrials+=1
                sr=summary[g['staircase_id']]
                thresholds.append(dict(condition_id=cid,participant_id=c['participant_id'],experiment=g['experiment'],staircase_id=g['staircase_id'],probe_position_deg=g['position'],distance_from_flanker_edge_deg=g['position']+.5,detected_reversals=sr['detected_reversals'],threshold_last_five=float(sr['final_five_geometric_threshold']),included_for_baseline_candidate=g['included_for_baseline_candidate']))
                examples[g['staircase_id']]=(contrast,g)
            active=[r for r in modes if r['condition_id']==cid and r['variant']=='thesis_percentiles' and r['retained_for_display']=='True']
            for m in active or [None]:
                endpoints.append(dict(condition_id=cid,participant_id=c['participant_id'],luminance_label_cd_m2=c['luminance_label_cd_m2'],included_in_thesis_main_isf=c['included_in_thesis_main_isf'],endpoint_status='Retained descriptive mode' if m else 'No mode passes display rule',candidate_mode=m['candidate_mode'] if m else '',isf_estimate_cpd=m['median_cpd'] if m else '',isf_lower_cpd=m['quantile_low_cpd'] if m else '',isf_upper_cpd=m['quantile_high_cpd'] if m else '',error_minus_cpd=float(m['median_cpd'])-float(m['quantile_low_cpd']) if m else '',error_plus_cpd=float(m['quantile_high_cpd'])-float(m['median_cpd']) if m else '',accepted_rows=c['thesis_accepted_count'],cluster_count=m['cluster_count'] if m else '',dominant_retained_mode=m['dominant_retained_mode'] if m else False,interval_definition='Within-candidate 0.5th–99.5th percentiles; not population CI',peak_locator=m['peak_locator'] if m else 'source_histogram_spline'))
            condition_rows.append(dict(condition_id=cid,participant_id=c['participant_id'],luminance_label_cd_m2=c['luminance_label_cd_m2'],included_in_thesis_main_isf=c['included_in_thesis_main_isf'],recorded_staircases=len(cg),included_staircases=sum(g['included_for_baseline_candidate'] for g in cg),excluded_staircases=sum(not g['included_for_baseline_candidate'] for g in cg),profile_positions=len({g['position'] for g in cg}),accepted_rows=c['thesis_accepted_count'],retained_endpoint_modes=len(active)))
            cached.clear()
    write_csv(ROOT/'data/processed/staircase_thresholds.csv',thresholds)
    write_csv(ROOT/'data/processed/isf_endpoints.csv',endpoints)
    write_csv(ROOT/'data/processed/condition_summary.csv',condition_rows)
    out=ROOT/'figures/measurements';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for c in cs:
        cg=[g for g in groups if g['condition_id']==c['condition_id']]
        fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        for j,role in enumerate(['base','flanker']):
            # Two actual staircases at the lowest recorded position for orientation.
            selected=sorted([g for g in cg if g['experiment']==role],key=lambda g:(g['position'],g['staircase_id']))[:2]
            for ax,g in zip(axes[:,j],selected):
                contrast,_=examples[g['staircase_id']];rv=[r for r in reversals if r['staircase_id']==g['staircase_id']]
                ax.plot(range(len(contrast)),contrast,color='#6d7c82',lw=1)
                ax.scatter([int(r['trial_within_staircase']) for r in rv],[float(r['normalized_contrast']) for r in rv],s=25,facecolors='none',edgecolors='#0076a8',label='Detected reversal')
                final=[r for r in rv if r['used_in_final_five']=='True']
                ax.scatter([int(r['trial_within_staircase']) for r in final],[float(r['normalized_contrast']) for r in final],s=22,color='#d65f00',label='Final five reversals')
                threshold=float(summary[g['staircase_id']]['final_five_geometric_threshold']);ax.axhline(threshold,color='#d65f00',ls='--',label='Geometric-mean threshold')
                ax.set_xlabel('Trial index within staircase (zero-based)');ax.set_ylabel('Normalized contrast')
                ax.set_title(f"{role.capitalize()}; x = {g['position']:.3f} deg\n{g['staircase_id']} | {'included' if g['included_for_baseline_candidate'] else 'excluded'}",fontsize=9)
        axes[0,0].legend(fontsize=8)
        fig.suptitle(f"Subject {c['thesis_subject']}, luminance label {c['luminance_label_cd_m2']} cd m$^{{-2}}$: staircase examples\nFirst two recorded staircases at lowest position in each condition; full records are in the tables",fontsize=11)
        fig.savefig(out/f"{c['condition_id'].lower()}_staircase_examples.png",dpi=140,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for pid,ax in zip(['P01','P02','P03','P04'],axes.ravel()):
        ee=[r for r in endpoints if r['participant_id']==pid]
        for mode,color,marker in [('1','#d65f00','o'),('2','#0076a8','s')]:
            mm=[r for r in ee if r['candidate_mode']==mode]
            ax.errorbar([r['luminance_label_cd_m2'] for r in mm],[float(r['isf_estimate_cpd']) for r in mm],yerr=np.array([[r['error_minus_cpd'] for r in mm],[r['error_plus_cpd'] for r in mm]]),fmt=marker,color=color,capsize=4,label=f'Candidate {mode}')
        for r in ee:
            if not r['candidate_mode']:ax.annotate('No retained mode',(r['luminance_label_cd_m2'],.5),xytext=(r['luminance_label_cd_m2'],1.3),ha='center',fontsize=8,arrowprops=dict(arrowstyle='-',color='gray'))
        ax.set_xscale('log');ax.set_xticks([10,20,50,100,200],[10,20,50,100,200]);ax.set_xlim(8,250);ax.set_ylim(0,10)
        ax.set_xlabel('Archived luminance label (cd m$^{-2}$)');ax.set_ylabel('ISF (cycles / deg)')
        ax.set_title(f"Subject {int(pid[1:])}"+(' — diagnostic only' if pid=='P03' else ''));ax.legend(fontsize=8)
    fig.suptitle('ISF endpoints from preserved RA-screened fits\nMedian and 0.5th–99.5th percentile range within retained candidate clusters\nSupplied histogram-spline cluster locator; descriptive estimator stability, not population confidence intervals',fontsize=11)
    fig.savefig(out/'isf_endpoints_with_percentile_ranges.png',dpi=160,bbox_inches='tight');plt.close(fig)
    (ROOT/'validation/measurement_export_summary.json').write_text(json.dumps(dict(staircase_trial_rows=ntrials,staircase_thresholds=len(thresholds),endpoint_rows=len(endpoints),retained_endpoints=sum(bool(r['candidate_mode']) for r in endpoints),conditions=17),indent=2)+'\n')
    print('Measurement tables, 17 staircase example figures and ISF endpoint overview created.')
if __name__=='__main__':main()
