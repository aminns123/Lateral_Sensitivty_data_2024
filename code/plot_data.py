"""Create clearly labelled diagnostics from archived/recomputed tables, without refitting."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from reproduce_raw import ROOT,read_csv,matrix
COLORS=['#d65f00','#0076a8']

def main():
    out=ROOT/'figures/diagnostics';out.mkdir(exist_ok=True)
    cs=json.loads((ROOT/'config/conditions.json').read_text())
    profiles=read_csv(ROOT/'data/processed/lateral_profiles_recomputed.csv')
    modes=read_csv(ROOT/'data/processed/isf_mode_summaries.csv')
    kde=read_csv(ROOT/'data/processed/isf_kde.csv')
    pts=read_csv(ROOT/'data/processed/representative_profile_points.csv')
    curves=read_csv(ROOT/'data/processed/representative_curves_recomputed.csv')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    for pid in ['P01','P02','P03','P04']:
        conditions=[c for c in cs if c['participant_id']==pid];n=len(conditions)
        status=' — diagnostic only; excluded from main thesis ISF interpretation' if pid=='P03' else ''
        fig,axs=plt.subplots(n,2,figsize=(10,2.35*n),squeeze=False,layout='constrained')
        for c,(left,right) in zip(conditions,axs):
            rows=[r for r in profiles if r['condition_id']==c['condition_id']]
            x=[float(r['distance_from_flanker_edge_deg']) for r in rows]
            left.plot(x,[float(r['base_sensitivity']) for r in rows],'o-',label='Base',color='#555555')
            left.plot(x,[float(r['flanker_sensitivity']) for r in rows],'s-',label='Flanker',color=COLORS[1])
            right.errorbar(x,[float(r['log_sensitivity_ratio']) for r in rows],yerr=[float(r['spread'])/2 for r in rows],fmt='o',color=COLORS[1],capsize=2)
            right.axhline(0,color='gray',ls='--',lw=.8)
            for ax in [left,right]:ax.set_title(f"Luminance label {c['luminance_label_cd_m2']} cd m$^{{-2}}$");ax.set_xlabel('Distance from flanker edge (deg)')
            left.set_ylabel('Sensitivity\n(1 / normalized contrast)');right.set_ylabel('Log sensitivity ratio')
        axs[0,0].legend(fontsize=9)
        fig.suptitle(f'Subject {int(pid[1:])}: lateral sensitivity{status}\nLast five reversals; bars = recorded pairwise spread / 2 (not confidence intervals)',fontsize=12)
        fig.savefig(out/f'{pid.lower()}_lateral_profiles.png',dpi=150,bbox_inches='tight');plt.close(fig)
        for selection,label in [('source_plotting_index_rule','Supplied plotting index rule'),('row_linked_max_fms','Maximum saved FMS within candidate cluster')]:
            fig,axs=plt.subplots(n,3,figsize=(12,2.5*n),squeeze=False,layout='constrained')
            for c,axes in zip(conditions,axs):
                cid=c['condition_id'];a=matrix(f'data/source/saved_fits/{cid.lower()}.tsv.gz');d=a[(a[:,0]>=.8)&(a[:,0]<=1),9]
                axes[0].hist(d,bins=31,density=True,histtype='step',color='gray')
                kd=[r for r in kde if r['condition_id']==cid]
                axes[0].plot([float(r['frequency_cpd']) for r in kd],[float(r['density_per_cpd']) for r in kd],color='black',lw=1.3)
                active=[r for r in modes if r['condition_id']==cid and r['variant']=='thesis_percentiles' and r['retained_for_display']=='True']
                axes[0].set_title(f"{c['luminance_label_cd_m2']} cd m$^{{-2}}$; accepted n = {len(d):,}")
                axes[0].set_xlabel('Intrinsic spatial frequency (cycles / deg)');axes[0].set_ylabel('Probability density')
                for k,ax in enumerate(axes[1:]):
                    if k>=len(active):
                        ax.set_facecolor('#f3f3f3');ax.text(.5,.5,'No additional mode retained',ha='center',va='center',transform=ax.transAxes,color='#666666');ax.set_xticks([]);ax.set_yticks([]);continue
                    m=active[k];mode=m['candidate_mode'];color=COLORS[k]
                    axes[0].axvline(float(m['spline_peak_cpd']),color=color,ls='--',lw=1)
                    pp=[r for r in pts if r['condition_id']==cid and r['candidate_mode']==mode and r['selection']==selection]
                    cc=[r for r in curves if r['condition_id']==cid and r['candidate_mode']==mode and r['selection']==selection]
                    ax.errorbar([float(r['distance_from_flanker_edge_deg']) for r in pp],[float(r['log_sensitivity_ratio']) for r in pp],yerr=[float(r['errorbar_halfwidth_source_convention']) for r in pp],fmt='o',ms=3,capsize=2,color=color)
                    ax.plot([float(r['distance_from_flanker_edge_deg']) for r in cc],[float(r['log_sensitivity_ratio']) for r in cc],color=color)
                    ax.axhline(0,color='gray',ls='--',lw=.7)
                    ax.set_title(f"Candidate {mode}\nSpline peak {float(m['spline_peak_cpd']):.2f} cycles / deg",fontsize=10)
                    ax.set_xlabel('Distance from flanker edge (deg)');ax.set_ylabel('Log sensitivity ratio')
            fig.suptitle(f'Subject {int(pid[1:])}: new diagnostic fits{status}\n{label}; assumed row-to-input link; 0.5th–99.5th percentile mode screening',fontsize=11)
            fig.savefig(out/f'{pid.lower()}_{selection}.png',dpi=150,bbox_inches='tight');plt.close(fig)
    print('Created 12 labelled diagnostic figures. Reference PNGs remain unchanged.')
if __name__=='__main__':main()
