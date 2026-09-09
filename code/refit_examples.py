"""Replay representative fits and audit a fixed sample; outputs are newly computed."""
from pathlib import Path
import json,warnings
import numpy as np
import fit_helpers as helper
from reproduce_raw import matrix,write_csv
ROOT=Path(__file__).resolve().parents[1]

def refit(x,y,spread,peak,frequency_bounds=None,source_exclusions=()):
    if source_exclusions:
        keep=np.array([i for i in range(len(x)) if i not in source_exclusions]);x,y,spread=x[keep],y[keep],spread[keep]
    low,high=frequency_bounds or (peak-1,peak+1)
    lower=[-1.,-np.pi,low*2*np.pi,.001];upper=[1.,np.pi,high*2*np.pi,1.]
    initial=[-.2,0.,peak*2*np.pi,.2]
    _,_,pars,_=helper.fit_to_ISF(x.tolist(),y.tolist(),helper.greenFUNC_phase,initial,lower,upper,1_000_000,'trf')
    residual=helper.goodnessFit_Residual_Analysis(x.tolist(),y.tolist(),1,helper.greenFUNC_phase,**pars)
    try:removed=helper.largestNvalues(residual[1].copy(),2)
    except IndexError:removed=[]
    # Preserve supplied behavior: these indices refer to the filtered residual array.
    keep=np.array([i for i in range(len(x)) if i not in removed]);xr,yr=x[keep],y[keep];sigma=np.clip(spread,.05,None)[keep]
    starts=np.arange(1,9,3)*2*np.pi if frequency_bounds else np.arange(peak-1,peak+2,1)*2*np.pi
    epsilon=max((upper[2]-lower[2])*1e-10,np.finfo(float).eps)
    starts=np.clip(starts,lower[2]+epsilon,upper[2]-epsilon)
    result=helper.MLE_range_frequency(helper.greenFUNC_phase,xr,yr,dict(zip(['x0','x1','f1','x3'],initial)),sigma,list(zip(lower,upper)),starts,'Powell',1_000_000,idx=['f1'])
    xx=np.linspace(float(x.min()),float(x.max()),400);yy=helper.greenFUNC_phase(xx,*result.x)
    score=helper.frequency_match_score(xr,yr,xx,yy)
    return result,xx,yy,removed,float(score)

def main():
    requests=json.loads((ROOT/'config/representative_requests.json').read_text());conditions=json.loads((ROOT/'config/conditions.json').read_text())
    all_results=[];points=[];curves=[];audit=[]
    for c in conditions:
        cid=c['condition_id'];saved=matrix(f'data/source/saved_fits/{cid.lower()}.tsv.gz')
        pending=[dict(r,kind='representative') for r in requests if r['condition_id']==cid]
        for index in [0,1,1234,25000,49999]:pending.append(dict(condition_id=cid,source_input_index=index,candidate_mode=0,selection='fixed_sample_current_helper',spline_peak_cpd=2.,kind='sample'))
        for task in pending:
            i=task['source_input_index'];stem=f'data/source/example_inputs/{cid.lower()}_{i:05d}'
            profile=matrix(stem+'_profile.tsv.gz');spread=matrix(stem+'_spread.tsv.gz');x,y=profile.T
            out=dict(condition_id=cid,candidate_mode=task['candidate_mode'],selection=task['selection'],assumed_input_index=i,saved_fit_row=i+1,saved_f_n_ra_cpd=saved[i,9],saved_fms_ra=saved[i,0],optimizer_success=False,refit_f_n_cpd='',refit_fms='',amplitude='',phase_rad='',k_n_rad_per_deg='',decay_per_deg='',removed_point_indices='',error='',warning_count=0)
            try:
                peak=task['spline_peak_cpd'];exclusions=[]
                if task['kind']=='sample':
                    guesses={'P01':[2,2,2,3,8],'P02':[2,2,3,3],'P03':[2,2,2,2],'P04':[2,2,3,5]}
                    same=[r for r in conditions if r['participant_id']==c['participant_id']];j=next(j for j,r in enumerate(same) if r['condition_id']==cid);peak=guesses[c['participant_id']][j]
                    if cid=='P01_L025':exclusions=[0,1,2]
                    if cid=='P01_L050':exclusions=[0,1,8,12]
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    fit,xx,yy,removed,score=refit(x,y,spread[:,1],peak,(.5,10) if task['kind']=='sample' else None,exclusions)
                out.update(optimizer_success=bool(fit.success),refit_f_n_cpd=float(fit.x[2]/(2*np.pi)),refit_fms=score,amplitude=float(fit.x[0]),phase_rad=float(fit.x[1]),k_n_rad_per_deg=float(fit.x[2]),decay_per_deg=float(fit.x[3]),removed_point_indices=';'.join(map(str,removed)),warning_count=len(caught))
                if task['kind']=='representative':
                    for j,(xp,yp,sd) in enumerate(zip(x,y,spread[:,1])):points.append(dict(condition_id=cid,candidate_mode=task['candidate_mode'],selection=task['selection'],assumed_input_index=i,point_index=j,probe_position_deg=xp,distance_from_flanker_edge_deg=xp+.5,log_sensitivity_ratio=yp,spread_recorded=sd,errorbar_halfwidth_source_convention=sd/2))
                    for xp,yp in zip(xx,yy):curves.append(dict(condition_id=cid,candidate_mode=task['candidate_mode'],selection=task['selection'],distance_from_flanker_edge_deg=xp+.5,log_sensitivity_ratio=yp))
            except Exception as error:out['error']=type(error).__name__+': '+str(error)
            if task['kind']=='sample':audit.append(out)
            else:all_results.append(out)
        print(cid,'example fits completed',flush=True)
    write_csv(ROOT/'validation/representative_refits.csv',all_results)
    write_csv(ROOT/'validation/sampled_fit_checks.csv',audit)
    write_csv(ROOT/'data/processed/representative_profile_points.csv',points)
    write_csv(ROOT/'data/processed/representative_curves_recomputed.csv',curves)
    matching=int(sum(r['refit_f_n_cpd']!='' and round(r['refit_f_n_cpd'],4)==r['saved_f_n_ra_cpd'] and round(r['refit_fms'],4)==r['saved_fms_ra'] for r in audit))
    summary={'sampled_fits':len(audit),'sampled_errors':sum(bool(r['error']) for r in audit),'sampled_optimizer_successes':sum(r['optimizer_success'] for r in audit),'sampled_frequency_and_fms_match_4dp':matching,'representative_refits':len(all_results),'representative_errors':sum(bool(r['error']) for r in all_results),'representative_optimizer_successes':sum(r['optimizer_success'] for r in all_results),'historical_row_to_input_identity_claimed':False,'residual_removal_behavior':'source filtered-residual-array indices applied to original profile; preserved and disclosed','parameter_mapping':'ordered names; original fitting script uses unordered set names'}
    (ROOT/'validation/refit_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
