"""Compare exported R PPI designs/fits with Nilearn's OLS and AR operators.

The design is the input to this oracle: this does not validate its construction.
"""
import hashlib, json, pathlib, sys
import nilearn, numpy as np, scipy
from nilearn.glm.regression import OLSModel, ARModel

fixture = pathlib.Path(sys.argv[1]); output = pathlib.Path(sys.argv[2]); tol=2e-8
cases = {}; failures=[]; source_hashes=None
for name in ('iid','censored_iid','rank_deficient','fixed_ar','scaled_columns'):
    d=json.loads((fixture / (name+'.json')).read_text())
    if source_hashes is None: source_hashes=d['source_hashes']
    elif source_hashes != d['source_hashes']: failures.append(name+': source hash mismatch')
    X=np.asarray(d['X'],float); Y=np.asarray(d['Y'],float); keep=np.asarray(d['keep'],bool)
    boundary_effect=None; scales=None; raw_rank=None
    if name=='fixed_ar':
        # Nilearn leaves first AR(1) row unscaled. fmriAR exact_first='ar1'
        # uses sqrt(1-phi^2), separately at each run boundary.
        runs=np.asarray(d['runs']); phi=float(d['phi']); Xw=X.copy(); Yw=Y.copy()
        Xplain=X.copy(); Yplain=Y.copy()
        for r in dict.fromkeys(runs):
            idx=np.flatnonzero(runs==r); op=ARModel(X[idx],rho=np.array([phi]))
            Xw[idx]=op.whiten(X[idx]); Yw[idx]=op.whiten(Y[idx])
            Xplain[idx]=Xw[idx]; Yplain[idx]=Yw[idx]
            Xw[idx[0]] *= np.sqrt(1-phi**2); Yw[idx[0]] *= np.sqrt(1-phi**2)
        boundary_effect=float(np.max(np.abs(OLSModel(Xplain).fit(Yplain).theta - OLSModel(Xw).fit(Yw).theta)))
        X,Y=Xw,Yw
    else: X,Y=X[keep],Y[keep]
    if name=='scaled_columns':
        raw_rank=int(np.linalg.matrix_rank(X)); scales=np.linalg.norm(X,axis=0); X=X/scales
    model=OLSModel(X); fit=model.fit(Y); rank=np.linalg.matrix_rank(X)
    # Nilearn OLS dispersion divides by n-p even if rank deficient.
    # Use its residuals with n-rank for the separately flagged aliased case.
    sigma2=np.sum(fit.whitened_residuals**2,axis=0)/(X.shape[0]-rank)
    rt=np.asarray(d['t_contrast'],float); rf=np.asarray(d['f_contrast'],float)
    # Independently verify condition-named A-B averages the two run terms.
    names=d["design_names"]; ia=[i for i,x in enumerate(names) if x.startswith("interaction:A:")]; ib=[i for i,x in enumerate(names) if x.startswith("interaction:B:")]
    expected=np.zeros((2,len(names))); expected[0,ia]=1/len(ia); expected[1,ib]=1/len(ib)
    if len(ia)!=2 or len(ib)!=2 or not np.allclose(rt[0],expected[0]-expected[1]) or not np.allclose(rf,expected): failures.append(name+": named A-B/run contrast mismatch")
    if name=='scaled_columns' and raw_rank >= int(rank): failures.append('scaled_columns: expected raw numerical rank loss')
    t=fit.Tcontrast(rt if scales is None else rt/scales,dispersion=sigma2)
    f=fit.Fcontrast(rf if scales is None else rf/scales,dispersion=sigma2)
    beta=np.asarray(d['coefficients'],float); se=np.asarray(d['standard_errors'],float)
    se_ref=np.sqrt(np.outer(np.diag(model.normalized_cov_beta),sigma2))
    beta_ref=fit.theta if scales is None else fit.theta/scales[:,None]
    if scales is not None: se_ref=se_ref/scales[:,None]
    expected_estimable=np.abs(1-np.diag(np.linalg.pinv(X) @ X)) < 1e-7
    if not np.array_equal(np.isfinite(beta).all(axis=1),expected_estimable): failures.append(name+': coefficient estimability mismatch')
    def err(a,b):
        a=np.asarray(a,float).ravel(); b=np.asarray(b,float).ravel()
        good=np.isfinite(a)&np.isfinite(b)
        return float(np.max(np.abs(a[good]-b[good]))) if np.any(good) else None
    if scales is None: beta_error=err(beta,beta_ref); se_error=err(se,se_ref)
    else:
        beta_error=float(np.max(np.abs(beta-beta_ref)/np.maximum(1,np.abs(beta_ref))))
        se_error=float(np.max(np.abs(se-se_ref)/np.maximum(1,np.abs(se_ref))))
    vals={'beta':beta_error,'variance':err(d['sigma2'],sigma2),'coefficient_se':se_error,
          'contrast_estimate':err(d['t_estimate'],t.effect),'contrast_se':err(d['t_se'],t.sd),
          't':err(d['t_stat'],t.t),'F':err(d['f_stat'],f.F),
          'small_contrast_t':err(d['small_t_stat'],t.t),'small_contrast_F':err(d['small_f_stat'],f.F)}
    if rank!=d['rank'] or X.shape[0]-rank != d['df']: failures.append(name+': rank/df mismatch')
    for k,v in vals.items():
        if v is None or v > tol: failures.append(f'{name} {k}: {v}')
    cases[name]={'n':int(X.shape[0]),'columns':int(X.shape[1]),'rank':int(rank),'df':int(X.shape[0]-rank),
                 'max_abs_error':vals,'nilearn_native_dispersion_max_difference':err(fit.dispersion,sigma2),
                 'df_adjusted_externally':name=='rank_deficient','phi':d['phi'] if name=='fixed_ar' else None,
                 'uncorrected_boundary_beta_max_difference':boundary_effect,
                 'raw_unscaled_numerical_rank':raw_rank,
                 'coefficient_error_metric':'max relative error with unit floor' if scales is not None else 'max absolute error'}
# Grouped estimated-noise fit: hold each exported group's estimated AR(1)
# operator fixed and compare each target independently with Nilearn.
d=json.loads((fixture/'grouped_ar.json').read_text())
if source_hashes != d['source_hashes']: failures.append('grouped_ar: source hash mismatch')
X=np.asarray(d['X'],float); Y=np.asarray(d['Y'],float); runs=np.asarray(d['runs']); phis=np.asarray(d['group_phi'],float)
rt=np.asarray(d['t_contrast'],float); rf=np.asarray(d['f_contrast'],float)
names=d['design_names']; ia=[i for i,x in enumerate(names) if x.startswith('interaction:A:')]; ib=[i for i,x in enumerate(names) if x.startswith('interaction:B:')]
expected=np.zeros((2,len(names))); expected[0,ia]=1/len(ia); expected[1,ib]=1/len(ib)
if len(ia)!=2 or len(ib)!=2 or not np.allclose(rt[0],expected[0]-expected[1]) or not np.allclose(rf,expected): failures.append('grouped_ar: named A-B/run contrast mismatch')
if abs(phis[0]-phis[1]) < .05 or not np.isclose(phis[0],phis[2]): failures.append('grouped_ar: expected two distinct group operators')
beta=[];se=[];sig=[]; te=[];ts=[];tt=[];ff=[]; ranks=[]
for j,phi in enumerate(phis):
    xw=X.copy(); yw=Y[:,j:j+1].copy()
    for r in dict.fromkeys(runs):
        idx=np.flatnonzero(runs==r); op=ARModel(X[idx],rho=np.array([phi]))
        xw[idx]=op.whiten(X[idx]); yw[idx]=op.whiten(Y[idx,j:j+1])
        xw[idx[0]]*=np.sqrt(1-phi**2); yw[idx[0]]*=np.sqrt(1-phi**2)
    model=OLSModel(xw); fit=model.fit(yw); rank=int(np.linalg.matrix_rank(xw)); ranks.append(rank)
    var=np.sum(fit.whitened_residuals**2,axis=0)/(len(xw)-rank)
    t=fit.Tcontrast(rt,dispersion=var); f=fit.Fcontrast(rf,dispersion=var)
    beta.append(fit.theta[:,0]);se.append(np.sqrt(np.diag(model.normalized_cov_beta)*var[0]));sig.append(var[0]);te.append(np.asarray(t.effect).item());ts.append(np.asarray(t.sd).item());tt.append(np.asarray(t.t).item());ff.append(np.asarray(f.F).item())
def group_err(x,y): return float(np.max(np.abs(np.asarray(x,float).ravel()-np.asarray(y,float).ravel())))
gvals={'beta':group_err(d['coefficients'],np.column_stack(beta)),'variance':group_err(d['sigma2'],sig),
       'coefficient_se':group_err(d['standard_errors'],np.column_stack(se)),
       'contrast_estimate':group_err(d['t_estimate'],te),'contrast_se':group_err(d['t_se'],ts),
       't':group_err(d['t_stat'],tt),'F':group_err(d['f_stat'],ff),
       'small_contrast_t':group_err(d['small_t_stat'],tt),'small_contrast_F':group_err(d['small_f_stat'],ff)}
if not np.array_equal(np.asarray(d['rank']),ranks) or not np.array_equal(np.asarray(d['df']),len(X)-np.asarray(ranks)): failures.append('grouped_ar: rank/df mismatch')
for k,v in gvals.items():
    if not np.isfinite(v) or v>tol: failures.append(f'grouped_ar {k}: {v}')
cases['grouped_ar']={'n':len(X),'columns':X.shape[1],'target_ranks':ranks,'max_abs_error':gvals,
                     'groups':d['groups'],'estimated_phi_by_target':phis.tolist(),
                     'qualification':'conditional comparison holding fitted group plans fixed; not noise calibration'}
result={'source_hashes':source_hashes,'versions':{'nilearn':nilearn.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
        'scope':'Regression and contrast given exported design; no construction validation.',
        'fixed_ar_boundary':'Nilearn ARModel first row left unchanged; corrected to sqrt(1-phi^2) per run before OLSModel fitting.',
        'rank_deficient_note':'Nilearn 0.14.1 OLSModel.fit divides residual sum of squares by n - number_of_columns; rank-corrected dispersion supplied to its contrast methods.',
        'tolerance':tol,'cases':cases,'failures':failures}
output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
if failures: raise SystemExit(1)
