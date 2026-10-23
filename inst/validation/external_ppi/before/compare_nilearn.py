"""Compare exported R PPI designs/fits with Nilearn's OLS and AR operators.

The design is the input to this oracle: this does not validate its construction.
"""
import hashlib, json, pathlib, sys
import nilearn, numpy as np, scipy
from nilearn.glm.regression import OLSModel, ARModel

fixture = pathlib.Path(sys.argv[1]); output = pathlib.Path(sys.argv[2]); tol=2e-8
cases = {}; failures=[]
for name in ('iid','censored_iid','rank_deficient','fixed_ar'):
    d=json.loads((fixture / (name+'.json')).read_text())
    X=np.asarray(d['X'],float); Y=np.asarray(d['Y'],float); keep=np.asarray(d['keep'],bool)
    boundary_effect=None
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
    model=OLSModel(X); fit=model.fit(Y); rank=np.linalg.matrix_rank(X)
    # Nilearn OLS dispersion divides by n-p even if rank deficient.
    # Use its residuals with n-rank for the separately flagged aliased case.
    sigma2=np.sum(fit.whitened_residuals**2,axis=0)/(X.shape[0]-rank)
    rt=np.asarray(d['t_contrast'],float); rf=np.asarray(d['f_contrast'],float)
    # Independently verify condition-named A-B averages the two run terms.
    names=d["design_names"]; ia=[i for i,x in enumerate(names) if x.startswith("interaction:A:")]; ib=[i for i,x in enumerate(names) if x.startswith("interaction:B:")]
    expected=np.zeros((2,len(names))); expected[0,ia]=1/len(ia); expected[1,ib]=1/len(ib)
    if len(ia)!=2 or len(ib)!=2 or not np.allclose(rt[0],expected[0]-expected[1]) or not np.allclose(rf,expected): failures.append(name+": named A-B/run contrast mismatch")
    t=fit.Tcontrast(rt,dispersion=sigma2); f=fit.Fcontrast(rf,dispersion=sigma2)
    beta=np.asarray(d['coefficients'],float); se=np.asarray(d['standard_errors'],float)
    se_ref=np.sqrt(np.outer(np.diag(model.normalized_cov_beta),sigma2))
    expected_estimable=np.abs(1-np.diag(np.linalg.pinv(X) @ X)) < 1e-7
    if not np.array_equal(np.isfinite(beta).all(axis=1),expected_estimable): failures.append(name+': coefficient estimability mismatch')
    def err(a,b):
        a=np.asarray(a,float).ravel(); b=np.asarray(b,float).ravel()
        good=np.isfinite(a)&np.isfinite(b)
        return float(np.max(np.abs(a[good]-b[good]))) if np.any(good) else None
    vals={'beta':err(beta,fit.theta),'variance':err(d['sigma2'],sigma2),'coefficient_se':err(se,se_ref),
          'contrast_estimate':err(d['t_estimate'],t.effect),'contrast_se':err(d['t_se'],t.sd),
          't':err(d['t_stat'],t.t),'F':err(d['f_stat'],f.F)}
    if rank!=d['rank'] or X.shape[0]-rank != d['df']: failures.append(name+': rank/df mismatch')
    for k,v in vals.items():
        if v is None or v > tol: failures.append(f'{name} {k}: {v}')
    cases[name]={'n':int(X.shape[0]),'columns':int(X.shape[1]),'rank':int(rank),'df':int(X.shape[0]-rank),
                 'max_abs_error':vals,'nilearn_native_dispersion_max_difference':err(fit.dispersion,sigma2),
                 'df_adjusted_externally':name=='rank_deficient','phi':d['phi'] if name=='fixed_ar' else None,
                 'uncorrected_boundary_beta_max_difference':boundary_effect}
result={'versions':{'nilearn':nilearn.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
        'scope':'Regression and contrast given exported design; no construction validation.',
        'fixed_ar_boundary':'Nilearn ARModel first row left unchanged; corrected to sqrt(1-phi^2) per run before OLSModel fitting.',
        'rank_deficient_note':'Nilearn 0.14.1 OLSModel.fit divides residual sum of squares by n - number_of_columns; rank-corrected dispersion supplied to its contrast methods.',
        'tolerance':tol,'cases':cases,'failures':failures}
output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
if failures: raise SystemExit(1)
