% Run in a private directory with pinned official SPM functions on path.
n = 40; SPM.pwd = '/tmp/multippi-spm'; SPM.xY.RT = 1; SPM.xBF.dt = 1; SPM.xBF.T0 = 1;
SPM.xX.W = eye(n); SPM.Sess.row = 1:n;
SPM.Sess.U(1).name = {'A'}; SPM.Sess.U(1).u = [zeros(32,1); repmat([0;1], n/2, 1)];
xY.Sess = 1; xY.u = sin((1:n)'/4); xY.X0 = ones(n,1); xY.name = 'synthetic';
try
  ppi = spm_peb_ppi(SPM,'ppi',xY,[1 1 1],'synthetic',0);
  PSY = ppi.psy.u * ppi.psy.w(:); PSY = spm_detrend(PSY);
  csvwrite('/tmp/multippi-spm/actual_ppi.csv', [xY.u(:), PSY(:), ppi.Y(:), ppi.P(:), ppi.xn(:), ppi.ppi(:)]);
  csvwrite('/tmp/multippi-spm/actual_hrf.csv', spm_hrf(SPM.xBF.dt));
  fprintf('SPM_PPI_SUCCESS %d\n',numel(ppi.ppi));
catch err
  fprintf('SPM_PPI_ERROR %s\n',err.message);
  for i=1:numel(err.stack), fprintf('  %s:%d\n',err.stack(i).name,err.stack(i).line); end
  rethrow(err);
end
