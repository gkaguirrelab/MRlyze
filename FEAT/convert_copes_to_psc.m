function convert_copes_to_psc(session_dir,func,gfeatNum,copes)

% Converts copes in feat/stats directories to percent signal change (psc)
%
%   Usage:
%   convert_copes_to_psc(session_dir,func,gfeatNum,copes)
%
%   Written by Andrew S Bock Nov 2015

%% set defaults
if ~exist('copes','var')
    copes = 1:5;
end
%% Find bold directories
b = find_bold(session_dir);
featDir = fullfile(session_dir,b{gfeatNum},[func '.gfeat']);
%% Convert copes
for cc = copes    
    invol = fullfile(featDir,['cope' num2str(cc) '.feat'],'stats','cope1.nii.gz');
    meanvol = fullfile(featDir,['cope' num2str(cc) '.feat'],'mean_func.nii.gz');
    ctmp = load_nifti(invol);
    mtmp = load_nifti(meanvol);
    pscout = ctmp;
    psctmp = (ctmp.vol./mtmp.vol)*100; % convert to percent signal change
    psctmp(psctmp==inf | psctmp==-inf) = nan; % set inf/-inf to nan
    pscout.vol = psctmp;
    save_nifti(pscout,fullfile(featDir,['cope' num2str(cc) '.feat'],'stats',...
        'cope1.psc.nii.gz'));
end