function save_surface_maps(session_dir,hemi,func,rr,srcind)
seedSig = .5:.5:10;
distmat = load(fullfile(session_dir,[hemi '.V1.dists.mat']));
trgdists = distmat.allDistances;
V1 = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
trgind = V1.vol<=1 & V1.vol >=-1;
load([hemi '.' func '.cfs.mat']);
% Load ecc and pol maps
eccsurf = load_nifti(fullfile(session_dir,[hemi '.ecc.nii.gz']));
Eccen = eccsurf.vol(trgind);
polsurf = load_nifti(fullfile(session_dir,[hemi '.pol.nii.gz']));
Polar = polsurf.vol(trgind);
% Save co and sig maps
surf = eccsurf;
surf.vol = nan(size(surf.vol));
co = surf;
varr = surf;
sig = surf;
ecc = surf;
pol = surf;
co.vol(srcind) = cfs.var_as_corr;
save_nifti(co,fullfile(session_dir,[hemi '_' func '_prf_cortical_run' num2str(rr) '_co_surf.nii.gz']));
varr.vol(srcind) = cfs.var;
save_nifti(varr,fullfile(session_dir,[hemi '_' func '_prf_cortical_run' num2str(rr) '_var_surf.nii.gz']));
sig.vol(srcind) = cfs.sig;
save_nifti(sig,fullfile(session_dir,[hemi '_' func '_prf_cortical_run' num2str(rr) '_sig_surf.nii.gz']));
% Save pol and ecc maps
seedX = 1:length(trgdists);
[xList,~] = ndgrid(seedX,seedSig);
xList = xList(:);
for v = 1:length(cfs.varseed);
    ecc.vol(srcind(v)) = Eccen(xList(cfs.varseed(v)));
    pol.vol(srcind(v)) = Polar(xList(cfs.varseed(v)))-pi/2;
end
save_nifti(ecc,fullfile(session_dir,[hemi '_' func '_prf_cortical_run' num2str(rr) '_ecc_surf.nii.gz']));
save_nifti(pol,fullfile(session_dir,[hemi '_' func '_prf_cortical_run' num2str(rr) '_pol_surf.nii.gz']));
