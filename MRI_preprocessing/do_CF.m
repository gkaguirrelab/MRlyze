function do_CF(session_dir,subject_name,runNum,srcfile,trgfile,srcind,trgind,srcROI,trgROI,seedSig,hemi,DoG)


%% Define initial variables
if ~exist('session_dir','var')
    error('No ''session_dir'' defined')
end
if ~exist('subject_name','var')
    error('No ''subject_name'' defined')
end
if ~exist('runNum','var')
    error('No ''runNum'' defined')
end
if ~exist('srcfile','var')
    error('No ''srcfile'' defined')
end
if ~exist('trgfile','var')
    error('No ''trgfile'' defined')
end
if ~exist('srcind','var')
    error('No ''srcind'' defined')
end
if ~exist('trgind','var')
    error('No ''trgind'' defined')
end
if ~exist('srcROI','var')
    srcROI = 'cortex';
end
if ~exist('trgROI','var')
    trgROI = 'V1';
end
if ~exist('seedSig','var')
    seedSig = [3 6 9]; % millimeters
end
if ~exist('hemi','var');
    hemi = 'lh';
end
if ~exist('DoG','var');
    DoG = 1;
end
if ~exist('SUBJECTS_DIR','var')
    SUBJECTS_DIR = getenv('SUBJECTS_DIR');
end
%% Add to log
SaveLogInfo(session_dir,mfilename,session_dir,subject_name,runNum,srcfile,trgfile,srcind,trgind,srcROI,trgROI,seedSig,hemi,DoG)

%% Find bold run directories
d = find_bold(session_dir);
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Found ' num2str(nruns) ' runs']);
disp(['Working on ' d{runNum}]);
%% Load distance info
disp('Loading distance values...');
distmat = load(fullfile(session_dir,[hemi '.' trgROI '.dists.mat']));
trgdists = distmat.allDistances;
disp('done.');

%% Load timecourses
disp('Loading timecourses...');
% load files
src = load_nifti(srcfile);
trg = load_nifti(trgfile);
srcdims = size(src.vol);
trgdims = size(trg.vol);
srctc = reshape(src.vol,srcdims(1)*srcdims(2)*srcdims(3),srcdims(4))';
trgtc = reshape(trg.vol,trgdims(1)*trgdims(2)*trgdims(3),trgdims(4))';
% Pull out relevant timecourses
trgtc = trgtc(:,trgind);
srctc = srctc(:,srcind);
% Set timecourses with very little variation (var<.1) to flat
trgtc = set_to_flat(trgtc);
srctc = set_to_flat(srctc);
disp('done.');
%% Run calc_CF
cd(fullfile(session_dir,d{runNum}));
[cfs] = calc_CF(srctc,trgtc,trgdists,seedSig,DoG);
% Save CF data
if ~exist(fullfile(session_dir,'CFs'),'dir')
    mkdir(fullfile(session_dir,'CFs'));
end
if ~exist(fullfile(session_dir,'CFs',d{runNum}),'dir')
    mkdir(fullfile(session_dir,'CFs',d{runNum}));
end
save(fullfile(session_dir,'CFs',d{runNum},[hemi '.' srcROI '.' trgROI '.cfs.mat']),'cfs');
%% Plot in src volume
disp('Making CF maps...');
% Load ecc and pol maps
if strcmp(trgROI,'V1');
    eccsurf = load_nifti(fullfile(session_dir,[hemi '.ecc.nii.gz']));
    polsurf = load_nifti(fullfile(session_dir,[hemi '.pol.nii.gz']));
elseif strcmp(trgROI,'prf_V1');
    eccsurf = load_nifti(fullfile(session_dir,[hemi '.ecc_pRF.nii.gz']));
    polsurf = load_nifti(fullfile(session_dir,[hemi '.pol_pRF.nii.gz']));
end
% Load in time to HRF peak data
trgpeak = load_nifti(fullfile(session_dir,[hemi '.cortex.avg.copeakt.prfs.nii.gz']));
% Extract values for V1
Eccen = eccsurf.vol(trgind);
Polar = polsurf.vol(trgind);
TrgPeak = trgpeak.vol(trgind);
% Load template file (we'll overwrite the values)
if strcmp(srcROI,'cortex');
    mri = eccsurf;
else
    mri = load_nifti(fullfile(session_dir,d{runNum},'single_TR.nii.gz'));
end
% Save maps
mri.vol = nan(size(mri.vol));
co = mri;
coshiftt = mri;
cotrgpeakt = mri;
varr = mri;
varshiftt = mri;
vartrgpeakt = mri;
var_as_co = mri;
cosig1 = mri;
cosig2 = mri;
cosig3 = mri;
cosig4 = mri;
varsig1 = mri;
varsig2 = mri;
varsig3 = mri;
varsig4 = mri;
disp('Saving correlation and variance explained maps...');
%co and copeakt
co.vol(srcind) = cfs.co;
save_nifti(co,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.co.cfs.nii.gz']));
coshiftt.vol(srcind) = cfs.copeakt;
save_nifti(coshiftt,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.coshiftt.cfs.nii.gz']));
cotrgpeakt.vol(srcind) = TrgPeak(cfs.cocenter);
save_nifti(cotrgpeakt,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.cotrgpeakt.cfs.nii.gz']));
%var and varpeakt
varr.vol(srcind) = cfs.var;
save_nifti(varr,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.var.cfs.nii.gz']));
varshiftt.vol(srcind) = cfs.varpeakt;
save_nifti(varshiftt,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varshiftt.cfs.nii.gz']));
vartrgpeakt.vol(srcind) = TrgPeak(cfs.varcenter);
save_nifti(vartrgpeakt,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.vartrgpeakt.cfs.nii.gz']));
if DoG
    %sig1
    cosig1.vol(srcind) = cfs.cosig1;
    save_nifti(cosig1,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.cosig1.cfs.nii.gz']));
    %sig2
    cosig2.vol(srcind) = cfs.cosig2;
    save_nifti(cosig2,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.cosig2.cfs.nii.gz']));
    %sig3
    cosig3.vol(srcind) = cfs.cosig3;
    save_nifti(cosig3,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.cosig3.cfs.nii.gz']));
    %sig4
    cosig4.vol(srcind) = cfs.cosig4;
    save_nifti(cosig4,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.cosig4.cfs.nii.gz']));
    % For variance, flip Gaussian params so all correlations are positive
    % Find var_as_co values that are negative (these we will flip)
    flip_sigs_ind = cfs.var_as_co;
    flip_sigs_ind(isnan(flip_sigs_ind)) = 0;
    flip_sigs_ind(flip_sigs_ind>0) = 0;
    flip_sigs_ind(flip_sigs_ind<0) = 1;
    flip_sigs_ind = logical(flip_sigs_ind);
    % Create the flipped values (i.e. if positive correlation)
    newsig3 = -cfs.varsig3(flip_sigs_ind);
    newsig4 = -cfs.varsig4(flip_sigs_ind);
    newvar_as_co = -cfs.var_as_co(flip_sigs_ind);
    % Replace with flipped values
    cfs.varsig3(flip_sigs_ind) = newsig3;
    cfs.varsig4(flip_sigs_ind) = newsig4;
    cfs.var_as_co(flip_sigs_ind) = newvar_as_co;
    % Save Gaussian params
    disp('Saving Gaussian parameter maps...');
    var_as_co.vol(srcind) = cfs.var_as_co;
    save_nifti(var_as_co,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.var_as_co.cfs.nii.gz']));
    %sig1
    varsig1.vol(srcind) = cfs.varsig1;
    save_nifti(varsig1,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varsig1.cfs.nii.gz']));
    %sig2
    varsig2.vol(srcind) = cfs.varsig2;
    save_nifti(varsig2,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varsig2.cfs.nii.gz']));
    %sig3
    varsig3.vol(srcind) = cfs.varsig3;
    save_nifti(varsig3,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varsig3.cfs.nii.gz']));
    %sig4
    varsig4.vol(srcind) = cfs.varsig4;
    save_nifti(varsig4,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varsig4.cfs.nii.gz']));
else
    cosig1.vol(srcind) = cfs.sig;
    save_nifti(cosig1,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.sig.cfs.nii.gz']));
end
% Save pol and ecc maps
disp('Saving eccentricity and polar angle maps...');
mri.vol = nan(size(mri.vol));
coecc = mri;
copol = mri;
varecc = mri;
varpol = mri;
coecc.vol(srcind) = Eccen(cfs.cocenter);
save_nifti(coecc,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.coecc.cfs.nii.gz']));
copol.vol(srcind) = Polar(cfs.cocenter);
save_nifti(copol,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.copol.cfs.nii.gz']));
varecc.vol(srcind) = Eccen(cfs.varcenter);
save_nifti(varecc,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varecc.cfs.nii.gz']));
varpol.vol(srcind) = Polar(cfs.varcenter);
save_nifti(varpol,fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.varpol.cfs.nii.gz']));
if ~strcmp(srcROI,'cortex');
    % Project to anatomical space
    disp('Projecting maps to anatomical space...');
    if DoG
        maps = {'co' 'coshiftt' 'cotrgpeakt' 'cosig1' 'cosig2' 'cosig3' 'cosig4' 'coecc' 'copol' ...
            'var' 'varshiftt' 'var_as_co' 'vartrgpeakt' 'varsig1' 'varsig2' 'varsig3' 'varsig4' 'varecc' 'varpol'};
    else
        maps = {'co' 'var' 'var_as_co' 'sig' 'ecc' 'pol'};
    end
    for m = 1:length(maps);
        savename = fullfile(session_dir,[hemi '.' srcROI '.' trgROI '.run' num2str(runNum) '.' maps{m} '.cfs']);
        [~,~] = system(['mri_vol2vol --mov ' savename '.nii.gz ' ...
            '--targ ' fullfile(SUBJECTS_DIR,subject_name,'mri','orig.mgz') ...
            ' --o ' savename '.orig.nii.gz --reg ' ...
            fullfile(session_dir,d{runNum},'brf_bbreg.dat') ...
            ' --interp nearest']);
    end
end
disp('done.');