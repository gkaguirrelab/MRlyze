function fsf = create_fsf_struct (EVs_number, Contrasts_number, Ftests_number, condition, templateDir, templateName)
% creates a struct file used to feed parameters into the fsf file generated
% in write_fsf.

% 'Condition' allows for autoselect standard parameters for the given analisys:
% - case blank creates an all zeros fsf struct to be filled in manually;
% - case MelanopsinMR creates a precompiled fsf struct (default values +
% fields to design template) for MelanopsinMR datasets

% Note that this function does not handle EVs, contrasts and Ftest
% generation, but will still need to know the number of each of those.


% June 2016 - written (GF)

%%

switch condition
    case 'blank'   %create an empty fsf_struct (to be filled in manually)
        fsf.version = 0;
        fsf.inmelodic = 0;
        fsf.level  = 0;
        fsf.analysis = 0;
        fsf.relative_yn = 0;
        fsf.help_yn = 0;
        fsf.featwatcher_yn = 0;
        fsf.sscleanup_yn = 0;
        fsf.outputdir = 0;
        fsf.tr = 0;
        fsf.npts = 0;
        fsf.ndelete = 0;
        fsf.tagfirst = 0;
        fsf.multiple = 0;
        fsf.inputtype = 0;
        fsf.filtering_yn = 0;
        fsf.brain_thresh = 0;
        fsf.critical_z = 0;
        fsf.noise = 0;
        fsf.noisear = 0;
        fsf.mc = 0;
        fsf.sh_yn = 0;
        fsf.regunwarp_yn = 0;
        fsf.dwell = 0;
        fsf.te = 0;
        fsf.signallossthresh = 10;
        fsf.unwarp_dir = 0;
        fsf.st = 0;
        fsf.bet_yn = 0;
        fsf.smooth = 0;
        fsf.norm_yn = 0;
        fsf.perfsub_yn = 0;
        fsf.temphp_yn = 0;
        fsf.templp_yn=0;
        fsf.melodic_yn = 0;
        fsf.stats_yn = 0;
        fsf.prewhiten_yn = 0;
        fsf.motionevs = 0;
        fsf.motionevbeta = 0;
        fsf.scriptevsbeta = 0;
        fsf.robust_yn = 0;
        fsf.mixed_yn = 0;
        fsf.evs_orig = 0;
        fsf.evs_real = 0;
        fsf.evs_vox = 0;
        fsf.ncon_orig = 0;
        fsf.ncon_real = 0;
        fsf.nftests_orig= 0;
        fsf.nftests_real= 0;
        fsf.constcol = 0;
        fsf.poststats_yn = 0;
        fsf.threshmask = 0;
        fsf.thresh = 0;
        fsf.prob_thresh = 0;
        fsf.z_thresh = 0;
        fsf.zdisplay = 0;
        fsf.zmin = 0;
        fsf.zmax = 0;
        fsf.rendertype = 0;
        fsf.bgimage = 0;
        fsf.tsplot_yn = 0;
        fsf.reginitial_highres_yn = 0;
        fsf.reginitial_highres_search = 0;
        fsf.reginitial_highres_dof = 0;
        fsf.reghighres_yn = 0;
        fsf.reghighres_search = 0;
        fsf.reghighres_dof = 0;
        fsf.regstandard_yn = 0;
        fsf.alternateReference_yn = 0;
        fsf.regstandard = 0;
        fsf.regstandard_search = 0;
        fsf.regstandard_dof = 0;
        fsf.regstandard_nonlinear_yn = 0;
        fsf.regstandard_nonlinear_warpres = 0;
        fsf.paradigm_hp = 0;
        fsf.totalVoxels = 0 ;
        fsf.ncopeinputs = 0;
        fsf.feat_files= 0;
        fsf.confoundevs = 0;
        % DEFINE EVS, contrast and F-tests
        fsf.alternative_mask = '';
        fsf.init_initial_highres = '';
        fsf.init_highres = '';
        fsf.init_standard = '';
        fsf.overwrite_yn = 0;
        
        
    case 'MelanopsinMR' % create a precompiled fsf_struct (default values + fields to design template)
        
        % FEAT version number
        fsf.version = 6.00;
        
        %Are we in MELODIC?' ...
        fsf.inmelodic = 0;
        
        % Analysis level
        % 1 : First-level analysis
        % 2 : Higher-level analysis
        fsf.level = 1;
        
        % Which stages to run
        % 0 : No first-level analysis (registration and/or group stats only)
        % 7 : Full first-level analysis
        % 1 : Preprocessing
        % 2 : Statistics
        fsf.analysis = 7;
        
        % Use relative filenames
        fsf.relative_yn = 0;
        
        % Balloon help
        fsf.help_yn = 1;
        
        % Run Featwatcher
        fsf.featwatcher_yn = 1;
        
        % Cleanup first-level standard-space images
        fsf.sscleanup_yn = 0;
        
        % Output directory
        fsf.outputdir = '';
        
        % TR(s)
        fsf.tr = 'DESIGN_TR'; %%%%%%%%%%%%%%%%%%%%% Specified in design
        
        % Total volumes
        fsf.npts = 'DESIGN_VOLS'; %%%%%%%%%%%%%%%%%%%%%%%%% specified in design
        
        % Delete volumes
        fsf.ndelete = 0;
        
        %  Perfusion tag/control order
        fsf.tagfirst = 1;
        
        % Number of first-level analyses
        fsf.multiple = 1;  % Not a flag - actual num of analyses
        
        % Higher-level input type
        % 1 : Inputs are lower-level FEAT directories
        % 2 : Inputs are cope images from FEAT directories
        fsf.inputtype = 2;
        
        % Carry out pre-stats processing?
        fsf.filtering_yn = 1;
        
        % Brain/background threshold %
        fsf.brain_thresh = 10;
        
        % Critical z for design efficiency calculation
        fsf.critical_z = 5.3;
        
        % Noise level
        fsf.noise = 0.66;
        
        % Noise AR(1)
        fsf.noisear = 0.34;
        
        % Motion correction
        % 0 : None
        % 1 : MCFLIRT
        fsf.mc = 0;
        
        % Spin-history (currently obsolete)
        fsf.sh_yn = 0;
        
        % B0 fieldmap unwarping?
        fsf.regunwarp_yn = 0;
        
        % EPI dwell time (ms)
        fsf.dwell = 0.7;
        
        % EPI TE (ms)
        fsf.te = 35;
        
        % Signal loss threshold
        fsf.signallossthresh = 10;
        
        % Unwarp direction
        fsf.unwarp_dir = 'y-';
        
        % Slice timing correction
        % 0 : None
        % 1 : Regular up (0, 1, 2, 3, ...)
        % 2 : Regular down
        % 3 : Use slice order file
        % 4 : Use slice timings file
        % 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
        fsf.st = 0;
        
        % Slice timings file
        fsf.st_file='';
        
        % BET brain extraction
        fsf.bet_yn = 0;
        
        % Spatial smoothing FWHM (mm)
        fsf.smooth = 0;
        
        % Intensity normalization
        fsf.norm_yn = 0;
        
        % Perfusion subtraction
        fsf.perfsub_yn = 0;
        
        % Highpass temporal filtering
        fsf.temphp_yn = 0;
        
        % Lowpass temporal filtering
        fsf.templp_yn = 0;
        
        % MELODIC ICA data exploration
        fsf.melodic_yn = 0;
        
        % Carry out main stats?
        fsf.stats_yn = 1;
        
        % Carry out prewhitening?
        fsf.prewhiten_yn = 1;
        
        % Add motion parameters to model
        % 0 : No
        % 1 : Yes
        fsf.motionevs = 0;
        fsf.motionevbeta = '';
        fsf.scriptevsbeta = '';
        
        % Robust outlier detection in FLAME?
        fsf.robust_yn = 0;
        
        % Higher-level modelling
        % 3 : Fixed effects
        % 0 : Mixed Effects: Simple OLS
        % 2 : Mixed Effects: FLAME 1
        % 1 : Mixed Effects: FLAME 1+2
        fsf.mixed_yn = 2;
        
        % Number of EVs
        fsf.evs_orig = EVs_number; %%%%%%%%%%%%%%%%%% find in design file
        fsf.evs_real = EVs_number; %%%%%%%%%%%%%%%%%%% find in design
        fsf.evs_vox = 0; %%%%%%%%%%%%%% find in design
        
        % Number of contrasts
        fsf.ncon_orig = Contrasts_number; %%%%%%%%%%%%%%%%%%%%%
        fsf.ncon_real = Contrasts_number; %%%%%%%%%%%%%%%%%%%%%
        
        % Number of F-tests
        fsf.nftests_orig= Ftests_number; %%%%%%%%%%%%%%%%%%%
        fsf.nftests_real= Ftests_number; %%%%%%%%%%%%%%%%%%%%
        
        % Add constant column to design matrix? (obsolete)
        fsf.constcol = 0;
        
        % Carry out post-stats steps?
        fsf.poststats_yn = 1;
        
        % Pre-threshold masking?
        fsf.threshmask = '';
        
        % Thresholding
        % 0 : None
        % 1 : Uncorrected
        % 2 : Voxel
        % 3 : Cluster
        fsf.thresh = 3;
        
        % P threshold
        fsf.prob_thresh = 0.05;
        
        % Z threshold
        fsf.z_thresh = 2.3;
        
        % Z min/max for colour rendering
        % 0 : Use actual Z min/max
        % 1 : Use preset Z min/max
        fsf.zdisplay = 0;
        
        % Z min in colour rendering
        fsf.zmin = 2;
        
        % Z max in colour rendering
        fsf.zmax = 8;
        
        % Colour rendering type
        % 0 : Solid blobs
        % 1 : Transparent blobs
        fsf.rendertype = 1;
        
        % Background image for higher-level stats overlays
        % 1 : Mean highres
        % 2 : First highres
        % 3 : Mean functional
        % 4 : First functional
        % 5 : Standard space template
        fsf.bgimage = 1;
        
        % Create time series plots' ...
        fsf.tsplot_yn = 1;
        
        % Registration to initial structural
        fsf.reginitial_highres_yn = 0;
        
        % Search space for registration to initial structural
        % 0   : No search
        % 90  : Normal search
        % 180 : Full search
        fsf.reginitial_highres_search = 90;
        
        % Degrees of Freedom for registration to initial structural
        fsf.reginitial_highres_dof = 3;
        
        % Registration to main structural
        fsf.reghighres_yn = 0;
        
        % Search space for registration to main structural
        % 0   : No search
        % 90  : Normal search
        % 180 : Full search
        fsf.reghighres_search = 90;
        
        % Degrees of Freedom for registration to main structural
        fsf.reghighres_dof = 'BBR';
        
        % Registration to standard image?
        fsf.regstandard_yn=1;
        
        % Use alternate reference images?
        fsf.alternateReference_yn = 0;
        
        % Standard image
        fsf.regstandard = 'DESIGN_STANDARD'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Search space for registration to standard space
        % 0   : No search
        % 90  : Normal search
        % 180 : Full search
        fsf.regstandard_search = 90;
        
        % Degrees of Freedom for registration to standard space
        fsf.regstandard_dof = 6;
        
        % Do nonlinear registration to standard space?
        fsf.regstandard_nonlinear_yn = 0;
        
        % Control nonlinear warp field resolution
        fsf.regstandard_nonlinear_warpres = 10;
        
        % High pass filter cutoff
        fsf.paradigm_hp = 100;
        
        % Total voxels
        fsf.totalVoxels = 'DESIGN_TOTAL_VOXELS' ; %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Number of lower-level copes feeding into higher-level analysis
        fsf.ncopeinputs=0;
        
        % 4D AVW data or FEAT directory
        fsf.feat_files= 'DESIGN_FEAT_DIR';
        
        % Add confound EVs text file
        fsf.confoundevs = 0;
        
        %%%%
        %%%% Here is the section where EVs, contrasts and F-test are defined
        %%%%
        
        % Alternative (to BETting) mask image
        fsf.alternative_mask = '';
        
        % Initial structural space registration initialisation transform
        fsf.init_initial_highres = '';
        
        % Structural space reigstration initialisation transform
        fsf.init_highres = '';
        
        % Standard space reigstration initialisation transform
        fsf.init_standard = '';
        
        % For full FEAT analysis: overwrite existing .feat output dir?
        fsf.overwrite_yn = 0;
        
    otherwise
        fprintf ('\nThis condition does not exist yet. Run create_fsf_struct(''blank'') to obtain a blank fsf stuct to fill manually.');
end

%% save out the structs as mat files
save(fullfile(templateDir, [templateName '_fsf.mat']), 'fsf');

