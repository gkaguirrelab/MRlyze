function contrasts = MakeRegressorsFSL(runName, TR, outputDir)
%  contrasts = MakeRegressorsFSL(runName, TR, outputDir)
%
% Generates contrasts depending on the passed 'runName' variable. This is
% currently hard-coded for a few select cases, including the Isochromatic
% localizer and square wave luxotonic data. For both, the temporal
% structure was always the same.
%
% Saves things for later use in a .txt file in the run directory

if ~isempty(strfind(runName, 'IsochromaticLocalizer'))
    nTRs = 150;
    t = 0:TR:(nTRs-1)*TR;
    
    % Construct a fine temporal model
    dt = 0.1;
    t_highres = 0:0.1:nTRs*TR-0.1;
    x_highres1 = square(2*pi*1/25*t_highres-pi);
    x_highres1(x_highres1<0) = 0;
    
    % Obtain the HRF at the high res model
    hrf_highres = spm_hrf(dt);
    contrast_highres=conv(x_highres1,hrf_highres);
    % Truncate
    contrast_highres = contrast_highres(1:length(x_highres1));
    
    % Resample for our TRs
    contrast = resample(contrast_highres,1,TR/dt)';
    
    % Save out
    dlmwrite(fullfile(outputDir,['cov_ON.txt']),contrast,'delimiter','\t');
    
    contrasts = contrast;
elseif ~isempty(strfind(runName, 'Isochromatic')) || ~isempty(strfind(runName, 'Melanopsin')) || ~isempty(strfind(runName, 'LMS'))
    nTRs = 156;
    t = 0:TR:(nTRs-1)*TR;
    
    bufferTRs = 10;
    % Construct a fine temporal model
    dt = 2;
    t_highres = 0:dt:(nTRs+bufferTRs)*TR-0.1;
    % Obtain the HRF at the high res model
    hrf_highres = spm_hrf(dt);
    
    %% Contrast 3 (sustained ON)
    x_highres1 = square(2*pi*1/48*t_highres);
    x_highres1(x_highres1<0) = 0;
    
    % Also delete the first 24 seconds
    tmp = find(diff(x_highres1));
    x_highres1(1:tmp(1)) = 0;
    contrast3_highres=conv(x_highres1,hrf_highres);
    
    % Resample for our TRs
    contrast3 = resample(contrast3_highres,1,TR/dt)';
    
    % Truncate based on the number of TRs. That way we avoid edge effects
    % from resampling at the end
    contrast3 = contrast3(1:nTRs);
    contrast3 = contrast3./max(contrast3);
    
    %% Contrast 1 (Transient ON)
    % Find 'on' response
    tmp = [0 diff(x_highres1)];
    x_highres2 = double(tmp == 1);
    contrast1_highres=conv(x_highres2,hrf_highres);
    
    % Resample for our TRs
    contrast1 = resample(contrast1_highres,1,TR/dt)';
    
    % Truncate based on the number of TRs. That way we avoid edge effects
    % from resampling at the end
    contrast1 = contrast1(1:nTRs);
    contrast1 = contrast1./max(contrast1);

    %% Contrast 2 (Transient OFF)
    % Find 'off' response
    tmp = [0 diff(x_highres1)];
    x_highres3 = double(tmp == -1);
    contrast2_highres=conv(x_highres3,hrf_highres);
    
    % Resample for our TRs
    contrast2 = resample(contrast2_highres,1,TR/dt)';
    
    % Truncate based on the number of TRs. That way we avoid edge effects
    % from resampling at the end
    contrast2 = contrast2(1:nTRs);
    contrast2 = contrast2./max(contrast2);
    
    dlmwrite(fullfile(outputDir,['cov_TransientON.txt']),contrast1,'delimiter','\t');
    dlmwrite(fullfile(outputDir,['cov_TransientOFF.txt']),contrast2,'delimiter','\t');
    dlmwrite(fullfile(outputDir,['cov_SustainedON.txt']),contrast3,'delimiter','\t');
    contrasts = [contrast1 contrast2 contrast3];
end
