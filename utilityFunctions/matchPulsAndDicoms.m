function [seriesPulsfile, seriesName, minTimeDiff] = matchPulsAndDicoms(dicomDir,pulsDir,OLD)
%
%
if ~exist('OLD','var')
    OLD = 0; % for the old format of .puls files (see below)
end

%% Get the acquistion time from the dicom files
series.List = dir(fullfile(dicomDir,'*'));
series.Names = {series.List(:).name}';
isHidden = strncmp('.',series.Names,1); % remove hidden files from the list
series.List(isHidden) = [];
numSeries = length(series.List);
display('Reading dicom files');
for s=1:numSeries
    thisSeries = series.List(s).name;
	display(['     series ' num2str(s) ': ' thisSeries]);
    thisSeries(thisSeries=='.') = '_';
    dicomList = dir(fullfile(dicomDir,series.List(s).name,'*'));
    dicomNames = {dicomList(:).name}';
    isHidden = strncmp('.',dicomNames,1); % remove hidden files from the list
    dicomList(isHidden) = [];
    nTRs = length(dicomList);
    dicom.(thisSeries).TR_all = nan(nTRs,1); % TR from each image header (msec)
    
    for t = 1:nTRs
        dicom.(thisSeries).info = dicominfo(fullfile(dicomDir,series.List(s).name,dicomList(t).name));
        dicom.(thisSeries).TR_all(t,1) = dicom.(thisSeries).info.RepetitionTime;
        dicom.(thisSeries).timeStr = dicom.(thisSeries).info.AcquisitionTime; % acquisition time
        % convert the timestamp from HHMMSS.SSSSSS to "msec since midnight"
        hrsInMsec = 60*60*1000*str2double(dicom.(thisSeries).timeStr(1:2));
        minsInMsec = 60*1000*str2double(dicom.(thisSeries).timeStr(3:4));
        secInMsec = 1000*str2double(dicom.(thisSeries).timeStr(5:end));
        dicom.(thisSeries).AT(t,1) = hrsInMsec + minsInMsec + secInMsec;
    end
    % verify that timestamps are non-decreasing across images
    % (i.e., that dicoms were loaded in the correct order)
    if ~all(diff(dicom.(thisSeries).AT)>=0)
        dicom.(thisSeries).AT_badorder = dicom.(thisSeries).AT;
        [dicom.(thisSeries).AT,dicom.(thisSeries).AT_trueorder] = sort(dicom.(thisSeries).AT);
    end
    % Bug with the first dicom AT
    if nTRs>1
        dicom.(thisSeries).AT(1) = dicom.(thisSeries).AT(2) - dicom.(thisSeries).TR_all(2);
    end
end

%% Get the start and stop times from the .puls
pulsFiles = dir(fullfile(pulsDir,'*.puls'));
pulsFiles = {pulsFiles(:).name}';
isHidden = strncmp('.',pulsFiles,1); % remove hidden files from the list
pulsFiles(isHidden) = [];
numPulsFiles = length(pulsFiles);
pulsDates = nan(numPulsFiles,1);
pulsTimes = cell(numPulsFiles,1);
display('Reading pulse files');
for f=1:numPulsFiles
    pulsDates(f) = str2double(pulsFiles{f}(12:19));  % Date in double format
    % load the .puls file
    fid = fopen(fullfile(pulsDir,pulsFiles{f}));
	%display(['     file ' num2str(f) ': ' fullfile(pulsDir,pulsFiles{f})]);
    % For the old format of .puls files
    if OLD
        textscan(fid,'%s',4); % Ignore first 4 values. % This was for the old, manual way of recording/saving .puls files
    else
        textscan(fid,'%s',19); % Ignore first 19 values, up to "PULS_SAMPLE_INTERVAL = 20000" % This is for the new, automated way of recording/saving .puls files
    end
    data = textscan(fid,'%u16'); % Read data until end of u16 data.
    footer = textscan(fid,'%s');   % Read in remaining footer info
    fclose(fid);
    % Get timing information
    footStampLabels = {'LogStartMDHTime', 'LogStopMDHTime'};
    for i = 1:length(footStampLabels)
        labelStr = footStampLabels{i};
        labelIdx = find(strcmp(footer{1},[labelStr,':']))+1; % position of this timestamp in 'footer'
        pulse.(labelStr) = str2double(footer{1}{labelIdx}); % store the timestamp in a more convenient data structure
    end
    pulse.all = double(data{1}); % read in all values from data vector in .puls file
    pulse.Idx = find(pulse.all<5000);% % Index of actual data values (i.e. not marker values >5000).
    pulse.data = pulse.all(pulse.Idx); % actual data values
    pulse.nSamps = length(pulse.Idx); % number of pulse values
    pulse.dur = pulse.LogStopMDHTime - pulse.LogStartMDHTime; % Duration of pulse data (msec)
    pulse.sampR = pulse.dur/(pulse.nSamps-1); % sampling rate in msec (~20 msec)
    pulse.Hz = round(1000/pulse.sampR); % Compute Hz from sampling rate (in msec) (round for PulseRegress later)
    assert(and(floor(pulse.sampR)>=19,ceil(pulse.sampR)<=22),'pulse sampling frequency is not 20 msec');
    % Determine time stamp (tstmp) for each Pulse value
    for i = 1:pulse.nSamps
        pulse.AT(i) = pulse.LogStartMDHTime + pulse.sampR*(i-1);
    end
    pulsTimes{f} = pulse.AT; % Time when pulse signal was received (msec since midnight)
    clear pulse
end

%% Get summary info for each series
seriesName = cell(numSeries,1);
seriesDates = nan(numSeries,1);
seriesTimes = nan(numSeries,1);
for s=1:numSeries
    seriesName{s} = series.List(s).name;
    seriesName{s}(seriesName{s}=='.') = '_';
    seriesDates(s) = str2double(dicom.(seriesName{s}).info.AcquisitionDate); % Date in double format
    seriesTimes(s) = dicom.(seriesName{s}).AT(1); % Time when first TR was received (msec since midnight)
end

%% Find the most likely puls file for each series
SCANxPULS_match = cell(numSeries,1);
timeDiff_match = cell(numSeries,1);
for s=1:numSeries
    try
        scanDate = str2double(dicom.(seriesName{s}).info.AcquisitionDate);
        scanStart = dicom.(seriesName{s}).AT(1);
        scanEnd = dicom.(seriesName{s}).AT(end);
        dateMatch = find(pulsDates==scanDate);
        % Loop through puls files searching for best match
        for f=1:length(dateMatch)
            fileIdx = dateMatch(f);
            pulsStart = min(pulsTimes{fileIdx});
            pulsEnd = max(pulsTimes{fileIdx});
            % Check if windows match
            [~,I] = min(abs(scanStart-pulsTimes{fileIdx})); % Pulse signal in this file that is closest to first TR
            timeDiff = scanStart-pulsTimes{fileIdx}(I);
            % Check conditions
            cond1 = and(  pulsStart<=scanStart  ,  pulsEnd>=scanEnd  ); % Pulse signals must have started before the start of the scan and ended after the end of the scan
            cond2 = abs(timeDiff) < 10000; % Some pulse signal must have been acquired within 10 seconds of the first TR
            if and(cond1,cond2)
                SCANxPULS_match{s} = [SCANxPULS_match{s},fileIdx];
                timeDiff_match{s} = [timeDiff_match{s},timeDiff];
            end
        end
    catch
    end
end

%% Remove duplicates
matchFound = find(~cellfun('isempty',SCANxPULS_match));
% Remove duplicates for a given scan
for i=1:length(matchFound)
    s1=matchFound(i); s1X = SCANxPULS_match{s1};
    if length(s1X)>1
        [~,bestmatch] = min(abs(timeDiff_match{s1}));
        SCANxPULS_match{s1} = SCANxPULS_match{s1}(bestmatch);
        timeDiff_match{s1} = timeDiff_match{s1}(bestmatch);
    end
end
% Remove duplicates between scans
for i=1:length(matchFound)
    s1=matchFound(i);
    for j=1:length(matchFound)
        s2=matchFound(j);
        if and(s1~=s2,any(SCANxPULS_match{s1}==SCANxPULS_match{s2}))
            if abs(timeDiff_match{s1})<abs(timeDiff_match{s2})
                SCANxPULS_match{s2} = 1000+rand;
                timeDiff_match{s2} = 1000+rand;
            else
                SCANxPULS_match{s1} = 1000+rand;
                timeDiff_match{s1} = 1000+rand;
            end
        end
    end
end

%% Find the most likely puls file for each series
seriesPulsfile = cell(size(seriesName));
minTimeDiff = nan(size(seriesName));
for s=1:numSeries
    if and(~isempty(SCANxPULS_match{s}),SCANxPULS_match{s}<1000)
        seriesPulsfile{s} = pulsFiles{SCANxPULS_match{s}};
        minTimeDiff(s) = timeDiff_match{s};
    end
end
