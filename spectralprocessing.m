%
% Script file to process the dynamics data from Northlan Power 
% to get spectrum and frequency shifting information
%

% Specify the processing parameters
fs = 24000;          %Sample rate in Hz
nfft = 16384;        %FFT length
scanToAvg = 10;      %Number of scans to average
windowParam = hanning(nfft);        % Windowing function
windowScaling = sum(windowParam);   % Scaling factor since we are windowing
noverlap = 0.75 * nfft;              % Overlap
fsVector = (0:nfft/2)*fs/nfft;      % Frequency vector
load canmapinfo2.mat;

% ToneBound=[980 1020;
%            600 650;
%            410 450;];
ToneBound=[980 1020;
           550 700;
           400 500;];
HighToneIndex = find(fsVector>=ToneBound(1,1) & fsVector<=ToneBound(1,2));
MidToneIndex = find(fsVector>=ToneBound(2,1) & fsVector<=ToneBound(2,2));
LowToneIndex = find(fsVector>=ToneBound(3,1) & fsVector<=ToneBound(3,2));

load id5LRch1.mat;
dynSignal = (id5LRch1/0.1)*139.167;  %Conversion rate is 100mv/psi, scaling factor is 139.167
dynSignal = dynSignal(:);
clear id5LRch1
    
%Specify the length of the data signal
ndata = length(dynSignal);
	
%Calculate the spectrogram information
numWindows = fix((ndata - noverlap) / (nfft - noverlap));
timeBias = 0;
fftInfoIndex = 1;
scanIndex = 1;
index = (1:nfft)';
	
% Now go through the FFT calculations
for j=1:numWindows
    tValue = max(index) / fs + timeBias;
    dataToFFT = dynSignal(index) .* windowParam;
    index = index + (nfft - noverlap);
    fftData = abs(fft(dataToFFT,nfft)) / windowScaling;
    if scanIndex ~= scanToAvg
        scanAvg(:,scanIndex) = fftData;
        scanIndex = scanIndex + 1;
    else
        scanAvg(:,scanIndex) = fftData;
        if scanToAvg ~= 1
            finalFFTData(:,fftInfoIndex) = mean(scanAvg')';
        else
            finalFFTData(:,fftInfoIndex) = scanAvg;
        end
        HighToneFFTData = finalFFTData(HighToneIndex,fftInfoIndex);
        HighToneDFVec(fftInfoIndex) = sum(fsVector(HighToneIndex)'.*(HighToneFFTData.^2)) / sum(HighToneFFTData.^2);
        MidToneFFTData = finalFFTData(MidToneIndex,fftInfoIndex);
        MidToneDFVec(fftInfoIndex) = sum(fsVector(MidToneIndex)'.*(MidToneFFTData.^2)) / sum(MidToneFFTData.^2);
        LowToneFFTData = finalFFTData(LowToneIndex,fftInfoIndex);
        LowToneDFVec(fftInfoIndex) = sum(fsVector(LowToneIndex)'.*(LowToneFFTData.^2)) / sum(LowToneFFTData.^2);
        timeVector(fftInfoIndex) = tValue;
        fftInfoIndex = fftInfoIndex + 1;
        scanIndex = 1;
    end
end
fileResults.timeVector = timeVector;
fileResults.fftData = 20*log10(2*finalFFTData(1:nfft/2+1,:));
fileResults.fsVector = fsVector;
freqShiftResults.HighTone = HighToneDFVec;
freqShiftResults.MidTone = MidToneDFVec;
freqShiftResults.LowTone = LowToneDFVec;
freqShiftResults.ToneBound=freqShiftResults;
save 'd:\liao\data\id2LRch1fft_1610.mat' fileResults freqShiftResults;
%save id10LRch1fft.mat fileResults freqShiftResults;

figure;
imagesc(timeVector,fsVector,fileResults.fftData);axis xy; colormap(map); %colorbar;    
xlabel('Time - sec');ylabel('Frequency - Hz');title('Spectrogram of id2LRch1');
axis([-inf inf 200 3000]);

clear all;

