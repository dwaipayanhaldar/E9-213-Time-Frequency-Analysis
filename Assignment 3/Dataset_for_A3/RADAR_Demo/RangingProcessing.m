function [] = RangingProcessing(matFile)
% Produces an RTI (range x time intensity) image of the
% radar recording. Also applies a simple two-pulse
% canceller to filter out clutter and CFAR normalization
% to improve visual detection. 
%
%    matFile = the filename of the .mat file to process
%
% MIT Build a Radar Course
% (c) 2019 Massachusetts Institute of Technology

% -------------------- Setup constants and parameters --------------------
c = 299792458;      % (m/s) speed of light
pw = 20e-3;         % (s) pulse length
fs = 40000;         % (Hz) sample rate
ovsRng = 10;         % oversampling factor applied in range
fStart = 2400e6;    % (Hz) LFM start frequency
fStop = 2480e6;     % (Hz) LFM stop frequency
nPulseCancel = 2;   % Number of pulses to use for canceller
maxRange = 50;     % (m) maximum range to display
% ------------------------------------------------------------------------

% Load example data
if(nargin==0)
    fprintf('Loading example data...');
    % matFile = 'RangeExample';
    matFile = 'Demo5';
end

% Calculate bandwidth and display
fprintf('Using %g MHz bandwidth\n', (fStop-fStart)*1e-6);
fprintf('Loading *.mat file...\n');

% Load mat file based on passed filename
a = load(matFile);
s = a.samples;

% Derived parameters
Np = round(pw * fs) - 1;    % of samples per pulse (-1 due to start flag)
bw = fStop - fStart;        % (Hz) transmit bandwidth

% Parse data and format as [pulse x sample]
fprintf('Parsing the recording...');
[parseData,numPulses] = parse_matFile(s,Np);
fprintf('Found %d pulses\n',numPulses);

% Taper applied to reduce range sidelobes
rngWin = repmat(hann_window(Np)',numPulses-1,1);

% Fourier transform to convert frequency to range
fprintf('Transforming from frequency to range...\n');
nfft = Np*ovsRng;
range = fftshift(fft(parseData.*rngWin,nfft,2),2);

% Compute the range axis
freqAxis = linspace(-fs/2,fs/2,nfft);
rangeAxis = c*pw*freqAxis/(2*bw);

% Compute the time axis
t = (0:numPulses-1)*pw*2;

% Trim data and range axis based on specified max range
samplesKeep = (rangeAxis<maxRange & rangeAxis>=0);
rangeAxisKeep = rangeAxis(samplesKeep);
range = range(:,samplesKeep);

% Apply N-pulse canceller
fprintf('Cancelling stationary background clutter...\n');
for i = 1:numPulses-nPulseCancel
    pulses = range(i:i+nPulseCancel-1,:);
    mask = -ones(nPulseCancel,length(rangeAxisKeep));
    mask(1,:) = 1;
    mask(2:end,:) = mask(2:end,:)./(nPulseCancel-1);
    rangeCancelled(i,:) = sum(pulses.*mask,1);
end

% Apply the median CFAR normalization
% rangeCancelled_dB = 20*log10(abs(rangeCancelled));
% rangeCancelled_dB = rangeCancelled_dB - repmat(median(rangeCancelled_dB,1),[size(rangeCancelled_dB,1) 1]); % over time
% rangeCancelled_dB = rangeCancelled_dB - repmat(median(rangeCancelled_dB,2),[1 size(rangeCancelled_dB,2)]); % over range

% ----------------------------- Plotting ---------------------------------
fprintf('Plotting...\n');

% Display RTI without clutter rejection.
figure, imagesc(rangeAxisKeep,t,20*log10(abs(range)));
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI without Clutter Rejection');
colormap(jet(256));
colorbar;

% Display the MTI results
figure, imagesc(rangeAxisKeep,t,20*log10(abs(rangeCancelled)));
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI with MTI Clutter Rejection');
colormap(jet(256));
colorbar;

% Display the CFAR results
% figure, imagesc(rangeAxisKeep,t,rangeCancelled_dB);
% ylabel('Time (s)');
% xlabel('Range (m)');
% title('RTI with MTI Clutter Rejection and CFAR');
% colormap(jet(256)); caxis([0 35]);
% colorbar;


% ---------------------------- Functions ---------------------------------

% Create a hann (cosine squared) window
function [w] = hann_window(N)
w = .5 + .5*cos(2*pi*((1:N).'/(N+1) - .5));

% Parses the digital data by looking for "start of pulse" flags. Converts
% data to voltages
function [parseData,numPulses] = parse_matFile(x,samplesPerPulse)
flag_pulseStart = 5000;
scale_factor = 3.3/2^12;
pulseStarts = find(x>flag_pulseStart);
numPulses = sum(x(:)>flag_pulseStart);
parseData = zeros(numPulses-1,samplesPerPulse); 
ctr = 1;
for i = 1:length(pulseStarts)-1
    parseData(ctr,:) = x(pulseStarts(i)+1:pulseStarts(i)+samplesPerPulse);
    ctr = ctr + 1;
end
parseData=(parseData*scale_factor)-3.3/2; % convert digital data to voltages (3.3V/12-bits) centered about 0V



