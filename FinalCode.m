%% Author: C. GULER
%% Adapted from Dr Abdul Gaffar
%% Date:   November 2021

% Code Written to Plot the range-Doppler Map, Model the Amplitude Distribution 
% and Generate the PSD of Land Clutter  

close all;
clear all;
set(0,'defaultfigurecolor',[1 1 1])

c = 3e8;                                                                    % Define constants: Speed of Light 
 
Load one dataset
DataStrc = load('ProcRL_south_east_into_ridge_25Hz_3500m.mat');
Data = DataStrc.RangeLines_AfterEq;                                         % Radar data
[N_pulses,N_bins] = size(Data); 

Data = Data - repmat(mean(Data),N_pulses,1);                                % Remove the DC component before performing any statistical analaysis 

Bandwidth_Hz = DataStrc.Info.Bandwidth_Hz;
PRF_Hz = DataStrc.Info.PRF_Hz;
RangeStart = DataStrc.Info.RangeStart_m; 


%% 1. Plot data: Range-Doppler Map
NumPulsesToProcess = N_pulses; 
RangeDoppler = fft(Data(1:NumPulsesToProcess,:), [], 1);                    % Only applying the FFT in the slow-time dimension 
                                                                            % Data(1:NumPulsesToProcess,:) takes all the data from row 1 to row 'N_pulses'        
clims = [-45 0];                                                            % Only plot values from -45dB to 0dB
RangeDoppler_normalised = abs(RangeDoppler)/max(max(abs(RangeDoppler)));    % Normalising so the peak value is 1 (linear) or 0dB
RangeDoppler_shifted = fftshift(RangeDoppler_normalised,1);                 % Shift Zero Frequency Comp to centre
RangeDoppler_dB = 20*log10(RangeDoppler_shifted);


%% Resolution Definitions
RangeResolution_m = c/(2 * Bandwidth_Hz);                                   % Range Resolution = 0.375
FrequencyResolution_m = PRF_Hz/N_pulses;                                    % Frequency Resolution = 1.5898
RangeAxis_m = ((0:1:(N_bins-1)) * RangeResolution_m);          % Range Axis
FrequencyAxis_m = (-N_pulses/2:1:(N_pulses/2 - 1)) * FrequencyResolution_m; % Frequency Axis

imagesc(RangeAxis_m, FrequencyAxis_m, RangeDoppler_dB,clims);               % Plot Range-Doppler map
xlabel('Range [km]');
ylabel('Doppler Frequency [Hz]');
ylim([-200 200])
grid on;
colorbar;     
colormap('jet');

ax = gca;
ax.XTick = 1500:500:7000;
ax.YDir = 'normal';

tick_scale_factor = 0.001;
ax.XTickLabel = ax.XTick * tick_scale_factor;

set(gca,'FontSize',14)
set(gca,'Xcolor', 'k','YColor','k')


% Identify Range Bins 
NoiseBin = round(30/FrequencyResolution_m);                                 % Find Doppler bin at 30Hz
Noise = RangeDoppler_shifted(NoiseBin:end,:);
Noiseindex = Noise == -Inf;                                                 % Remove any instances of infinite values 
Mean_Noise = mean(Noise(Noiseindex==0),'all');                              % Noise Power Average past the 30Hz Doppler bin

NewThreshold = Mean_Noise * 3;                                              % Land Clutter Power Threshold
[row, col] = find(RangeDoppler_shifted > NewThreshold);                                    

Bin_m = unique(col,'rows');                                                 % Find Unique Cols possessing Land Clutter


%% Remove Noise Data from Data Matrix  
New_Data = RangeDoppler_shifted(:,Bin_m);


%% Convert Data to Column Vector  
NewDataCol = New_Data(:);


%% Plot Histogram of Data
figure;
[rownum,colnum] = size(NewDataCol);

histfit(NewDataCol, rownum)
xlabel('Amplitude')
ylabel('Frequency')

set(gca,'FontSize',12)
set(gca,'Xcolor', 'k','YColor','k')
grid on;

xlim([0 0.15])
ylim([0 130])


%% 2. Plot Amplitude PDF 
figure;
ksdensity(NewDataCol,'NumPoints',800);
xlabel('Amplitude')
ylabel('PDF')
grid on;

set(gca, 'YScale', 'log')

xlim([0 0.5])
ylim([0.0001 100])
set(gca,'FontSize',12)
set(gca,'Xcolor', 'k','YColor','k')

legend('10 MHz')


%% Plot Amplitude PDF with Parameter Values

BinSize = 0.004;                                                          % Variable Bin Size
BinAxis = 0:BinSize:max(NewDataCol);
PDF = hist(NewDataCol,BinAxis)./(rownum*BinSize);

AreaUnderPDF = sum(PDF*BinSize);                                          % Ensure Area of PDF = 1 
figure; semilogy(BinAxis,PDF,'*b')                                        % Plot Amplitude PDF
xlabel('Ampitude');
ylabel('PDF');
grid on;

xlim([0 0.1]);
ylim([1e-4 100]);

% Plot Closed form expression of Gaussian PDF
norm = normpdf(BinAxis,0.0152, 0.0164);
hold on;
plot(BinAxis,norm,'c','LineWidth',1);

% Plot Closed form expression of Rayleigh PDF
rayl = raylpdf(BinAxis,0.0159);
hold on;
plot(BinAxis,rayl,'LineWidth',1);

% Plot Closed form expression of Weibull PDF
weibull = wblpdf(BinAxis,0.0164,1.24);
hold on;
plot(BinAxis,weibull,'LineWidth',1);

% Plot Closed form expression of Lognormal PDF
LogPDF = makedist('Lognormal','mu',-4.51,'sigma',1.60);
logn = pdf(LogPDF,BinAxis);
hold on;
plot(BinAxis,logn,'r','LineWidth',1);

% Plot Closed form expression of Exponential PDF
exp = exppdf(BinAxis,0.0152);
hold on;
plot(BinAxis,exp,'LineWidth',1);

legend('Land Clutter Data', 'Normal Distribution','Rayleigh Distribution','Weibull Distribution', 'Lognormal Distribution', 'Exponential Distribution');
set(gca,'FontSize',12)
set(gca,'Xcolor', 'k','YColor','k')
grid on;
set(0,'defaultfigurecolor',[1 1 1])


%% 3. Generate PSD of Data

DataStrc = load('ProcRL_south_east_25Hz_3500m.mat');                        % load dataset 
Data = DataStrc.RangeLines_AfterEq;                                         % Radar data
[N_pulses,No_bins] = size(Data);
Data = Data - repmat(mean(Data),N_pulses,1); 
PRF_Hz = DataStrc.Info.PRF_Hz;                                              % Sampling rate of clutter [Hz]

% Setup parameters for the computation of the Power Spectral Density 
WinSize = 256; Nonoverlap = 32; Nfft = 1024;                                                                                                            
Data_PSD = zeros(Nfft,No_bins);

% Create Periodograms for all range cells 
for RangeBinNum = 1:No_bins
Data_PSD(:,RangeBinNum) = fftshift(pwelch(Data(:,RangeBinNum),WinSize,Nonoverlap,Nfft)) * PRF_Hz/2;
end

% New Frequency Bin Size
FreqBinSize = PRF_Hz/Nfft;

% Frequency Axis
FreqAxis = ((-Nfft/2):1:(Nfft/2-1)) * FreqBinSize;
PSD = 10*log10(abs(Data_PSD));   

% Plot the PSD for all range bins
Sum = PSD(:,1);

for Range_Bin = 2:No_bins
        Sum = Sum + PSD(:,Range_Bin); 
end


%% Plot Average PSD of All Data

Mean_Graph = Sum / No_bins;

Norm_Graph = Mean_Graph - max(Mean_Graph);
plot(FreqAxis,Norm_Graph)

ylim([-9 9]);
xlim([-30 30]);

xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')

legend('Past Southeast Ridge');
set(gca,'FontSize',12);
set(gca,'Xcolor', 'k','YColor','k');
set(0,'defaultfigurecolor',[1 1 1])
grid on;

ax = gca;
ax.YTick = -9:3:9;



%% 4. Categorise Range Bins according to -3dB Bandwidth

DataStrc = load('ProcRL_Urban_300m_400M_700m_900m.mat');            % load dataset 
Data = DataStrc.RangeLines_AfterEq;                                         % Radar data
[N_pulses,No_bins] = size(Data);
Data = Data - repmat(mean(Data),N_pulses,1); 
PRF_Hz = DataStrc.Info.PRF_Hz;                                              % Sampling rate of clutter [Hz]

% Setup parameters for the computation of the Power Spectral Density 
WinSize = 64; Nonoverlap = 2; Nfft = 1024;                                                                                                            
Data_PSD = zeros(Nfft,No_bins);

% Create Periodograms for all range cells 
for RangeBinNum = 1:No_bins
Data_PSD(:,RangeBinNum) = fftshift(pwelch(Data(:,RangeBinNum),WinSize,Nonoverlap,Nfft)) * PRF_Hz/2;
end

% New Frequency Bin Size
FreqBinSize = PRF_Hz/Nfft;

% Frequency Axis
FreqAxis = ((-Nfft/2):1:(Nfft/2-1)) * FreqBinSize;
PSD = 10*log10(abs(Data_PSD));   

AllBWs = [];

for Range_Bin = 1:No_bins

    x = PSD(:,Range_Bin);
    Max = max(x);
    thresh = Max-3;
    [rows, cols] = find(x > thresh);
    Bandwidth_Hz = numel(rows) * FreqBinSize;
    
    AllBWs = [AllBWs; Bandwidth_Hz];  
    
end

plot(FreqAxis,PSD(:,1));

BinsHalfHz = find(0.5 > AllBWs);
BinsOneHz = find(1 > AllBWs & AllBWs > 0.5);
BinsTwoHz = find(2 > AllBWs & AllBWs > 1);
BinsGreaterHz = find(AllBWs > 2);

plot(FreqAxis,PSD(:,100));
xlim([-40 40])
ylim([-56 -0])

