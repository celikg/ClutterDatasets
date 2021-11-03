close all;
clear all;

c = 3e8;                                                                    % Define constants: Speed of Light 
 
% Load one dataset
DataStrc = load('ProcRL_Urban_300m_400M_700m_900m.mat');
Data = DataStrc.RangeLines_AfterEq;                                         % Radar data
[N_pulses,N_bins] = size(Data); 

Data = Data - repmat(mean(Data),N_pulses,1);                                % Remove the DC component before performing any statistical analaysis 

Bandwidth_Hz = DataStrc.Info.Bandwidth_Hz;
PRF_Hz = DataStrc.Info.PRF_Hz;
RangeStart = DataStrc.Info.RangeStart_m; 


%% Plot data: Range-Doppler plot
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
RangeAxis_m = ((0:1:(N_bins-1)) * RangeResolution_m) + 4965;          % Range Axis
FrequencyAxis_m = (-N_pulses/2:1:(N_pulses/2 - 1)) * FrequencyResolution_m; % Frequency Axis

% imagesc(RangeAxis_m, FrequencyAxis_m, RangeDoppler_dB,clims);               % Plot Range-Doppler map
% xlabel('Range (km)');
% ylabel('Doppler Frequency [Hz]');
% ylim([-200 200])
% grid on;
% colorbar;     
% colormap('jet');
% ax = gca;
% ax.XTick = 4000:100:6000;
% ax.YDir = 'normal';
% tick_scale_factor = 0.001;
% ax.XTickLabel = ax.XTick * tick_scale_factor;
% set(0,'defaultfigurecolor',[1 1 1])
% set(gca,'FontSize',14)
% set(gca,'Xcolor', 'k','YColor','k')


%% Identify Range Bins 
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


%% Plot PSD of Data #1

newRes = PRF_Hz/4096;                                                       % Frequency Resolution = 0.98 or 0.49

% Find range bins for Bandwidth 0.5Hz, 1Hz, 2Hz, 10Hz

BinHalfHz = round(0.5/newRes);  
BinOneHz = round(1/newRes);
BinTwoHz = round(2/newRes);
BinTenHz = round(10/newRes);

% % Alter the Row of Data by Zero Padding
% 
% PSD_Matrix1 = fft(Data(BinHalfHz,:),4096); 
% PSD_Matrix2 = fft(Data(BinOneHz,:),4096); 
% PSD_Matrix3 = fft(Data(BinTwoHz,:),4096); 
% PSD_Matrix4 = fft(Data(BinTenHz,:),4096); 
% 
% PSD = [PSD_Matrix1 
%        PSD_Matrix2 
%        PSD_Matrix3 
%        PSD_Matrix4];
%     
% PSD_normalised = abs(PSD)/max(max(abs(PSD)));                               % Normalising so the peak value is 1 (linear) or 0dB
% PSD_shifted = fftshift(PSD_normalised,1);                                   % Shift Zero Frequency Comp to centre
% PSD_dB = 20*log10(PSD_shifted);
% 
% 
% % FrequencyAxis_m = (-4096/2:1:(4096/2 - 1)) * FrequencyResolution_m*0.01;    % Frequency Axis
% % plot(FrequencyAxis_m,PSD_dB(1,:))
% 
% plot(PSD_dB(4,:))
% 


[row, col] = find(0.5 > RangeDoppler_normalised > 0);                                    

RangeBin = unique(col,'rows'); 


% Setup parameters for the computation of the Power Spectral Density 
WindowSize = 256;     % 256 FFT window
Nonoverlap = 32;     % 128 samples of overlap from section to section in FFT
Nfft = 4096;          % Number of FFT. Amer used 1024 

Data_PSD = zeros(Nfft,N_bins);


for RangeBinNum = 1:N_bins
       
         % periodograms of all range cells
         Data_PSD(:,RangeBinNum) = fftshift(pwelch(Data(:,RangeBinNum),WindowSize,Nonoverlap,Nfft)) * PRF_Hz/2;

end


% New Frequency Resolution
newRes = PRF_Hz/Nfft;


% Frequency Axis
FreqAxis = ((-Nfft/2):1:(Nfft/2-1)) * newRes;


% Plot the PSD for a single range bin 
% RangeBinToPlot = 950;
% plot(FreqAxis,10*log10(abs(Data_PSD(:,RangeBinToPlot))),'LineWidth',2);   
% hold on;

% RangeBinToPlot = 2310;
% plot(FreqAxis,10*log10(abs(Data_PSD(:,RangeBinToPlot) )),'-*','LineWidth',2);   
% hold on;
%   
% RangeBinToPlot = 1990;
% plot(FreqAxis,10*log10(abs(Data_PSD(:,RangeBinToPlot) )),'-*','LineWidth',2);   
% hold on;

% RangeBinToPlot = 2000;
% plot(FreqAxis,10*log10(abs(Data_PSD(:,RangeBinToPlot) )),'LineWidth',2);   
% hold on;

% ylim([-35 0]);
% xlim([-50 50]);
% 
% title('Power Spectral Density')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (dB)')
% set(0,'defaultfigurecolor',[1 1 1])

%% Plot Histogram of Data #1
% 
% figure;
% [rownum,colnum] = size(NewDataCol);
% 
% histfit(NewDataCol, rownum)
% xlabel('Amplitude')
% ylabel('Frequency')
% 
% set(gca,'FontSize',12)
% set(gca,'Xcolor', 'k','YColor','k')
% grid on;
% 
% xlim([0 0.15])
% ylim([0 130])

%% Plot Amplitude PDF #1

% figure;
% ksdensity(NewDataCol,'NumPoints',800);
% xlabel('Amplitude')
% ylabel('PDF')
% title('Amplitude PDF of Urban Vegetation Clutter')
% grid on;
% 
% set(gca, 'YScale', 'log')
% 
% xlim([0 0.15])
% ylim([0.0001 100])
% set(gca,'FontSize',12)
% set(gca,'Xcolor', 'k','YColor','k')
% 
% legend('400 MHz')


%% Plot Amplitude PDF #2

% 
% BinSize = 0.004;                                                          % Variable Bin Size
% BinAxis = 0:BinSize:max(NewDataCol);
% PDF = hist(NewDataCol,BinAxis)./(rownum*BinSize);
% 
% AreaUnderPDF = sum(PDF*BinSize);                                          % Ensure Area of PDF = 1 
% figure; semilogy(BinAxis,PDF,'*b')                                        % Plot Amplitude PDF
% xlabel('Ampitude');
% ylabel('PDF');
% grid on;
% 
% xlim([0 0.15]);
% ylim([1e-4 100]);
% 
% % Plot Closed form expression of Gaussian PDF
% norm = normpdf(BinAxis,0.0152434, 0.0164434);
% hold on;
% plot(BinAxis,norm,'c','LineWidth',1);
% 
% % Plot Closed form expression of Rayleigh PDF
% rayl = raylpdf(BinAxis,0.0158548);
% hold on;
% plot(BinAxis,rayl,'LineWidth',1);
% 
% % Plot Closed form expression of Weibull PDF
% weibull = wblpdf(BinAxis,0.0164112,1.23535);
% hold on;
% plot(BinAxis,weibull,'LineWidth',1);
% 
% % Plot Closed form expression of Lognormal PDF
% LogPDF = makedist('Lognormal','mu',-4.50971,'sigma',1.59461 );
% logn = pdf(LogPDF,BinAxis);
% hold on;
% plot(BinAxis,logn,'r','LineWidth',1);
% 
% 
% % Plot Closed form expression of Exponential PDF
% exp = exppdf(BinAxis,0.0152434);
% hold on;
% plot(BinAxis,exp,'LineWidth',1);
% 
% legend('Land Clutter Data', 'Normal Distribution','Rayleigh Distribution','Weibull Distribution', 'Lognormal Distribution', 'Exponential Distribution');
% set(gca,'FontSize',12)
% set(gca,'Xcolor', 'k','YColor','k')
% grid on;
% set(0,'defaultfigurecolor',[1 1 1])


%% Fit Distributions to Amplitude PDF
% %   Number of datasets:  1
% %   Number of fits:  5
% 
% %   Data from dataset "Land Clutter":
% %   Y = NewDataCol
% 
% % Ensure  all inputs to be column vectors
% NewDataCol = NewDataCol(:);
% 
% % Prepare figure
% clf;
% hold on;
% LegHandles = []; LegText = {};
% 
% % --- Plot data originally in dataset "Land Clutter"
% [CdfF,CdfX] = ecdf(NewDataCol,'Function','cdf');                            % compute empirical cdf
% BinInfo.rule = 5;
% BinInfo.width = 0.01;
% BinInfo.placementRule = 1;
% [~,BinEdge] = internal.stats.histbins(NewDataCol,[],[],BinInfo,CdfF,CdfX);
% [BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
% hLine = bar(BinCenter,BinHeight,'hist');
% set(hLine,'Marker','x','MarkerSize',10,'LineWidth',0.1);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Land Clutter';
% 
% % Create grid where function will be computed
% XLim = get(gca,'XLim');
% XLim = XLim + [-1 1] * 0.01 * diff(XLim);
% XGrid = linspace(XLim(1),XLim(2),1000);
% 
% % Create fit "Normal Distribution"
% % Fit this distribution to get parameter values
% % pd1 = ProbDistUnivParam('normal',[0.007769446497716, 0.006702028502241])
% 
% pd1 = fitdist(NewDataCol, 'normal');
% YPlot = pdf(pd1,XGrid);
% hLine = plot(XGrid,YPlot,'Color',[1 0 0],'LineWidth',2);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Normal Distribution';
% 
% % Create fit "Lognormal Distribution"
% % Fit this distribution to get parameter values
% % pd2 = ProbDistUnivParam('lognormal',[ -5.072674904059, 1.118590200882])
% 
% pd2 = fitdist(NewDataCol, 'lognormal');
% YPlot = pdf(pd2,XGrid);
% hLine = plot(XGrid,YPlot,'Color',[0 0 1],'LineWidth',2);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Lognormal Distribution';
% 
% % Create fit "Exponential Distribution"
% % Fit this distribution to get parameter values
% % pd3 = ProbDistUnivParam('exponential',[0.007769446497716])
% pd3 = fitdist(NewDataCol, 'exponential');
% YPlot = pdf(pd3,XGrid);
% hLine = plot(XGrid,YPlot,'Color',[0.6 0.3 1],'LineWidth',2);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Exponential Distribution';
% 
% % Create fit "Rayleigh Distribution"
% % Fit this distribution to get parameter values
% % pd4 = ProbDistUnivParam('rayleigh',[0.007255393006031])
% pd4 = fitdist(NewDataCol, 'rayleigh');
% YPlot = pdf(pd4,XGrid);
% hLine = plot(XGrid,YPlot,'LineWidth',2);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Rayleigh Distribution';
% 
% % Create fit "Weibull Distribution"
% % Fit this distribution to get parameter values
% % pd5 = ProbDistUnivParam('weibull',[ 0.008580929344244, 1.475081133578])
% pd5 = fitdist(NewDataCol, 'weibull');
% YPlot = pdf(pd5,XGrid);
% hLine = plot(XGrid,YPlot,'LineWidth',2);
% LegHandles(end+1) = hLine;
% LegText{end+1} = 'Weibull Distribution';
% 
% % Create legend from accumulated labels, Intervals and Axes
% hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 8);
% set(hLegend);
% 
% xlabel('Amplitude');
% ylabel('PDF');
% set(gca,'YScale','log');
% 
% xlim([0 0.15]);
% ylim([1e-5 100]);
% 
% set(gca,'FontSize',12);
% set(gca,'Xcolor', 'k','YColor','k');
% grid on;