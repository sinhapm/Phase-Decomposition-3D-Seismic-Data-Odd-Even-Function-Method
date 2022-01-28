% The Program performs the 3D Seismic Data Phase Decomposition based on 
% Odd-Even function. It computes only the odd components (-90 degree, 
% +90degree) and the even components (+/-180 degree, 0 degree)
% from the 3D seismic data

% The program reads the sample 3D seismic data set and the test the utility
% of the Phase Decomposition function "pdoe3dfilter2d()"

% Input 3D Seismic data in mat file , Dimension: (TWT,Xline,Inline)
% Sample Seismic data = "sample_3Dseismic.mat"

% Input Parameters
% Sampint= Sampling Interval (in msec) of Input 3D Seismic, sample seismic sampint=4msec 
% MinTimeGap = Minimum Time Gap (in msec) for two consecutive Local Maxima Points
% Windowsize = Time Window (in msec) to compute Frequency-Phase over Amplitude
% - Apex from the Odd-Even Segment trace data

% PhaseComp = Either -90 or 0 or 90 or 180 Degree Phase Components are the only valid input

% CombinationFlag = 0 for no combination or 1 for combination of such as
% e.g. if PhaseComp=-90 or +90, for CombinationFlag = 1, the single 
% output will have combine result of -90 and +90 degree, i.e. Odd Function
% e.g. if PhaseComp=0 or 180 or -180, for CombinationFlag = 1, the single 
% output will have combine result in even function phase component

% Normalization = 1 for yes or 0 for no, Seismic Volume Normalization by
% Mean & Standard Deviation Values

% InLine2DSmoothingFlag = Smoothing Filter Flag for Inline
% SmoothSizeILine = Smoothing Filter Size along InLine   

% XLine2DSmoothingFlag = Smoothing Filter Flag for CrossLine
% SmoothSizeXLine = Smoothing Filter Size along CrossLine

% Program tested on MATLAB (R2021)
% Additional Module requirement -> "Signal Processing"

% Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

% Clear the variable & terminal
clear;
clc;

%% Input Seismic & Parameters 
% Load the 3D Seismic data 
data_in=load('sample_3Dseismic.mat'); % Change the Path and name of Seismic data in mat format
seis3d_in1=data_in.seis3d_in1;
clear data_in;

% Parameters Input
sampint = 4; % Sampint
timeGap = 4; % MinTimeGap
ws = 16; % Windowsize
phasecomp=[-90 +90 0 180]; % PhaseComp
inline2Dsmooth_flag=1; % InLine2DSmoothingFlag
inline2Dsmoothingsize = 3; % SmoothSizeILine
xline2Dsmoothing_flag = 1; % XLine2DSmoothingFlag
xline2Dsmoothingsize = 3; % SmoothSizeXLine
norm_flag=1; % Normalization
combinationflag=0; % CombinationFlag

%% Call Phase Decomposition Function  

[out_phasedecom1] = pdoe3dfilter2d( seis3d_in1, sampint, timeGap,...
            ws, phasecomp(1),combinationflag,...
            inline2Dsmooth_flag, inline2Dsmoothingsize, ...
            xline2Dsmoothing_flag ,xline2Dsmoothingsize , norm_flag);

[out_phasedecom2] = pdoe3dfilter2d( seis3d_in1, sampint, timeGap,...
            ws, phasecomp(2),combinationflag,...
            inline2Dsmooth_flag, inline2Dsmoothingsize, ...
            xline2Dsmoothing_flag ,xline2Dsmoothingsize , norm_flag);

[out_phasedecom3] = pdoe3dfilter2d( seis3d_in1, sampint, timeGap,...
            ws, phasecomp(3),combinationflag,...
            inline2Dsmooth_flag, inline2Dsmoothingsize, ...
            xline2Dsmoothing_flag ,xline2Dsmoothingsize , norm_flag);

[out_phasedecom4] = pdoe3dfilter2d( seis3d_in1, sampint, timeGap,...
            ws, phasecomp(4),combinationflag,...
            inline2Dsmooth_flag, inline2Dsmoothingsize, ...
            xline2Dsmoothing_flag ,xline2Dsmoothingsize , norm_flag);

%% QC Plot Display
inlineno=10; % Enter the Inline Number to Display

seis3d_in1_2D=shiftdim(seis3d_in1(:,:,inlineno));
phase_decom1_2D=shiftdim(out_phasedecom1(:,:,inlineno));
phase_decom2_2D=shiftdim(out_phasedecom2(:,:,inlineno));
phase_decom3_2D=shiftdim(out_phasedecom3(:,:,inlineno));
phase_decom4_2D=shiftdim(out_phasedecom4(:,:,inlineno));

nfr=2; nfc=3;
figure;
subplot(nfr,nfc,1);imagesc(seis3d_in1_2D,"Interpolation","bilinear"); colormap gray;
title(strcat("Input Seismic - Inline# ",num2str(inlineno)));
ax1 = gca; axis tight; colorbar('eastoutside');

subplot(nfr,nfc,2);imagesc(phase_decom1_2D,"Interpolation","bilinear"); colormap gray;
title(strcat('Phase Decomposed:',num2str(phasecomp(1)),' Degree'));
ax2 = gca; axis tight;colorbar('eastoutside');

subplot(nfr,nfc,3);imagesc(phase_decom2_2D,"Interpolation","bilinear"); colormap gray;
title(strcat('Phase Decomposed:',num2str(phasecomp(2)),' Degree'));
ax3 = gca; axis tight;colorbar('eastoutside');

subplot(nfr,nfc,4);imagesc(seis3d_in1_2D,"Interpolation","bilinear"); colormap gray;
title(strcat("Input Seismic - Inline# ",num2str(inlineno)));
ax4 = gca; axis tight; colorbar('eastoutside');

subplot(nfr,nfc,5);imagesc(phase_decom3_2D,"Interpolation","bilinear"); colormap gray;
title(strcat('Phase Decomposed:',num2str(phasecomp(3)),' Degree'));
ax5 = gca; axis tight;colorbar('eastoutside');

subplot(nfr,nfc,6);imagesc(phase_decom4_2D,"Interpolation","bilinear"); colormap gray;
title(strcat('Phase Decomposed:',num2str(phasecomp(4)),' Degree'));
ax6 = gca; axis tight;colorbar('eastoutside');

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6]);
