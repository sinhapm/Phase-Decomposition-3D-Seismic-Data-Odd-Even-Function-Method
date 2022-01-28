function pars = od_getparameters() 

% Function call to initialize the opendtect get parameters

% Both function files od_doprocess() and od_getparametrs() are the initial 
% requisite function files to get the input seismic data (single 3D seismic
% volume only) and other input parameters via method 
% "OpendTect"s main menu > Analysis > Volume Builder" step of Opendtect
% Files are required to be converted to shared c library dll files through
% Matlab compiler & Compiler SDK (c shared library dll files)
% For further check the reference link as mentioned along with attached pdf
% for Matlab Environment link setup
% https://doc.opendtect.org/6.0.0/doc/od_userdoc/Default.htm#appendix_f-matlab_link_plugin_for_opendtect/usage_through_opendtects_gui/accessing_compiled_matlab_functions_in_opendtects_gui.htm%3FTocPath%3D16%2520Appendix%2520F%2520-%2520MATLAB%2520Link%2520Plugin%2520for%2520OpendTect%2520v6.0%7C16.3%2520Option%25201%253A%2520Usage%2520Through%2520OpendTect's%2520GUI%7C_____5
% https://doc.opendtect.org/6.0.0/doc/od_userdoc/Default.htm#appendix_f-matlab_link_plugin_for_opendtect/usage_through_opendtects_gui/compilation_of_matlab_functions.htm%3FTocPath%3D16%2520Appendix%2520F%2520-%2520MATLAB%2520Link%2520Plugin%2520for%2520OpendTect%2520v6.0%7C16.3%2520Option%25201%253A%2520Usage%2520Through%2520OpendTect's%2520GUI%7C_____4
% https://doc.opendtect.org/6.0.0/doc/od_userdoc/Default.htm#appendix_f-matlab_link_plugin_for_opendtect/background.htm%3FTocPath%3D16%2520Appendix%2520F%2520-%2520MATLAB%2520Link%2520Plugin%2520for%2520OpendTect%2520v6.0%7C_____1

% Opendtect tested version: 6.4
% https://dgbes.com/index.php/download

% Matlab Runtime 9.8 or later, to run the program without MATLAB License 
% https://in.mathworks.com/products/compiler/matlab-runtime.html

% Input Parameters details

   % Sampint= Sampling Interval (in msec) of Input 3D Seismic 
   % MinTimeGap = Minimum Time Gap (in msec) for two consecutive Local Maxima Points
   % Windowsize = Time Window (in msec) to compute Frequency-Phase over Amplitude
   % - Apex from the Odd-Even Segment trace data
   
   % PhaseComp= Either -90 or 0 or 90 or 180 Degree Phase Components are the only valid input
   
   % CombinationFlag = 0 for no combination or 1 for combination of such as
   % e.g. if PhaseComp=-90 or +90, for CombinationFlag = 1, the single 
   % output will have combine result of -90 and +90 degree, i.e. Odd
   % Function
   % e.g. if PhaseComp=0 or 180 or -180, for CombinationFlag = 1, the single 
   % output will have combine result in even function phase component
   
   % Normalization = 1 for yes or 0 for no
   
   % InLine2DSmoothingFlag = Smoothing Filter Flag for Inline
   % SmoothSizeILine = Smoothing Filter Size along InLine   
   
   % XLine2DSmoothingFlag = Smoothing Filter Flag for CrossLine
   % SmoothSizeXLine = Smoothing Filter Size along CrossLine

   % Opendtect Version 6.4, Uncheck Batch Processing while running the volume builder
   % Single Seismic 3D Volume as input
   % Minimum Matlab Runtime Version 9.8
   
   % Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

    pars.nrinputs = 1; 
    
    pars.Sampint = 4;
    pars.MinTimeGap = 4; 
    pars.Windowsize = 16;
    
    pars.PhaseComp=-90; % Only valid inputs are -90, 0, 90, 180
    pars.CombinationFlag=0;
    
    pars.NormalizationFlag=1;
    
    pars.InLine2DSmoothingFlag=1;
    pars.InLine2DSmoothingSize=5;
    pars.XLine2DSmoothingFlag=1;
    pars.XLine2DSmoothingSize=5;

end