function out = od_doprocess(pars,seismicin) 
    
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


    sampint = pars.Sampint;
    timegap= pars.MinTimeGap;
    ws = pars.Windowsize;
    
    phasecomp = pars.PhaseComp;
    combinationflag=pars.CombinationFlag;
    
    norm_flag = pars.NormalizationFlag;
    
    inline2Dsmooth_flag=pars.InLine2DSmoothingFlag;
    inline2Dsmoothingsize = pars.InLine2DSmoothingSize;
    xline2Dsmoothing_flag = pars.XLine2DSmoothingFlag;
    xline2Dsmoothingsize = pars.XLine2DSmoothingSize;
    
    
    % Input Format of 3D Seismic from Opendtect
    % Dimension: (TWT,Xline,Inline)
    seis3d_in1 = cell2mat( seismicin(1) );
    
    % Call function to carry out the phase decomposition (PD)
    % Method: Odd-Even Function
    % Phase Component Call valid only for -90, 0, +90, +/-180 degree 
    out = pdoe3dfilter2d( seis3d_in1, sampint, timegap,ws, ...
        phasecomp, combinationflag,...
        inline2Dsmooth_flag,inline2Dsmoothingsize, ...
        xline2Dsmoothing_flag, xline2Dsmoothingsize,norm_flag );
 
end