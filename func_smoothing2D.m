function out_smooth=func_smoothing2D(data2d_in,filt_window)

% Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

    % data3d_in = 2D Seismic Data Traces Input Dimension: TWT x Traces
    % filt_window = Filter Window Size such 3X3 or 5X5 or 7X7 etc
        
    %% 2D Smoothing Filter from filter length such as 3X3X3, 5X5X5 etc 
   
        out_smooth = filter2(filt_window, data2d_in);
    
end