function norm_out=func_norm_mean_sd(data3d_in,data_flag_IN_XL,ns)
    %% Normalization of 3D data Cube w.r.t. Mean and Standard Deviation
    
    % data3d_in = 3D data (Seismic) Dimension: TWT x Xline x Inline
    
    %  Discarding the position values where data flag is zero
    % data_flag_IN_XL = 2D planner data flag Dimension: Xline x Inline
    % data_flag_IN_XL = 0 or 1 (0=Bad Data position, 1= Data Position to
    % - accept)
    
    % ns= number of samples along twt or 1st axis i.e. rows
    
    % Find the Index where data flag is not zero
    TF=find(data_flag_IN_XL~=0);

    % Normalization of the output w.r.t. mean & standard deviation
    % Mean of the 3D data volume
    mu_along_twt=mean(data3d_in);
    mu=mean(mu_along_twt(TF));

    % Standard Deviation of the 3D data volume
    rho_along_twt=std(data3d_in);
    rho=std(rho_along_twt(TF));

    % Normalizing the Volume w.r.t. mean & std
    norm_out=(data3d_in-mu)./rho;

    % Replace NaN element with zero value
    TF_nan=isnan(norm_out);
    norm_out(TF_nan)=0;

    % Data_flag=0 to be replaced as Zero Values
    for ii=1:ns
        norm_out(ii,:,:)=shiftdim(norm_out(ii,:,:)).*data_flag_IN_XL;
    end
end