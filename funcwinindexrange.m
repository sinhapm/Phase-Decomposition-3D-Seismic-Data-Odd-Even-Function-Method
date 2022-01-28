function [indx_range]=funcwinindexrange(inputdata,ws,dt,N_2)

% Function call to find the index range of window size
% Author: Prashant M. Sinha, sinha.pm@gmail.com, 2021

    [~,apos_max_indx]=max(abs(inputdata));
    indx_range=(apos_max_indx-floor(floor(ws/dt)/2)):(apos_max_indx+floor(floor(ws/dt)/2));

    nelement1=numel(indx_range(indx_range<=0));
    nelement2=numel(indx_range(indx_range>N_2));
    if nelement1>0
        indx_range=indx_range+nelement1;
    end

    if nelement2>0
        indx_range=indx_range-nelement2;
    end
    
    indx_range=indx_range(indx_range>0);
    indx_range=indx_range(indx_range<=N_2);
    
end
