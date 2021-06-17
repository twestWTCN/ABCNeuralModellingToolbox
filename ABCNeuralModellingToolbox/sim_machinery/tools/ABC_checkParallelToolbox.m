function ABC_checkParallelToolbox
        try
            a = gcp;
            if ~a.Connected
            parpool
            end
        catch
            warning('You should try and get access to the parallel computing toolbox to speed things up!')
        end