function [ output_args ] = numapprx_LLlogwrapper_DDM_RTlps(Para,nbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin,savefile)
% to be called by optimization function, to record every step of calling
% input par for non-negative guys are log -- need to be transferred back to
% normal space
persistent funtable timeRec
if isempty(Para)
    output_args = [funtable,timeRec];
    timeRec = [];
    funtable = []; % so next time calling you'll get a new table
else
    feature accel on
    tstart = tic;
    Para(2:end) = exp(Para(2:end));
    
    output_args = Fun_LL_DDM_welbullRT(Para,nbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin);
    if isinf(output_args)
        a=1;
    end
    funtable = [funtable;Para,output_args];
    Tcalc = toc(tstart);
    timeRec = [timeRec;Tcalc];
    if ~isempty(savefile)
        if exist(savefile,'file')
            save(savefile,'funtable','timeRec','-append')
        else
            save(savefile,'funtable','timeRec')
        end
    end
end

end

