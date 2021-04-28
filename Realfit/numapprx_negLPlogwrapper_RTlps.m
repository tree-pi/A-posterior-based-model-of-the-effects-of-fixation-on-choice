function [ output_args ] = numapprx_negLPlogwrapper_RTlps(Para,Ufun,nUbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin,savefile)
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
    Para([1,3,5,6,7,8]) = exp(Para([1,3,5,6,7,8]));
    
    output_args = -Fun_LL_PUC_welbullRT(Para,Ufun, nUbin, FixNumLNR, LRating, RRating, Choice,ReactionTime,RTbin);
  
    funtable = [funtable;Para,-output_args(1)];
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

