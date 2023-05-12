function [eodind]=get_eod_times(t,amp,ops)
    if(ops.eodchan==0)            
        if(strcmp(ops.eoddiff,'on'))
            a=[0;max(abs(diff(amp)),[],2)];
        else
            a=max(amp,[],2);
        end
    else
        if(strcmp(ops.eoddiff,'on'))
            a=[0;abs(diff(amp(:,ops.eodchan)))];
        else
            a=amp(:,ops.eodchan);
        end
    end            
    ind=find(diff(a>ops.eodth)==1); %find threshold posedge
    j=2;
    while(j<=numel(ind))
        if((t(ind(j))-t(ind(j-1)))<=ops.eodref/1e3) %within ref period
            ind(j)=[]; %remove element
        else
            j=j+1;
        end
    end   
    eodind=ind;
end
