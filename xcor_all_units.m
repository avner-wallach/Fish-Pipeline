function [xcor_strct] = xcor_all_units(spk_t,spk_c,ops,segnum)
%XCOR_ALL_UNITS compute autocorrelation and crosscorrelation for all units
%in segment
Tmax=50; %ms
bins=[];
spk_t=double(spk_t)/ops.samplerate*1e3;
U=unique(spk_c);
for i=1:numel(ops.seg(segnum).Spkgroups) %go over all groups
    indu=find(floor(U/10)==ops.seg(segnum).Spkgroups(i) & mod(U,10)>1); %all units in this group
    Ug=U(indu);  
    xcor_strct.unames{i}=Ug;
    for u1=1:numel(Ug)
        ind1=find(spk_c==Ug(u1));
        isi=diff(spk_t(ind1));
        ind1(isi<1)=[];
        [c,bins]=MyXcor(spk_t(ind1),spk_t(ind1),Tmax);
        xcor_strct.val{i}(:,u1,u1)=c;    
        for u2=(u1+1):numel(Ug)
            ind2=find(spk_c==Ug(u2));
            [c,bins]=MyXcor(spk_t(ind1),spk_t(ind2),Tmax);
            xcor_strct.val{i}(:,u1,u2)=c;    
        end
    end
end
xcor_strct.bins=bins;
end

