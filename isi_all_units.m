function [isi_strct] = isi_all_units(spk_t,spk_c,ops,segnum)
%XCOR_ALL_UNITS compute autocorrelation and crosscorrelation for all units
%in segment
Tmax=50; %ms
bin_size=1; %ms
bthreshold=7; %ms
edges=[0:bin_size:Tmax];
bins=edge2bin(edges);
spk_t=double(spk_t)/ops.samplerate*1e3;
U=unique(spk_c);
for i=1:numel(ops.seg(segnum).Spkgroups) %go over all groups
    indu=find(floor(U/10)==ops.seg(segnum).Spkgroups(i) & mod(U,10)>1); %all units in this group
    Ug=U(indu);  
    isi_strct.unames{i}=Ug;
    for u1=1:numel(Ug)
        ind1=find(spk_c==Ug(u1));
        isi=diff(spk_t(ind1));
        isi(isi<1)=[];
        h=histcounts(isi,edges,'Normalization','pdf');        
        isi_strct.val{i}(:,u1)=h;
        isi_strct.mode{i}(u1)=mode(round(isi));
        isi_strct.refractory{i}(u1)=find_ref(h,bins);
        isi_strct.burst{i}(u1)=burst_index(isi);
    end
end
isi_strct.bins=bins;

    function r=find_ref(h,bins)
        b=[bins(1):bin_size/10:bins(end)];
        hh=interp1(bins,h,b);
        r=b(find(hh>=(nanmax(h)/2),1));
    end

    function b=burst_index(isi)
        b=(sum(isi(1:end-1)<bthreshold | isi(2:end)<bthreshold))/(numel(isi)-1);        
    end

end

