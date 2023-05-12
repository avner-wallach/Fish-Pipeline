for i=1:numel(ops.seg)
    sesspath=[ops.datapath,'\',num2str(ops.seg(i).dates),'\']
    ops.seg(i).eodnum=0;
    ops.seg(i).framenum=0;
    for j=1:numel(ops.seg(i).files)
        filenum=num2str(ops.seg(i).files(j));
        load([sesspath,'data_',filenum]);
        ops.seg(i).framenum=ops.seg(i).framenum+numel(data.FRAME.t);
        ops.seg(i).eodnum=ops.seg(i).eodnum+numel(data.EOD.t); 
    end
    save('output.mat','ops','-append');
end