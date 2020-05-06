function [c,bins]=MyXcor(x,y,T)
bin_size=0.1; %ms
x=x(:); y=y(:);
edges=[-T+bin_size/2:bin_size:T-bin_size/2];
bins=edges+bin_size/2;
h=zeros(size(edges));
K=1e6;
z=nan(1,K);
k=1;
p=1;
P=1e4;
if(numel(x))
    while((p+P)<numel(x))
        xx=x(p:p+P-1);
        yy=y(y>=xx(1)+edges(1) & y<xx(end)+edges(end));
        for i=1:length(xx)
            yyy=yy(yy>=(x(i)+edges(1)) & yy<(x(i)+edges(end)) & yy~=x(i))'-x(i);
            z(k:k+numel(yyy)-1)=yyy;        
            k=k+numel(yyy);
            if(k>K-1e3)
                h=h+histc(z,edges);
                z=nan(1,K);
                k=1;
            end
        end
        p=p+P;
    end 
    if(k>1)
        h=h+histc(z,edges);
    end    
    c=h/numel(x);    
    c(abs(bins)<bin_size/2)=nan;
else
    c=nan(size(bins));
end