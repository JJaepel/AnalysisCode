

function plotcycRes2(cyc,cycn)
[~,nt,nl] = size(cyc);
inds = [1:nt];

for ii = cycn
    xt = (1:nl) + (5+nl)*(ii-cycn(1));
    if ii==cycn(1)
        hold off
    end
    x = (ones(length(inds),1)*xt);
    y = nanmean(squeeze(cyc(ii,:,:)),1);
    errBar = nanstd(squeeze(cyc(ii,:,:)),1)./sqrt(nt); %%%
    shadedErrorBar(xt,y,errBar)
    hold on
end
end
