function burstinds = edgeCorrect(burstinds,winmark)
edgeflag = []; % This marks if burst is at edge of window
for seg = 1:numel(burstinds)
    if any(intersect(burstinds{seg},find(winmark)))
        edgeflag(seg) = 1;
    else
        edgeflag(seg) = 0;
    end
end
burstinds(edgeflag==1) = [];