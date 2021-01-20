figure(3)
for i=1:2
    for j=1:2
        subplot(2,2,sub2ind([2 2],i,j))
       if i ==j
           ylim([0 6])
       else
           ylim([0 1])
       end
       grid off
       box off
    end
end