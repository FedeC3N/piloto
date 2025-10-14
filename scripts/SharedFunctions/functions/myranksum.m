function cpval = myranksum(group1,group2)

cpval = nan(1,size(group1,2));

for icomp = 1:size(group1,2)
  
    dummy_comp1 = group1(:,icomp);
    dummy_comp2 = group2(:,icomp);
    
    cpval(icomp) = ranksum(dummy_comp1,dummy_comp2);
    
end



end