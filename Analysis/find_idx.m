%%
function ii=find_idx(list_val,val)

for ii=1:length(list_val)
    
    if val==list_val(ii)
        return;
    end
    
end

error('could not find index for requested value.')

end