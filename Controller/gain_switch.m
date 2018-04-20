function [temp_in, new_val]= gain_switch(temp_in, old_temp_in, old_val_in)
    if temp_in > old_temp_in
        new_val = 1.69122;
    else
        if temp_in < old_temp_in
            new_val = 6.98702;
        else 
            new_val = old_val_in;
        end
    end

