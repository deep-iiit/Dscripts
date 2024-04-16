%% FUNCTION DEFINITION
% Matlab function to add prefix as per id number
% Author : Deepak Sharma
% Date : October, 2019

function pre = pref_set(i, pref, suf)
    i = double(i);
    if (i < 10)
        pre = strcat(pref, '00', num2str(i), suf);
    elseif (i >= 10 && i < 100)
        pre = strcat(pref, '0', num2str(i), suf);
    elseif (i >= 100 && i < 1000)
        pre = strcat(pref, '', num2str(i), suf);
    end
    return
end
%%