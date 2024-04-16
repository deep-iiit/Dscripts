% To run this file provide a realignment parameter file with 6 motion parameters
% The script/function will convert standard 6 motion parameter/realignnment file to 24 Friston regressors

% Generates file with 24 regressors.
% Column no.  1 to 6  are the motion parameters.
% Column no.  7 to 12 are the square of the motion parameters.
% Column no. 13 to 18 are the temporal difference of motion parameters.
% Column no. 19 to 24 are the square of the temporal difference values.

function motion_parameters_to_friston_regressors(x)
    in = readmatrix(x);
    format('long','e') % Standard SPM numeric format for motion parameters
    rows = size(in, 1);
    out = zeros(rows, 24);
    for r = 1:rows
        if r ~= 1
            out(r, 1)  = in(r, 1) ;
            out(r, 2)  = in(r, 2) ;
            out(r, 3)  = in(r, 3) ;
            out(r, 4)  = in(r, 4) ;
            out(r, 5)  = in(r, 5) ;
            out(r, 6)  = in(r, 6) ;
            out(r, 7)  = in(r, 1) ^ 2 ;
            out(r, 8)  = in(r, 2) ^ 2 ;
            out(r, 9)  = in(r, 3) ^ 2 ;
            out(r, 10) = in(r, 4) ^ 2 ;
            out(r, 11) = in(r, 5) ^ 2 ;
            out(r, 12) = in(r, 6) ^ 2 ;
            out(r, 13) = in(r, 1) - in(r - 1, 1) ; % Motion difference between t and (t-1) time point (t-(t-1))
            out(r, 14) = in(r, 2) - in(r - 1, 2) ;
            out(r, 15) = in(r, 3) - in(r - 1, 3) ;
            out(r, 16) = in(r, 4) - in(r - 1, 4) ;
            out(r, 17) = in(r, 5) - in(r - 1, 5) ;
            out(r, 18) = in(r, 6) - in(r - 1, 6) ;
            out(r, 19) = out(r, 13) ^ 2 ;
            out(r, 20) = out(r, 14) ^ 2 ;
            out(r, 21) = out(r, 15) ^ 2 ;
            out(r, 22) = out(r, 16) ^ 2 ;
            out(r, 23) = out(r, 17) ^ 2 ;
            out(r, 24) = out(r, 18) ^ 2 ;
        end   
    end
    out(1, :) = 0;
    save regressors.txt out -ascii -tabs % Save regressors as matlab variable
    present_path = pwd; 
    splits = strsplit(present_path, '\');
    subj_folder = splits{end}; % Get the immediate subject directory name
    reg_filename = strcat(subj_folder, '_regressors.txt');
    movefile('regressors.txt', reg_filename); % Rename the regressor file as per subject name
end


