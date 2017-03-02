function resMatrix = getFreq(f0,r,resMatrix)
%%
% resMatrix = getFreq(f0,r,resMatrix)
%
% resMatrix may be an array of structure containing the results
% f0 and r are vector of the same size as resMatrix respectively containing
% the initial frequency and sweep rate
% The output is the initial cell array in which a field 'f' containg the
% frequencies had been added in each structure
%%

    for i = 1:length(resMatrix)
        freq = f0(i) + (r(i)/60).*resMatrix{i}.t;
        resMatrix{i}.f = freq;
    end
end