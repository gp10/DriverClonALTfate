function [Freq_truncated] = freq2squeeze(freqs,freqs_dim_lim,timepoints,vartype)
%% REFORMATS CLONE SIZE FREQUENCIES ASSIGNING HUGE OUTLIER CLONES TO MAXIMUM BINNED SIZE CATEGORY SPECIFIED:
% Due to the difficulty in computing likelihood prob. of a few (<1%) huge
% outlier clones, bins can be defined up to a maximum size, so that all
% clones with size above or equal to that threshold are assigned to the same
% category.

% from Herms et al, 2020

%% Input:
% freqs: cell array {1,n}(:,1) (experimental data) or matrix [m,n] (simulated data) of clone size frequencies binned in ranges increasing in powers of 2 (n = No. of time points, m = No. of clone size categories)
% freqs_dim_lim: max. category of clone sizes to be considered (all clones of size >= [2^(freqs_dim_lim-3)+1] are assigned to that category)
% timepoints: vector of time points (expressed in weeks)
% vartype: (1= Experimental data | 0= Simulated data)

%% Output:
% Freq_truncated: cell array {1,n}(:,1) (experimental data) or matrix [m,n] (simulated data) of clone size frequencies binned in ranges increasing in powers of 2 after outlier clones are reassigned into max. category

%% Calculations:
% cell No. are already categorized into ranges increasing in size in powers of 2
switch vartype

    case 1
        % Experimental data:
        Freq_truncated = freqs;
        for aa = 1:size(timepoints,2)
            if size(freqs{:,aa},1)>=freqs_dim_lim
                cap_nx_total = [];
                cap_nx_total = sum(freqs{:,aa}(freqs_dim_lim:end,1),1);
                Freq_truncated{:,aa}(freqs_dim_lim:end,1) = 0;
                Freq_truncated{:,aa}(freqs_dim_lim,1) = cap_nx_total;
            end
        end

    case 2
        % Simulated data (ifneedbe):
        Freq_truncated = freqs;
        if size(freqs,1)>=freqs_dim_lim
            cap_nx_total = sum(freqs(freqs_dim_lim:end,:),1);
            Freq_truncated(freqs_dim_lim:end,:) = 0;
            Freq_truncated(freqs_dim_lim,:) = cap_nx_total;
        end

end
