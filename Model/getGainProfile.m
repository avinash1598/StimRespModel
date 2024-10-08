function gainProfiles = getGainProfile(nNeurons, timebins, stimDur)
    x = linspace(0, stimDur, length(timebins));
    cumulativeGaussians = zeros(nNeurons, length(x));
    
    % Parameters adjusted so that cum distribution is between 0 and 1
    means = 0.4 + 0.1*randn(1, nNeurons);   % Random means from normal distribution
    stdDevs = 0.08*rand(1, nNeurons); % Random std devs (from uniform, shifted to avoid near-zero values)

    for i = 1:nNeurons
        cumulativeGaussians(i, :) = normcdf(x, means(i), stdDevs(i));
    end

    gainProfiles = cumulativeGaussians;
end