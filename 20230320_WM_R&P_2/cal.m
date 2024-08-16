%re_cal_12,output = 24 ~[0,180]
%re_cal_19,output = 25 ~[-90,90]
% By YuanGao.24-Mar-2023 09:01:49

for trial_i = 1:72
    %% re_cal_12,output = 24 ~[0,180]

    while results(trial_i,12) > 180 || results(trial_i,12) < 0
        if results(trial_i,12) > 180
            results(trial_i,12) = results(trial_i,12) - 180;
        elseif results(trial_i,12) < 0
            results(trial_i,12) = results(trial_i,12) + 180;
        end
    end
    results(trial_i,24) = results(trial_i,12);

%% re_cal_19,output = 25 ~[-90,90]

    if results(trial_i,22) == results(trial_i,2)
        error = results(trial_i,24) - results(trial_i,8);
    elseif results(trial_i,22) == results(trial_i,3)
        error = results(trial_i,24) - results(trial_i,9);
    elseif results(trial_i,22) == results(trial_i,4)
        error = results(trial_i,24) - results(trial_i,10);
    end

    if error < (-90)
        results(trial_i,25) = error + 180;
    elseif error > 90
        results(trial_i,25) = error - 180;
    else
        results(trial_i,25) = error;
    end

    results(trial_i,26) = results(trial_i,9) - results(trial_i,8);
    results(trial_i,27) = (results(trial_i,10) - 180) - results(trial_i,8);

end


% % results(results(:, 11) == 1, 25) % error 1
% % results(results(:, 11) == 2, 25) % error 2
% % results(results(:, 11) == 3, 25) % error 3

data.errors = [results(results(:, 11) == 1, 25)];
data.distractors = [results(results(:, 11) == 1, 26); results(results(:, 11) == 1, 27)];
% data.distractors = [results(results(:, 11) == 1, 27)];

MemFit(data, Orientation(SwapModel(), 3))

% % MemFit(data, Orientation(WithBias(StandardMixtureModel), [1,3]))

% e.g.,
%  model = Orientation(StandardMixtureModel(), 2)
%  or
%  model = Orientation(WithBias(StandardMixtureModel), [1 3])
%  or
%  model = Orientation(SwapModel(), 3)