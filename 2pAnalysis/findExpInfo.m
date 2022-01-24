function exp_info = findExpInfo(xls_txt, xls_all)

animalCol = find(contains(xls_txt(1,:),'animal'),1);
expCol = find(contains(xls_txt(1,:),'expNumber'),1);
spk2Col = find(contains(xls_txt(1,:),'spk2'),1);
volCol = find(contains(xls_txt(1,:),'volume'),1);
nameCol = find(contains(xls_txt(1,:),'name'),1);
runCol = find(contains(xls_txt(1,:),'Run'),1);
flagCol = find(contains(xls_txt(1,:),'Flag'),1); 
try
    specCol = find(contains(xls_txt(1,:),'special'),1);
end
try
    regCol = find(contains(xls_txt(1,:),'region'), 1);
end
try
    stimCol = find(contains(xls_txt(1,:),'stimulus'),1);
end
try
    comCol = find(contains(xls_txt(1,:), 'comments'),1);
end
try
    EOCol = find(contains(xls_txt(1,:), 'EO'),1);
end
try
    IncludedCol = find(contains(xls_txt(1,:), 'Include'),1);
end

exp_info = struct;

k = 1;
for i = 2:size(xls_txt,1)
    exp_info.animal{k} = xls_all(i,animalCol);
    exp = xls_all(i,expCol);
    if iscell(exp)
        exp = cell2mat(exp);
        if ischar(exp)
            exp = str2double(exp);
        end
    end
    if exp > 9
        if exp > 99
            exp_info.exp_id{k} = ['t00' num2str(exp)];
        else
            exp_info.exp_id{k} = ['t000' num2str(exp)];
        end
    else
        exp_info.exp_id{k} = ['t0000' num2str(exp)];
    end
    exp_info.exp_series{k} = ['tseries_' num2str(exp)];
    sp2 = xls_all(i,spk2Col);
    if iscell(sp2)
        sp2 = cell2mat(sp2);
        if ischar(sp2)
            sp2 = str2double(sp2);
        end
    end
    if sp2 > 9
        if sp2>99
            exp_info.sp2_id{k} = ['t00' num2str(sp2)];
        else
            exp_info.sp2_id{k} = ['t000' num2str(sp2)];
        end
    else
        exp_info.sp2_id{k} = ['t0000' num2str(sp2)];
    end
    vol = xls_all(i,volCol);
    try
        if contains(vol{1}, 'yes')
            exp_info.vol{k} = 1;
        elseif contains(vol{1}, 'no')
            exp_info.vol{k} = 0;
        end
    end
    exp_info.name{k} = xls_all(i,nameCol);
    runInd = xls_all{i,runCol};
    if ischar(runInd)
        runInd = str2double(runInd);
    end
    exp_info.run(k) = runInd;
    flagInd = xls_all{i,flagCol};
    if ischar(flagInd)
        flagInd = str2double(flagInd);
    end
    exp_info.flag(k) = flagInd;
    try
        specInd = xls_all{i,specCol};
        if ischar(specInd)
            specInd = str2double(specInd);
        end
        exp_info.special(k) = specInd;
    end
    try
        regInd = xls_all{i, regCol};
        if ischar(regInd)
            regInd = str2double(regInd);
            if isnan(regInd)
               regInd = xls_all{i, regCol};
            end
        end
        exp_info.region{k} = regInd;
    end
    try 
       exp_info.stimulus{k} = xls_all{i, stimCol};
    end
    try
       exp_info.comments{k} = xls_all{i, comCol};
    end
    try
       exp_info.EO{k} = xls_all{i, EOCol};
    end
    try
       exp_info.included{k} = xls_all{i, IncludedCol};
    end
    k = k+1;
end