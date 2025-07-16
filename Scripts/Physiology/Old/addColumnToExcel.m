monkey = 'Gilligan';
singleOrAll = 'Single';
cond = 'ESS';
PSTH_type = 'Trial';
PSTHdir = dir(['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\', singleOrAll,'\']);
PSTHdir = PSTHdir(cellfun(@(a) ~isnan(str2double(a)),{PSTHdir.name}));
sessionInfoPath = ['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHS\', singleOrAll,'\FRs\','M1_',cond,'_',PSTH_type,'_Summary.xlsx'];
spsd = spreadsheetDatastore(sessionInfoPath,'Sheet',1,'TextType', 'string');
spsd.ReadSize = 'sheet';
sessionInfo = read(spsd);

for d = 1:length(PSTHdir)
    m = matfile([PSTHdir(d).folder,'/',PSTHdir(d).name,'/',PSTHdir(d).name,'.mat']);
    if(~isempty(m.Corrs))
        condCorr = m.Corrs(:,:,1);
        belowDig = tril(true(size(condCorr)),-1);
        sessionInfo(sessionInfo.('Site #') == str2double(PSTHdir(d).name),:).Xcorr = nanmedian(abs(condCorr(belowDig)));
    end
end
writetable(sessionInfo, ['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHS\', singleOrAll,'\FRs\','M1_',cond,'_',PSTH_type,'_Summary.xlsx'],...
    'WriteVariableNames', 1)