function saveOutcome = setWWPPPath()
    
    if ispc()
        root = '[A-Za-z]\:\\';
    else
        root = '\/';
    end
    filePath            = mfilename('fullpath')                     ;
    packagePath         = [root,'.+?WisconsinWaterPropertyPackage'] ;
    topDirectory        = char(regexp(filePath,packagePath,'match'));
    excludeTopDirectory = true;
    excludedDirectories = {'_Misc';'Diagnostics';'Documentation';'FixUpGenerators';'Fortran90+';'Tests'};
    saveOutcome         = updatePath(topDirectory,excludeTopDirectory,excludedDirectories);

end