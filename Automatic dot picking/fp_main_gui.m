function fp_main_gui()
% GUI to run the FRET pipeline
% Requires DOTTER
% Check fp_main line 9 and 233
% Check fp_getFret line 88

close all
dbstop error

s = fp_main();
f = fields(s);

conf = [];
for kk = 1:numel(f)
    conf(kk).name = f{kk};
    v = s.(f{kk});
    conf(kk).value = v;
    name = conf(kk).name;
    conf(kk).type = 'string';
    switch class(v)
        case 'double'
            conf(kk).type = 'numeric';
        case 'logical'
            conf(kk).type = 'logical';
    end
    if contains(name, 'folder')
        conf(kk).type = 'dir';
    end
    if contains(name, 'dir')
        conf(kk).type = 'dir';
    end
    if contains(name, 'file')
        conf(kk).type = 'file'; 
    end        
end    

    gs = df_guisettings('name', 'Fret Pipeline', 'settings', conf, 'help', @showHelp);
    if numel(gs) == 0
        return;
    end

    s = df_structput(s, gs);
    fp_main(s);    
end

function showHelp(varargin)
    p = mfilename('fullpath');
    [p1, ~] = fileparts(p);
    helpfile = [p1 filesep() 'README.html'];
    web(helpfile);
end