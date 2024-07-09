function paths_out = rework_paths(paths_in) % with this the tool should run on any operationg systems
    string_in = false;
    struct_in = false;
    if isstring(paths_in)
        paths_in = {paths_in};
        string_in = true;
    elseif isstruct(paths_in)
        paths_fieldnames = fieldnames(paths_in);
        paths_in = struct2cell(paths_in);
        struct_in = true;
    end
    paths_out = cell(size(paths_in)); 
    for i_dir = 1:numel(paths_in)
        path_components = {};
        cur_dir = paths_in{i_dir};
        split1 = strsplit(cur_dir,'\'); % split at every backslash
        for ii = 1:numel(split1)
            path_components = [path_components, strsplit(split1{ii},'/')]; % split at every slash
        end
        paths_out{i_dir} = path_components{1};
        for ii = 2:numel(path_components)
            paths_out{i_dir} = fullfile(paths_out{i_dir}, path_components{ii});
        end
    end
    if string_in
        paths_out = paths_out{1};
    end
    if struct_in
        paths_out = cell2struct(paths_out,paths_fieldnames,1);
    end
end