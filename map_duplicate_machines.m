function map_duplicate = map_duplicate_machines(G, G_init)
% Find a mapping between initial graph G_init, without duplicates, and processed
% graph G that contains duplicated machines.
map_duplicate = [];

    for m=1:max(max(G_init))
        map_duplicate(m,:) = [m m];
    end
    if max(max(G)) > max(max(G_init))
        for m=max(max(G_init))+1:max(max(G))
            [row_dup, col_dup] = find(G == m);
            map_duplicate(m,:) = [m G_init(row_dup,col_dup)];
        end   
    end
end