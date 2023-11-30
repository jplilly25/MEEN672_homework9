function a = print_output(var, vals)

    fprintf('%s = \n', var)

    for i = 1:size(vals, 1)
        for j = 1:size(vals, 2)
            fprintf('\t%-.4f', vals(i, j))
        end
        fprintf('\n')
    end

end