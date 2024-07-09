function swcdata = read_swcdata(fname)

swcdata = zeros(0,7);

file = fopen(fname);
while(~feof(file))
    fline = fgetl(file);
    
    if fline(1) == '#'
        continue;
    end
    
    A = sscanf(fline, '%d %d %f %f %f %f %d');
    swcdata = [swcdata; A'];
end
fclose(file);
end