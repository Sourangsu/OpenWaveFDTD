function [ f, eps, mu ] = readMatPropsX( filename, p, ip, ipp )
    fid = fopen(filename);
    for i = 1:max([ip ipp]);
        c = textscan(fid,'%f %f','headerlines',2);
        if i == ip
                f = c{1,1}(1:p:end);
                epsp = c{1,2}(1:p:end);
        end
        if i == ipp
                epspp = c{1,2}(1:p:end);
        end
    end
    eps = epsp - 1i * epspp;
    mu = ones(size(eps),'like',eps);
    fclose(fid);
end