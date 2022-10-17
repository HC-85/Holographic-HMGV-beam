function vy = Spm_f(kf,v,mc,nmax, n)
    temp_y = Spm(kf,v,mc,nmax);
    vy = temp_y(n);
end