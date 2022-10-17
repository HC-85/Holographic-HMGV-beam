function vy = Jpm_f(kf,u,q,mc,nmax, n)
    temp_y = Jpm(kf,u,q,mc,nmax);
    vy = temp_y(n);
end