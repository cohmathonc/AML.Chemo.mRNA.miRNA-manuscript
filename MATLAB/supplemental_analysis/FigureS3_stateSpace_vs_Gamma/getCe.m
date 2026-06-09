function Ce = getCe(ke0,t,dt,C)
        for n = 1:length(t)
        Ce(n) = exp(-ke0*t(n))*ke0.*((exp(ke0*t(1:n))*C(1:n)')).*dt;
        end
end