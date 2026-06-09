function E = getE(Em,EC50,Ce)
        E = (Em*Ce) ./ (EC50+Ce);
end