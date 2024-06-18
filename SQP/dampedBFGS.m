function B = dampedBFGS(H, p, q)
    % Damped BFGS update
    theta = 1;
    if p'*q < 0.2*p'*H*p
        theta = 0.8*p'*H*p/(p'*H*p - p'*q);
    end

    r = theta*q + (1-theta)*H*p;

    B = H - (H*p)*(p'*H)/(p'*H*p) + (r*r')/(p'*r);
end
