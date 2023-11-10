function  z = loglikelihood(observed, predicted, varp, mode, para)
    
    N = length(observed);

    if mode == "no"
        err = predicted - observed;
        z = -0.5*N*log(2*pi*varp) - 0.5*sum(err.^2)/varp;
    elseif mode == "log"
        err = log(predicted + para(1)) - log(observed + para(1));
        z = -0.5*N*log(2*pi*varp) - 0.5*sum(err.^2)/varp - sum(log(observed+lambda_a));
    end