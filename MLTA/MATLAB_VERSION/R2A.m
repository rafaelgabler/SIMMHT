function R2A = R2A(y, y_hat, p)
    % Calculate R-squared
    ss_total = sum((y - mean(y)).^2);
    ss_residual = sum((y - y_hat).^2);
    r_squared = 1 - (ss_residual / ss_total);

    % Number of observations
    n = length(y);

    % Calculate adjusted R-squared
    R2A = 1 - ((1 - r_squared) * (n - 1) / (n - p - 1));

end
