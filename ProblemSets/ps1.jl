#PROBLEM SET 1

using Plots, Printf, LinearAlgebra

# This is a system of non-linear equations
function myf(x)
    F = zeros(size(x));
    F[1] = x[1]^2 + x[2]^2 - 1;
    F[2] = x[1] - x[2];
    return F
end


# This is my function for Broyden's method
function mybroyden(f, x0; tol=1e-4, max_iter=20)
    iter = 0;
    dif  = 1;
    A0   = I(length(x0));

    while dif>tol && iter<max_iter
        iter   += 1;
        x1      = x0 -inv(A0)*f(x0);
        delta_x = x1 - x0;
        delta_f = f(x1) - f(x0);
        A1      = A0 + (delta_f - A0*delta_x)/delta_x;
        dif     = norm(x1-x0);  

        #Update
        x0 = x1;
        A0 = A1;
    end 
    return x0
end

#Initial guess
x_ini=[1, 1];

sol = mybroyden(myf, x_ini);

println("Solution: ", sol)
