function p = bisection(f,a,b,r,theta)
% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code. But may come handy
if f(a,r,theta)*f(b,r,theta)>0 
    disp('Wrong choice bro')
else
    p = (a + b)/2;
    err = abs(f(p,r,theta));
    while err > 1e-7
      if f(a,r,theta)*f(p,r,theta)<0 
          b = p;
      else
          a = p;          
      end
      p = (a + b)/2; 
      err = abs(f(p,r,theta));
    end
end