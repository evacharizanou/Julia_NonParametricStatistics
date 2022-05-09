Pkg.add("Distributions")
using Distributions
function chisq_test(O; a = 0.05)
    if ndims(O) ==1
        d = length(O)
        n = sum(O)
        df = d-1
        E = repeat([n/d] , outer=d)
    else
      ni = sum(O, dims=2)
      nj = sum(O, dims=1)
      n = sum(O)
      E = (ni*nj)./n
      df = (size(O)[1]-1)*(size(O)[2]-1)
   end
   q = sum(((O-E).^2)./E)
   qa = quantile(Chisq(df), 1-a)
   p = 1-cdf(Chisq(df), q)
   (q, qa, p)
end
#EXERCISE 2
O = [ 34 39 19 20 28;
       20 15 23 19 40;
       41 20 42 29 38]
chisq_test(O, a = 0.01)

#EXERCISE 3
id = 233
using Random
Random.seed!(id)

function gf_test(O; p0=nothing, a=0.05)
  O = collect(O)
  n = sum(O)
  d = length(O)
  df = d-1
  if isnothing(p0)
    E = repeat([n/d], outer=d)
  else
   E = n.*p0
  end
  q = sum(((O.-E).^2)./E)
  qa = quantile(Chisq(df), 1-a)
  p = 1-cdf(Chisq(df), q)
  (q, qa, p)
end

n = 200
O = n * rand(Dirichlet(6,2)) #sum(O)=200
p_0 = (0.18, 0.07, 0.32, 0.19, 0.11, 0.13)
gf_test(O, p0=p_0)
