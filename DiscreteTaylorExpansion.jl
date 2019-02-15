using SymPy

'''
Expands a taylor series about a point at the boundary in order to infer its derivative
(The program works wonderfully - the concept not so much.)

n is the number of terms, order is the desired order of the derivative
'''

n = 7
order = 2
x₀ = 0
@vars x Δx
@symfuns f
u = zeros(Sym,n)
c = zeros(Sym,n)
M = zeros(Sym,(n,n))


u = [Sym("u$(i-1)") for i in 1:n]
c = [Sym("c$(i-1)") for i in 1:n]
f = [x^(i-1) for i in 1:n]
g(x) = sum(f.*c)
G = diff(g(x), x, order)
out =  G |> subs(x => x₀)

for i in 1:n, j in 1:n
    M[i,j] = f[j] |> subs(x => (i-1)*Δx)
end
v = M\u
result = out |> subs([c[i] => v[i] for i in 1:n]...)
#really = as_coefficients_dict(result, u...)
