import numpy as np 
import galois

x = galois.GF(2**7, galois.Poly([1, 0, 0, 0, 1, 0, 0, 1], galois.GF(2)))


msg = galois.Poly(np.random.randint(2, size=113), x)

gen = galois.Poly([1,0,0,0,0,1,1,0,1,1,1,0,1,1,1], x)
codeword = msg * gen
codeword.coeffs
cf = codeword.coeffs.copy()
r = galois.Poly(cf, x)
r.coeffs
r.coeffs == codeword.coeffs
r.coeffs[1] = 0


t = 2

syndromes = []

for i in range(2*t):
    syndromes.append(r(gen.roots()[0]**(i+1)))
syndromes

for v in range(t, 0, -1):
    m = []

    print(v)

    for i in range(v):
        m.append(syndromes[i:v+i])

    print(np.linalg.det(x(m)))

    

    if np.linalg.det(x(m)) != 0:
        svs = x([[x] for x in syndromes[v:2*v]])
        svs
        np.linalg.inv(x(m))
        locs = np.linalg.inv(x(m)) @ x(svs)
        locs
        p = np.reciprocal(locs)
        p
        break

syndromes

x.display('power')
x.Elements()



p = galois.Poly([1, 1, 1, 1, 1, 2, 1, 1, 0, 1, 1, 0], field=x)
