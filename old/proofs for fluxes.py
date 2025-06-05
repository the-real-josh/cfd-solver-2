import sympy as sym

rho, u, v, E, p, cp, cv, gamma, R, T, e = sym.symbols('rho, u, v, E, p, cp, cv, gamma, R, T, e')
q0 = rho
q1 = rho*u
q2 =rho*v
q3 = rho*E

relations = {
    R: cp-cv,
    p: rho*R*T,
    gamma: cp/cv,
    E: 0.5*(u**2 + v**2) + e,
    e: cv*T
}
# other transformations - everything goes to cp, cv, u, v, and T
# p = rho*R*T
# R = cp-cv
# gamma = known
# cp/cv = gamma
# E = 1/2 * (u**2 + v**2) + e
# e = cp*T

f0 = rho*u
f1 = rho*u**2 + p
f2 = rho*u*v
f3 = rho*E*u + p*u

g0 = rho*v
g1 = rho*u*v
g2 = rho*v**2 + p
g3 = rho*E*v + p*v


boil = lambda a: sym.simplify(a.subs(relations).subs(relations).subs(relations)) 

# need to represnet f0...f3 in terms of q
# need to represnet g0...g3 in terms of q

# proofs for f
f0s = q1
print(f'proof: {sym.simplify(f0s.subs(relations)/f0.subs(relations))}=1')

f1s = q1*q1/q0 + (q3-0.5*(q1**2 + q2**2)/q0)*(gamma-1)
print(f'proof: {boil(f1s)/boil(f1)}=1')

f2s = q1*q2/q0
print(f'proof: {boil(f2s)/boil(f2)}=1')

f3s = q3*q1/q0 + (q3-0.5*(q1**2 + q2**2)/q0)*(gamma-1)*q1/q0
print(f'proof: {boil(f3s)/boil(f3)}=1')


# proofs for g
g0s = q2
print(f'proof: {sym.simplify(g0s/g0)}=1')

g1s = q1*q2/q0
print(f'proof: {boil(g1s)/boil(g1)}=1')

g2s = q2*q2/q0 + (q3-0.5*(q1**2 + q2**2)/q0)*(gamma-1)
print(f'proof: {boil(g2s)/boil(g2)}=1')

g3s = q3*q2/q0 + (q3-0.5*(q1**2 + q2**2)/q0)*(gamma-1)*q2/q0
print(f'proof: {boil(g3s)/boil(g3)}=1')

