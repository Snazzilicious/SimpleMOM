
import sympy as sp


#R = np.linalg.norm( x - y )
#G = np.exp( 1j * k * R ) / ( 4 * np.pi * R )

x1,x2,x3, y1,y2,y3 = sp.symbols("x1 x2 x3 y1 y2 y3", real=True)

R = sp.sqrt( (x1-y1)**2 + (x2-y2)**2 + (x3-y3)**2 )

k = sp.symbols("k",real=True)

G = sp.exp( sp.I*k*R ) / (4*sp.pi*R)

gradG = [ sp.diff(G,x) for x in [x1,x2,x3] ]

gradgradG = [ [ sp.diff(dG,y) for y in [y1,y2,y3] ] for dG in gradG ]

r = sp.symbols("R",real=True)

simple_gradgradG = [ [ sp.simplify( exp.subs( R, r ) ) for exp in rowG ] for rowG in gradgradG ]

# diagonal term
""" ( (xi - yi)**2 * (R**2*k**2 + 3*I*R*k - 3) + R**2 + -I*R**3*k ) * exp(I*R*k)/(4*pi*R**5) """
# Off diagonal
""" (x1 - y1)*(x2 - y2)*(R**2*k**2 + 3*I*R*k - 3)*exp(I*R*k)/(4*pi*R**5) """
# note both real and imag parts of exp(I*R) are even functions because R>0
