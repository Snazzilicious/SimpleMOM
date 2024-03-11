
import sympy as sp


#R = np.linalg.norm( x - y )
#G = np.exp( 1j * k * R ) / ( 4 * np.pi * R )

x1,x2,x3, y1,y2,y3, dx = sp.symbols("x1 x2 x3 y1 y2 y3 dx", real=True)

R = sp.sqrt( (x1-y1)**2 + (x2-y2)**2 + (x3-y3)**2 )

k = sp.symbols("k",real=True)

G = sp.exp( -sp.I*k*R ) / (4*sp.pi*R)

gradxG = [ sp.diff(G,x) for x in [x1,x2,x3] ]
gradyG = [ sp.diff(G,y) for y in [y1,y2,y3] ]

gradygradxG = [ [ sp.diff(dG,y) for y in [y1,y2,y3] ] for dG in gradxG ]

r = sp.symbols("R",real=True)

simple_gradxG = [ sp.simplify( exp.subs( R, r ) ) for exp in gradxG ]
simple_gradxG[0] = sp.simplify( simple_gradxG[0].subs( x1, dx+y1 ) )
# Grad w.r.t. X
""" -dx*(I*R*k + 1)*exp(-I*R*k)/(4*pi*R**3) """

simple_gradygradxG = [ [ sp.simplify( exp.subs( R, r ) ) for exp in rowG ] for rowG in gradygradxG ]
# diagonal term
""" ( dx**2 * (R**2*k**2 - 3*I*R*k - 3) + I*R**3*k + R**2 )*exp(-I*R*k)/(4*pi*R**5) """
# Off diagonal
""" (x1 - y1)*(x2 - y2)*(R**2*k**2 - 3*I*R*k - 3)*exp(-I*R*k)/(4*pi*R**5) """
# note both real and imag parts of exp(I*R) are even functions because R>0



# Integrating the Green's integral
a,b,r = sp.symbols("a b r",real=True)
sp.integrate( sp.cos(r)*sp.sin(r), (r,0,2*sp.pi) )
sp.integrate( sp.exp(-sp.I*a*r)*sp.besselj(0,b*r) , (r,0,sp.oo) )
# -I*b / sqrt( a**2 - b**2 ), a**2 > b**2 // -b/sqrt(a**2 - b**2), o/w
sp.hankel_transform( sp.exp(-sp.I*a*r)/r, r,b,0)
# exp_polar(3*I*pi/2)/(b*sqrt(a**2/b**2 - 1)), a**2/b**2 > 1
