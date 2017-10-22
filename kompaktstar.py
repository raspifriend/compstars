import numpy as np
import pylab
from scipy.integrate import quad
from matplotlib import pyplot as plt

# Definiere Konstanten (Einheiten im cgs-System)
rho_c = 2.0*10**14	# zentrale Dichte
G = 6.674*10**(-8)	# Gravitationskonstange
c = 29979245800.0	# Lichtgeschwindigkeit
K = 5.3802*10**9	# aus Shapiro, Teukolsky, nicht-relativistisches Neutronengas (nrel-ngas)
A = 2./3.*np.pi*G 	# mehrmals verwendete Konstante
g = 5./3.		# Gamma der Zustandsgleichung (EOS) p=rho^Gamma für nrel-ngas	
b = 3./5.		# 1/Gamma

# Berechne Druck im Zentrum mit EOS
rho = rho_c
p_c = K*rho_c**g

# Schrittweite der Integration nach r (in cm)
dr = 100.0

# Definiere Funktion zur Berechnung von dm und dp (Sternstruktur-Gleichungen)
def mass(r, d):
    return 4.0*np.pi*r**2.0*d
def press(r, m, d):
    return (1/(r**2.0))*(-G*m*d)
# dimensionslose Grösse
def dimless(M, r):
    return 1-(2.0*G*M)/(c**2.0*r)

# 1. Integrationsschritt von r=0 nach r=dr mittels Näherung rho=rho_c (aus Vorlesung [26.9.17])
p = p_c-rho_c**2.0*A*dr**2.0
rho = (p/K)**b
M = 0.0

# Definiere Vektoren r, p, m und rho zum Plotten
rvec = [0.0, dr]
pvec = [p_c, p]
mvec = [0.0, 4./3.*np.pi*dr**3.0*rho_c]
rhovec = [rho_c]

n = 2.0	# n. Integrationsschritt

# numerische Integration zur schrittweise Berechnung von rho(r) und p(r)
while p>=0:
    # berechne Dichte rho
    rho = (p/K)**b
    # berechne Masse der n-ten Massenschale durch Integration
    m = quad(mass,(n-1.0)*dr,n*dr,args=(rho))
    # addiere die neuberechnete Massenschale zur Gesamtmasse M
    M += m[0]
    # berechne Druckänderung bei der n-ten Massenschale durch Integration
    P = quad(press,(n-1.0)*dr, n*dr, args=(M, rho))
    # neuer Druck berechnen
    p = p+P[0]
    # neu berechneter Druckwert p an Vektor anhängen
    pvec = np.append(pvec, p)
    # r-Vektor um r+dr erweitern
    rvec = np.append(rvec, rvec[-1]+dr)
    # m-Vektor um aktuelle Masse erweitern
    mvec = np.append(mvec, M)
    # neu berechnete Dichte rho an Vektor anhängen
    rhovec = np.append(rhovec, rho)
    #Integrationsschritt n um 1 erhöhen
    n += 1	

# negative Druckwerte wegschneiden (ebenso dazugehörige Werte von r und M)
rvec = np.delete(rvec,-1)
pvec = np.delete(pvec,-1)
mvec = np.delete(mvec,-1)

# Dimensionslose Grösse 1-2Gm/(c^2*r) berechnen für r = 0-19 km in 1km-Schritten
dlvec = [1]
rrvec = [0]
for i in range(1, 20):
    a = i*1000
    dl = dimless(mvec[a], rvec[a])
    dlvec = np.append(dlvec, dl)
    rrvec = np.append(rrvec, rrvec[-1]+1)

# Daten in txt-Datei abspeichern
data = np.array([rvec, rhovec, pvec, mvec])
data = data.T
np.savetxt('data.txt', data)
np.savetxt('dimless.txt', dlvec)

# Plot des Dichte-Profils
plt.figure(0)
plt.xlabel(r"Radius $r$ [km]")
plt.ylabel(r"Dichte $\rho$ [$\cdot 10^{14}\frac{g}{cm^3}$]")
plt.title('Dichte-Profil des Neutronensterns')
plt.plot(1.0*10**(-5)*rvec, 1.0*10**(-14)*rhovec) # von cm --> km
pylab.savefig('dichteprofil.png', dpi = 300)

# Plot die dimensionslose Grösse 1-2Gm/(c^2*r)
plt.figure(1)
plt.xlabel(r"Radius $r$ [km]")
plt.ylabel(r"$1-\frac{2Gm}{c^2r}$")
plt.title('Dimensionslose Grösse')
plt.plot(rrvec, dlvec, 'x')
pylab.savefig('dimless.png', dpi = 300)

# Felix Läderach


    
    
