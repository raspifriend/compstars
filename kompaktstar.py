import numpy as np
import matplotlib.pyplot as plt

rho_c=2*10**14;
G=6.674*10**(-8);
K = 5.3802*10**9;

dr = 1000;
rho=rho_c;
p_c=K*rho_c**(5./3.);
rho_tot = rho_c;
print(p_c, rho)

p=p_c-rho_c**2*2./3.*np.pi*G*dr**2
print(p)
rho=p**(3./5.)/((5.3802*10**9)**(3./5.))
print(rho)
rho_tot += rho
print(rho_tot)
n=2

while True:
	dp = 2./3.*np.pi*G*rho_tot*rho**2*dr**2*(n**3-(n-1)**3)/(n*(n-1))
	p = p-dp
	print(p,dp)
	rho=(p**(3./5.))/((5.3802*10**9)**(3./5.))
        rho_tot += rho
	n += 1

    
    
