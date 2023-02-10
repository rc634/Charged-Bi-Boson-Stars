import numpy as np
import matplotlib.pyplot as plt
import csv

def load_array(name):
	array = []
	file = open(name)
	file = csv.reader(file, delimiter=',')
	for thing in file:
		array += [float(thing[0])]
	array = np.array(array)
	return array

p_array = load_array("p.csv")
f_array = load_array("f.csv")
z_array = load_array("q.csv")
psi_array = load_array("psi.csv")
omega_array = load_array("omega.csv")
r_array = load_array("r.csv")
m_array = load_array("adm.csv")
ham_array = load_array("ham.csv")
N_array = load_array("n.csv")



plt.plot([r_array[0],r_array[-1]],[0,0], linestyle='dashed', color='grey', linewidth='0.8')
plt.plot(r_array,psi_array, 'r', label = r'$\sqrt{g_{rr}}$')
plt.plot(r_array,omega_array, 'b', label = r'$\sqrt{-g_{tt}}$')
plt.plot(r_array,ham_array, 'orange', label = r'$\mathcal{H}$')
plt.plot(r_array,N_array, 'darkorchid', label = r'${\rm d} {N}/{\rm d} r$')
#plt.plot(r_array,f_array**2 + g_array**2, 'brown', label = r'$\rho$')
plt.plot(r_array,m_array, 'g', label = r'$ADM(r)$')
plt.plot(r_array,z_array, 'limegreen',  label = r'$A_t$')#color = [119./255.,1,0])
plt.plot(r_array,p_array, 'k', label = r'$\phi$')#color = [119./255.,0,1])
plt.plot(r_array,f_array, 'gold', label = r'$\psi$')#color = [119./255.,0,1])
plt.legend()

plt.xlim(0,r_array[-1])
plt.ylim(-0.5,2.5)
plt.xlabel('Radius')

plt.show()
