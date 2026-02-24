Generate Figure 2(a)  
'''
axs[0].plot(x, y2, color='orange', linewidth=1.6, linestyle='-')
axs[0].scatter(30, g1(30), color='red', s=20, zorder=5)
axs[0].axhline(y=threshold, color='red', linewidth=1.6, linestyle='--', label='Bottoming')   
if len(unstable_indices) > 0:
    axs[0].fill_between(x, 0, 5.1, where=stability_results['is_unstable'],
                        color='red', alpha=0.2, label='Rayleigh Unstable')  
axs[0].set_xlabel('$I/I_0$', fontsize=18)
axs[0].set_ylabel('$h/R_0$', fontsize=18)
axs[0].set_title('(a)', loc='left', fontsize=17)
axs[0].legend(loc='upper left', fontsize=13)
axs[0].grid(True, linestyle='--', alpha=0.3)
axs[0].set_ylim(0, 1.4)
axs[0].set_xlim(-3, 70)
axs[0].xaxis.set_major_locator(plt.MaxNLocator(6))
axs[0].yaxis.set_major_locator(plt.MaxNLocator(6))
axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
'''     
Generate Figure 2(b)  
'''
axs[1].plot(x, y1 / np.pi, color='blue', linewidth=1.6, linestyle='-')
axs[1].axvline(x=64, color='red', linewidth=1.6, linestyle='--',label='Bottoming')
if len(unstable_indices) > 0:
    axs[1].fill_between(x, 0, 1.2, where=stability_results['is_unstable'],
                        color='red', alpha=0.2, label='Unstable')
axs[1].set_xlabel('$I/I_0$', fontsize=18)
axs[1].set_ylabel('$\\theta/\\pi$', fontsize=18)
axs[1].set_title('(b)', loc='left', fontsize=17)
axs[1].legend(loc='upper left', fontsize=13)
axs[1].grid(True, linestyle='--', alpha=0.3)
axs[1].set_ylim(0.75, 1)
axs[1].set_xlim(-3, 70)
axs[1].xaxis.set_major_locator(plt.MaxNLocator(6))
axs[1].yaxis.set_major_locator(plt.MaxNLocator(6))
axs[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) 
'''
Generate Figure 2(c) 
'''
axs[2].plot(x, y3 / np.pi, color='green', linewidth=1.6, linestyle='-') 
axs[2].axvline(x=64, color='red', linewidth=1.6, linestyle='--',label='Bottoming')
if len(unstable_indices) > 0:
    axs[2].fill_between(x, 0, 1.2, where=stability_results['is_unstable'],
                        color='red', alpha=0.2, label='Unstable')
axs[2].set_xlabel('$I/I_0$', fontsize=18)
axs[2].set_ylabel('$\\beta/\\pi$', fontsize=18)
axs[2].set_title('(c)', loc='left', fontsize=17)
axs[2].legend(loc='upper right',fontsize=13)  
axs[2].grid(True, linestyle='--', alpha=0.3)
axs[2].set_ylim(0.5, 1.05)
axs[2].set_xlim(-3, 70)
axs[2].xaxis.set_major_locator(plt.MaxNLocator(6))
axs[2].yaxis.set_major_locator(plt.MaxNLocator(6))
axs[2].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
'''
4.Generate the illustration in Figure 2(a)  
'''
V0 = 1.5e-9  
rho = 1e3 
g = 9.8  
r = 2.75e-4  
d = 2e-4  
R = 8.314  
T = 298.15  
gamma_1 = 7.197e-2  
n0 = 1400  
mu0 = np.pi * 4e-7  
m0 =2.7e-7  
m=8e-4
P0=1.013e5
a = 2.8e-3  
D = 3.7e-3  
H = 1e-3  
R0 = (3 * V0 / (4 * np.pi)) ** (1 / 3)
theta0 = np.pi - np.arcsin(r / R0)
Pl=gamma_1*2/R0+P0
def gg2(g2):
    idx = np.argmin(np.abs(x - g2))  
    b1=ARG[:, 0]%(2*np.pi)
    b3=ARG[:, 2]%(2*np.pi)
    theta = b1[idx]  
    beta = b3[idx]   
    I=g2*I00
    N=10000
    h=np.linspace(0,1e-3,N) 
    def f(h):
        term1=-m0*g*(h+d*np.cos(beta))
term2=gamma_1*((2*np.pi*r**2)/(1+np.cos(theta))+np.pi*(r+d*np.sin(beta))*np.sqrt((r-d*np.sin(beta))**2+h**2))
        term3=-(n0*mu0*I*a**2*m)/(2*((D-h-d*np.cos(beta))**2+a**2)**(3/2)) 
        term4_1=np.pi*r**2*H-np.pi*h*(r**2+(d*np.sin(beta))**2+r*d*np.sin(beta))/3-np.pi*d**3*(2+3*np.cos(beta)-np.cos(beta)**3)/3
        term4_2=np.pi*r**2*H-np.pi*r**3*(2+3*np.cos(theta0)-np.cos(theta0)**3)/(3*np.sin(theta0)**3)
        n=(np.pi*r**2*H-np.pi*r**3*(2+3*np.cos(theta0)-np.cos(theta0)**3)/(3*np.sin(theta0)**3))/0.0245
        term4_3=-Pl*(np.pi*h*(r**2+(d*np.sin(beta))**2+r*d*np.sin(beta))/3+np.pi*d**3*(2+3*np.cos(beta)-np.cos(beta)**3)/3-np.pi*r**3*(2+3*np.cos(theta0)-np.cos(theta0)**3)/(3*np.sin(theta0)**3))
        term4=-n*R*T*np.log(term4_1/term4_2)+term4_3
        return term1 + term2 + term3  + term4
    y1=f(h)
    min_index = np.argmin(y1)
    min_h = h[min_index] / R0
    min_energy = y1[min_index]
    plt.plot(h/R0, y1, 'black', label='$I/I_0=30$', linewidth=1.5)
    plt.plot(min_h, min_energy, 'ro', markersize=5)  
    plt.rc('legend', fontsize=14)  
    plt.ylabel('E(J)', fontsize=17, labelpad=0)
    plt.xlabel('h/R0', fontsize=17, labelpad=0)
    plt.tick_params(axis='both', which='major', pad=0)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.gca().yaxis.get_offset_text().set_fontsize(16)
    plt.gca().yaxis.get_offset_text().set_fontweight('bold')
    plt.gca().yaxis.get_offset_text().set_color('black')
    plt.subplots_adjust(bottom=0.5, left=0.5)
    plt.legend()
    plt.tight_layout()
gg2(30)
plt.savefig('E.png', dpi=150)
'''
5.Generate the illustration in Figure 2(c)  
'''
def gg1(g1):
    idx = np.argmin(np.abs(x - g1))  
    b1=ARG[:, 0]%(2*np.pi)
    b2=ARG[:, 1]*1e-4
    b3=ARG[:, 2]%(2*np.pi)
    theta_val = b1[idx]  
    h_val = b2[idx]  
    beta_val = b3[idx]   
    R00 = r / np.sin(theta_val)  
    h0 = (3*r*np.sin(theta_val)**3/(4*(2-3*np.cos(theta_val)+np.cos(theta_val)**3))-r/np.tan(theta_val))  
    def getx1(yy1):
        return np.sqrt(R00**2-(yy1-h0)**2)
    yy1=np.arange(0,h0+R00,1e-8)
    xx1=getx1(yy1)
    plt.plot(xx1,yy1,color="blue", linewidth=1.5)
    plt.plot(-xx1,yy1,color="blue", linewidth=1.5)
    def getx2(yy2):
        d1=2e-4
        return r-(d1*np.sin(beta_val)-r)*yy2/h_val
    yy2=np.arange(-h_val,0,1e-8)
    xx2=getx2(yy2)
    plt.plot(xx2,yy2,color="blue", linewidth=1.5)
    plt.plot(-xx2,yy2,color="blue", linewidth=1.5)
    def getx3(yy3):
        d1=2e-4
        hs=-h_val-d1*np.cos(beta_val)
        return np.sqrt(d1**2-(yy3-hs)**2)
    yy3=np.arange(-H,H,1e-8)
    xx3=getx3(yy3)
    plt.plot(xx3,yy3,color="red", linewidth=1.5)
    plt.plot(-xx3,yy3,color="red", linewidth=1.5)
    x8=np.ones(100)*r
    y9=np.linspace(0,-H,100)
    plt.plot(x8,y9,color="black",lw=1.5)    
    plt.plot(-x8,y9,color="black",lw=1.5) 
    x88=np.linspace(-r,r,100)
    y99=np.ones(100)*(-H)
    plt.plot(x88,y99,color="black",lw=1.5)
    x888=np.linspace(r,4*r,100)
    plt.plot(x888,np.zeros(len(x888)),color="black",lw=1.5)  
    plt.plot(-x888,np.zeros(len(x888)),color="black",lw=1.5) 
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.gca().set_aspect('equal')  
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.gca().xaxis.get_offset_text().set_fontsize(11)
    plt.gca().xaxis.get_offset_text().set_fontweight('bold')
    plt.gca().xaxis.get_offset_text().set_color('black')
    plt.gca().yaxis.get_offset_text().set_fontsize(11)
    plt.gca().yaxis.get_offset_text().set_fontweight('bold')
    plt.gca().yaxis.get_offset_text().set_color('black')
gg1(30)
plt.savefig('I0.png', dpi=150, bbox_inches='tight')
'''
