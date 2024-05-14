import numpy as np
import matplotlib.pyplot as plt
#définition des paramètres
C0=0
Sm = 0.5
S0=Sm 
E0=1
Lambda0=0.1
d=0.2
q=0.8
p=0.5
k=(p/q)*Sm
h=p*Sm
b=2*d/q
Lambda=1.5
Rm=7
mE=0.01
mLambda=0.01

# initialisation de C,S et E
C = np.zeros(300)  
S = np.zeros(300)  
E = np.zeros(300)  
A = np.zeros(300)   
v = np.zeros(300) 
Psi = np.zeros(300) 

S[0]=Sm
E[0]=E0

# itérations 
for i in range(300):
    R = min(np.random.poisson(Lambda), Rm)
    Psi[i]= C[i] -E[i] - S[i]
    v[i]=max(0,min(Psi[i],1))
    A[i]=q*v[i] + (R/Rm)*q*(1-v[i])

    if i<299:
        C[i+1]=(1-d)*C[i]+b*min(1,1-C[i])*A[i]
        S[i+1]=S[i]+p*max(0,Sm-S[i])-h*C[i]-k*A[i]
        E[i+1]=E[i]-mE
    Lambda+=0.0001
        



    
    

T=np.linspace(0,300,300)
# Tracer les graphiques
plt.figure(figsize=(10, 6))
plt.plot(T, A, label='A')
plt.plot(T, S, label='S')
plt.plot(T, C, label='C')
plt.plot(T, E, label='E')
plt.plot(T, v, label='v')
plt.plot(T, Psi, label='Psi')
plt.xlabel('Temps')
plt.ylabel('Valeurs')
plt.title('Évolution des variables au fil du temps')
plt.legend()
plt.grid(True)
plt.show()



print("A:", A)
print("C:", C)
print("S:", S)
print("R:", R)

print("E:", E)
print("Psi:", Psi)  
print("v:", v)
