import numpy as np
import matplotlib.pyplot as plt
import math
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
Lambda=0.8
Rm=7
mE=0.01
mLambda=0.001
alfa = 0.5
beta=0.5
gama=10

# initialisation de C,S et E
C1 = np.zeros(300)  
C2= np.zeros(300)
S1 = np.zeros(300)
S2= np.zeros(300)  
E = np.zeros(300)  
A1 = np.zeros(300)
A2 = np.zeros(300)  
v1 = np.zeros(300) 
v2 = np.zeros(300) 
Psi1 = np.zeros(300) 
Psi2 = np.zeros(300) 

S1[0]=5
E[0]=E0
S2[0]=-5
C2[0]=5
C1[0]=-5




# itérations 
for i in range(300):
    R = min(np.random.poisson(Lambda), Rm)
    Psi1[i]= C1[i] -E[i] - S1[i]
    v1[i]=max(0,min(Psi1[i],1))
    A1[i]=q*v1[i] + (R/Rm)*q*(1-v1[i])

    Psi2[i]= C2[i] -E[i] - S2[i]
    v2[i]=max(0,min(Psi2[i],1))
    A2[i]=q*v2[i] + (R/Rm)*q*(1-v2[i])
    if i<299:
        C1[i+1]=(1-d)*C1[i]+alfa*A2[i]*C1[i]+b*min(1,1-C1[i])*A1[i]
        S1[i+1]=S1[i]+(p+beta*math.exp(-gama*A2[i]))*max(0,Sm-S1[i])-h*C1[i]-k*A1[i]
        C2[i+1]=(1-d)*C2[i]+alfa*A1[i]*C2[i]+b*min(1,1-C2[i])*A2[i]
        S2[i+1]=S2[i]+(p+beta*math.exp(-gama*A1[i]))*max(0,Sm-S2[i])-h*C2[i]-k*A2[i]
        E[i+1]=E[i]-mE
    Lambda+=0.0001
        



    
    

T=np.linspace(0,300,300)
# Tracer les graphiques
plt.figure(figsize=(10, 6))
#plt.plot(T, A1, label='A1')
plt.plot(T, S1, label='S1')
plt.plot(T, C1, label='C1')
#plt.plot(T, E, label='E')
#plt.plot(T, v1, label='v1')
#plt.plot(T, Psi1, label='Psi1')
#plt.plot(T, A2, label='A2')
plt.plot(T, S2, label='S2')
plt.plot(T, C2, label='C2')
#plt.plot(T, v2, label='v2')
#plt.plot(T, Psi2, label='Psi2')
plt.xlabel('Temps')
plt.ylabel('Valeurs')
plt.title('Évolution des variables au fil du temps')
plt.legend()
plt.grid(True)
plt.show()


"""
print("A:", A)
print("C:", C)
print("S:", S)
print("R:", R)

print("E:", E)
print("Psi:", Psi)  
print("v:", v)
"""
