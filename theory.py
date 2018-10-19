from math import *
from scipy.integrate import *
import numpy as np
#CONSTANT TIME PERCENTAGE INFECTED
#VARYING LAMBDA
'''
rho = 4
mu = 5
writefile = open('theory/constrho4.csv','w')
writefile.write('lambda,rho,t\n')
for biglambd in range(1,100):
        lambd = biglambd/10
        tau = 1 - exp(-lambd)
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        alpha = 1 - taur/tau
        A = tau*(alpha+mu*(1-alpha))
        B = tau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#CONSTANT TIME RHO=4 LOWER BOUND
writefile = open('theory/constrho4lower.csv','w')
writefile.write('lambda,rho,z\n')
rho = 4
mu = 5
for biglambd in range(1,100):
        lambd = biglambd/10
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        z = 0
        for i in range(1000):
                z = exp(-mu*taur*(1-z))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#VARYING RHO
lambd = 1
mu = 5
writefile = open('theory/constlambda1.csv','w')
writefile.write('lambda,rho,t\n')
for bigrho in range(1,80):
        rho = bigrho/20
        tau = 1 - exp(-lambd)
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        alpha = 1 - taur/tau
        A = tau*(alpha+mu*(1-alpha))
        B = tau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#CONSTANT TIME LAMBDA=1 LOWER BOUND
writefile = open('theory/constlambda1lower.csv','w')
writefile.write('lambda,rho,z\n')
lambd = 1
mu = 5
for bigrho in range(1,80):
        rho = bigrho/20
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        z = 0
        for i in range(1000):
                z = exp(-mu*taur*(1-z))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#EXPONENTIAL TIME SURVIVAL PROBABILITY
#VARYING LAMBDA
writefile = open('theory/exprho4prob.csv','w')
rho = 4
mu = 5
def integrand(t):
        return exp(-t)*exp(mur*(1-z)*exp(-(lambd+rho)*t))
writefile.write('lambda,rho,u\n')
for biglambd in range(1,100):
        lambd = biglambd/10
        mur = mu*lambd/(lambd+rho)
        z = .5
        for i in range(100):
                z = exp(-mur*(1-z))*quad(integrand,0,np.inf)[0]
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#VARYING RHO
writefile = open('theory/explambda1prob.csv','w')
lambd = 1
mu = 5
def integrand(t):
        return exp(-t)*exp(mur*(1-z)*exp(-(lambd+rho)*t))
writefile.write('lambda,rho,u\n')
for bigrho in range(1,80):
        rho = bigrho/20
        mur = mu*lambd/(lambd+rho)
        z = .5
        for i in range(1000):
                z = exp(-mur*(1-z))*quad(integrand,0,np.inf)[0]
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#EXPONENTIAL TIME PERCENTAGE INFECTED
#VARYING LAMBDA
rho = 4
mu = 5
writefile = open('theory/exprho4size.csv','w')
writefile.write('lambda,rho,t\n')
for biglambd in range(1,200):
        lambd = biglambd/20
        Etau = lambd/(lambd+1)#
        taur = lambd/(lambd+1+rho)
        #is the taur calculation the same?
        alpha = 1 - taur/Etau
        A = Etau*(alpha+mu*(1-alpha))
        B = Etau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#VARYING RHO
lambd = 1
mu = 5
writefile = open('theory/explambda1size.csv','w')
writefile.write('lambda,rho,t\n')
for bigrho in range(1,80):
        rho = bigrho/10
        Etau = lambd/(lambd+1)
        taur = lambd/(lambd+1+rho)
        alpha = 1 - taur/Etau
        A = Etau*(alpha+mu*(1-alpha))
        B = Etau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#NO REWIRING CONSTANT
rho = 0
mu = 5
writefile = open('theory/constrho0.csv','w')
writefile.write('lambda,rho,t\n')
for biglambd in range(1,50):
        lambd = biglambd/10
        tau = 1 - exp(-lambd)
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        alpha = 1 - taur/tau
        A = tau*(alpha+mu*(1-alpha))
        B = tau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#NO REWIRING EXPONENTIAL
rho = 0
mu = 5
writefile = open('theory/exprho0size.csv','w')
writefile.write('lambda,rho,t\n')
for biglambd in range(1,50):
        lambd = biglambd/10
        Etau = lambd/(lambd+1)
        taur = lambd/(lambd+1+rho)
        alpha = 1 - taur/Etau
        A = Etau*(alpha+mu*(1-alpha))
        B = Etau*alpha
        t = 1
        for i in range(1000):
                t = 1 - A/(B+(A-B)*exp(A*t))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(t)+ '\n')
writefile.close()
#NO REWIRING CONSTANT LOWER BOUND
writefile = open('theory/constrho0lower.csv','w')
writefile.write('lambda,rho,z\n')
rho = 0
mu = 5
for biglambd in range(1,100):
        lambd = biglambd/10
        taur = (lambd/(lambd+rho))*(1-exp(-(lambd+rho)))
        z = 0
        for i in range(1000):
                z = exp(-mu*taur*(1-z))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#NO REWIRING EXPONENTIAL LOWER BOUND
writefile = open('theory/exprho0prob.csv','w')
rho = 0
mu = 5
def integrand(t):
        return exp(-t)*exp(mur*(1-z)*exp(-(lambd+rho)*t))
writefile.write('lambda,rho,u\n')
for biglambd in range(1,100):
        lambd = biglambd/10
        mur = mu*lambd/(lambd+rho)
        z = .5
        for i in range(100):
                z = exp(-mur*(1-z))*quad(integrand,0,np.inf)[0]
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#NO REWIRING EXPONENTIAL GENERATING FUNCTION
writefile = open('theory/exprho0genfunc.csv','w')
rho = 0
mu = 5
#def integrand(t):
#        return exp(-t)*exp(mur*(1-z)*exp(-
#        (lambd+rho)*t))
writefile.write('lambda,rho,z\n')
for biglambd in range(1,100):
        lambd = biglambd/10
        mur = mu*lambd/(lambd+rho)
        Etau = lambd/(lambd+1)
        z = .5
        for i in range(100):
                z = exp(-mu*Etau*(1-z))
                #z = exp(-mur*(1-z))*quad(integrand,0,np.inf)[0]
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#RHO=4 EXPONENTIAL GENERATING FUNCTION
writefile = open('theory/exprho4genfunc.csv','w')
rho = 4
mu = 5
writefile.write('lambda,rho,z\n')
for biglambd in range(1,200):
        lambd = biglambd/20
        Etau = lambd/(lambd+1)
        z = .5
        for i in range(1000):#Make this the GF for all exponential sizes
                z = exp(-mu*(lambd/(lambd+1+rho))*(1-z))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()
#LAMBDA=1 EXPONENTIAL GENERATING FUNCTION
writefile = open('theory/explambda1genfunc.csv','w')
lambd = 1
mu = 5
writefile.write('lambda,rho,z\n')
for bigrho in range(1,160):
        rho = bigrho/40
        Etau = lambd/(lambd+1)
        z = .5
        for i in range(100):#Make this the GF for all exponential sizes
                z = exp(-mu*(lambd/(lambd+1+rho))*(1-z))
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(z)+ '\n')
writefile.close()

#NO REWIRE REAL EXP TIME THEORY
writefile = open('theory/realtimeexptheory.csv','w')
writefile.write('S,I,R,psi,t\n')
mu = 5
lambd = 1
rho = 0
n = 10000
Etau = lambd/(lambd+rho+1)
divs = 10000
U = (n-mu-1)
E = mu
I = 1
R = 0
numits = 20*divs
for i in range(numits):
        U += -lambd*E*mu*Etau*U/n/divs
        E += lambd*E*(mu*Etau*U/n-1)/divs
        I += (lambd*E - I)/divs
        R += I/divs
        S = E+U
        if((10*i) % divs == 0):
                writefile.write(str(S/n) + ',' + str(I/n) + ',' + str(R/n) + ',' + str(E) + ',' + str(i/divs) + '\n')
writefile.close()
#NO REWIRE IMPROVED ODE
writefile = open('theory/realtimeexptheoryimproved.csv','w')
writefile.write('S,I,R,t\n')
mu = 5
lambd = 1
rho = 0
n = 10000
Etau = lambd/(lambd+rho+1)
J = 10
divs = 100
I = 1
R = 0
E = [0]*(J+1)
E[0] = n-mu-1
E[1] = mu
F = 0
for j in range(len(E)):
        F += j*E[j]
numits = 20*divs
for i in range(numits):
        F = 0
        for j in range(len(E)):
                F += j*E[j]
        E[0] += (-lambd*F*mu*Etau*E[0]/n)/divs
        for j in range(1,J):
                E[j] += (-lambd*j*E[j] + lambd*F*(mu*Etau*E[j-1]/n-mu*Etau*E[j]/n))/divs
        E[J] += (-lambd*J*E[J]+lambd*F*mu*Etau*E[J-1]/n)/divs
        I += (lambd*F-I)/divs
        R += I/divs
        S = sum(E)
        if((10*i) % divs == 0):
                writefile.write(str(S/n) + ',' + str(I/n) + ',' + str(R/n) +  ',' + str(i/divs) + '\n')
writefile.close()
#LAMBDA=1 IMPROVED ODE
writefile = open('theory/exprho4improvedode.csv','w')
rho = 4
mu = 5
n = 10000
writefile.write('lambda,rho,R\n')
for biglambd in range(1,200):
        lambd = biglambd/20
        Etau = lambd/(lambd+rho+1)
        J = 10
        divs = 1000
        I = 1
        R = 0
        E = [0]*(J+1)
        E[0] = n-mu-1
        E[1] = mu
        F = 0
        for j in range(len(E)):
                F += j*E[j]
        numits = 40*divs
        for i in range(numits):
                F = 0
                for j in range(len(E)):
                        F += j*E[j]
                E[0] += (-lambd*F*mu*Etau*E[0]/n)/divs
                for j in range(1,J):
                        E[j] += (-lambd*j*E[j] + lambd*F*(mu*Etau*E[j-1]/n-mu*Etau*E[j]/n))/divs
                E[J] += (-lambd*J*E[J]+lambd*F*mu*Etau*E[J-1]/n)/divs
                I += (lambd*F-I)/divs
                R += I/divs
        S = sum(E)
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(R/n)+ '\n')
writefile.close()
#Miller11 model
writefile = open('theory/miller11.csv','w')
writefile.write('S,I,R,theta,t\n')
mu = 5
lambd = 1
rho = 0
n = 10000
Etau = lambd/(lambd+rho+1)
divs = 10000
S = (n-1)/n
I = 1/n
R = 0/n
beta = lambd
theta = .999999
gamma = 1
def psi(z):
        return exp(-mu*(1-z))
numits = 20*divs
for i in range(numits):
        theta += (gamma - (beta + gamma)*theta + beta*exp(-mu*(1-theta)))/divs
        S = psi(theta)
        I = 1 - R - S
        R += gamma*I/divs
        if((10*i) % divs == 0):
                writefile.write(str(S) + ',' + str(I) + ',' + str(R) + ',' + str(theta) + ',' + str(i/divs) + '\n')
writefile.close()
#third ode
writefile = open('theory/ODE3.csv','w')
writefile.write('S,I,R,t\n')
mu = 5
lambd = 1
rho = 0
n = 10000
Etau = lambd/(lambd+rho+1)
J = 20
divs = 100
I = 1
R = 0
E = [0]*(J+1)
E[0] = n-mu-1
E[1] = mu
F = 0
for j in range(len(E)):
        F += j*E[j]
numits = 20*divs
for i in range(numits):
        F = 0
        for j in range(len(E)):
                F += j*E[j]
        E[0] += (-lambd*F*mu/n*E[0]+E[1])/divs
        for j in range(1,J):
                E[j] += (-lambd*j*E[j]+lambd*F*mu/n*(E[j-1]-E[j])+(j+1)*E[j+1]-j*E[j])/divs
        E[J] += (-lambd*J*E[J]+lambd*F*mu/n*(E[J-1]-E[J])-J*E[J])/divs
        I += (lambd*F-I)/divs
        R += I/divs
        S = sum(E)
        if((10*i) % divs == 0):
                writefile.write(str(S/n) + ',' + str(I/n) + ',' + str(R/n) +  ',' + str(i/divs) + '\n')
writefile.close()
#Constant time star graph
writefile = open('theory/consttimestar.csv','w')
writefile.write('lambda,rho,prop\n')
def firststep(lambd,rho,r):
        return lambd*exp(-(lambd+rho)*r)
def secondstep(lambd,rho,r1,r2,r3):
        return lambd*exp(-(lambd+rho)*r1)*rho*exp(-(lambd+rho)*r2)*lambd*exp(-(lambd+rho)*r3)
totalprop = 0
rho = 4
for biglambd in range(100):
        totalprop = 0
        lambd = biglambd/10
        totalprop += quad(lambda r:firststep(lambd,rho,r),0,1)[0]
        totalprop += tplquad(lambda r3,r2,r1:secondstep(lambd,rho,r1,r2,r3),0,1,0,lambda r1: r1,0,1)[0]
        #totalprop += tplquad(lambda r3,r2,r1:secondstep(lambd,rho,r1,r2,r3),0,1,lambda r1: r1,lambda r1: r1,0,lambda r1,r2:r1+1-r2)[0]
        totalprop += tplquad(lambda r3,r2,r1:secondstep(lambd,rho,r1,r2,r3),0,1,lambda r1: r1,1,0,lambda r1,r2:r1+1-r2)[0]
        writefile.write(str(lambd)+','+str(rho)+','+str(totalprop)+'\n')
lambd = 1
for bigrho in range(100):
        totalprop = 0
        rho = bigrho/10
        totalprop += quad(lambda r:firststep(lambd,rho,r),0,1)[0]
        totalprop += tplquad(lambda r3,r2,r1:secondstep(lambd,rho,r1,r2,r3),0,1,0,lambda r1: r1,0,1)[0]
        totalprop += tplquad(lambda r3,r2,r1:secondstep(lambd,rho,r1,r2,r3),0,1,lambda r1: r1,lambda r1: r1+1,0,lambda r1,r2:r1+1-r2)[0]
        writefile.write(str(lambd)+','+str(rho)+','+str(totalprop)+'\n')
writefile.close()
'''
#real time ode with rewiring
writefile = open('theory/realtimeoderewire.csv','w')
writefile.write('S,I,R,mu,t\n')
mu = 5
lambd = 2
rho = 3
n = 100000
J = 30
divs = 100
I = 1
R = 0
S = [0]*J
S[0] = float(n-I-R-mu)
S[1] = float(mu)
dS = [0]*len(S)
numits = 20*divs
for i in range(numits):
        t = i/divs
        F = 0
        for i in range(len(S)):
                F += i*S[i]
        dS[0] =  lambd*F*(mu/n)*(-S[0]) + S[1] + rho*(1-I/n)*S[1] 
        for k in range(1,len(S)-1):
                dS[k] = -lambd*k*S[k] + lambd*F*(mu/n)*(S[k-1]-S[k]) + (k+1)*S[k+1] - k*S[k] + rho*(1-I/n)*((k+1)*S[k+1]-k*S[k]) 
        dmu = rho*F/n*(1 - (I+R)/n)
        dI = lambd*F-I
        dR = I
        for k in range(len(S)):
                S[k] += dS[k]/divs
        mu += dmu/divs
        I += dI/divs
        R += dR/divs
        writefile.write(str(sum(S)/n)+','+str(I/n)+','+str(R/n)+','+str(mu)+','+str(t)+'\n')
writefile.close()
#real time ode with rewiring
writefile = open('theory/ode3withrewiring.csv','w')
writefile.write('lambda,rho,S,I,R,mu\n')
n = 1000
J = 30
divs = 1000
numits = 20*divs
rho = 4
for biglambd in range(0,100):
        lambd = biglambd/10
        mu = 5
        I = 1
        R = 0
        S = [0]*J
        S[0] = float(n-I-mu)
        S[1] = float(mu)
        dS = [0]*len(S)
        for i in range(numits):
                t = i/divs
                F = 0
                for i in range(len(S)):
                        F += i*S[i]
                dS[0] =  lambd*F*(mu/n)*(-S[0]) + S[1] + rho*(1-I/n)*S[1] 
                for k in range(1,len(S)-1):
                        dS[k] = -lambd*k*S[k] + lambd*F*(mu/n)*(S[k-1]-S[k]) + (k+1)*S[k+1] - k*S[k] + rho*(1-I/n)*((k+1)*S[k+1]-k*S[k]) 
                dmu = rho*F/n*(1 - (I+R)/n)
                dI = lambd*F-I
                dR = I
                for k in range(len(S)):
                        S[k] += dS[k]/divs
                mu += dmu/divs
                I += dI/divs
                R += dR/divs
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(sum(S)/n)+','+str(I/n)+','+str(R/n)+','+str(mu)+'\n')
lambd = 1
for bigrho in range(0,100):
        rho = bigrho/10
        mu = 5
        I = 1
        R = 0
        S = [0]*J
        S[0] = float(n-I-R-mu)
        S[1] = float(mu)
        dS = [0]*len(S)
        for i in range(numits):
                t = i/divs
                F = 0
                for i in range(len(S)):
                        F += i*S[i]
                dS[0] =  lambd*F*(mu/n)*(-S[0]) + S[1] + rho*(1-I/n)*S[1] 
                for k in range(1,len(S)-1):
                        dS[k] = -lambd*k*S[k] + lambd*F*(mu/n)*(S[k-1]-S[k]) + (k+1)*S[k+1] - k*S[k] + rho*(1-I/n)*((k+1)*S[k+1]-k*S[k]) 
                dmu = rho*F/n*(1 - (I+R)/n)
                dI = lambd*F-I
                dR = I
                for k in range(len(S)):
                        S[k] += dS[k]/divs
                mu += dmu/divs
                I += dI/divs
                R += dR/divs
        writefile.write(str(lambd) + ',' + str(rho) + ',' + str(sum(S)/n)+','+str(I/n)+','+str(R/n)+','+str(mu)+'\n')
writefile.close()
