'''
2D Frame Python Code with Member Loads
Units: kN and mm
'''
from __future__ import division #fixes floating point and integer division 
######## CLEAR VARS ########
from IPython import get_ipython
get_ipython().magic('reset -sf') 
######## IMPORT MODULES ########
import numpy as np # import matrix operators
######## FORMAT FLOATS ########
float_formatter = lambda x: "%.4f" % x
np.set_printoptions(formatter={'float_kind':float_formatter})


######## DEFINE FUNCTIONS ########

solve = np.linalg.solve               # set up a shorthand for factorization
###
def Create_TKF(E,I,A,LX,LY,Qf):
    L=np.sqrt(LX**2+LY**2)
    c=LX/L #cosine
    s=LY/L #sine
    
    k = (E*I)/(L**3)*np.array([[A*L**2/I  ,0            ,0         ,-A*L**2/I ,0            ,0         ],
                               [0         ,12           ,6*L       ,0         ,-12          ,6*L       ],
                               [0         ,6*L          ,4*L**2    ,0         ,-6*L         ,2*L**2    ],
                               [-A*L**2/I ,0            ,0         ,A*L**2/I  ,0            ,0         ],
                               [0         ,-12          ,-6*L      ,0         ,12           ,-6*L      ],
                               [0         ,6*L          ,2*L**2    ,0         ,-6*L         ,4*L**2    ]])

    T = np.array([[c ,s ,0 ,0 ,0 ,0 ],
                  [-s,c ,0 ,0 ,0 ,0 ],
                  [0 ,0 ,1 ,0 ,0 ,0 ],
                  [0 ,0 ,0 ,c ,s ,0 ],
                  [0 ,0 ,0 ,-s,c ,0 ],
                  [0 ,0 ,0 ,0 ,0 ,1 ]])

    Ff = np.transpose(T).dot(Qf)

    K = np.transpose(T).dot(k).dot(T)
  
    return K,Ff,T
###
def Assemble_S_Pf(S,Pf,k,Qf,MemNum,NEdof,NSdof,CodeNum):
    for KRow in range(0,NEdof,1):
        SRow=CodeNum[MemNum-1,KRow]
        if SRow <= NSdof:
            Pf[SRow-1]=Pf[SRow-1]+Qf[KRow]
            for KCol in range(0,NEdof,1):
                SCol=CodeNum[MemNum-1,KCol]
                if SCol <= NSdof:
                    S[SRow-1,SCol-1]=S[SRow-1,SCol-1]+k[KRow,KCol]
    return S, Pf
###

######## DEFINE SYSTEM FOR ANALYSIS ########

# Define the number of element dofs (always 2 for the uniaxial element)
NEdof=6

# Define the number of Structural dofs
NSdof=7

# Initialize {P} and {S} 
P = np.zeros((NSdof,1),dtype=float)
S = np.zeros((NSdof,NSdof),dtype=float)
Pf = np.zeros((NSdof,1),dtype=float)

E = 4500 #ksi

#Define nonzero terms in P
P[0,0] = 0

print 'P ='
print P
print

# Define the code numbers
CodeNum=np.array([[8,9,1,2,3,4],
                  [2,3,4,5,6,7],
                  [5,6,7,10,11,12]]) # Here we write out each row explicitly
print 'CodeNum ='
print CodeNum
print 

# Define geometry and loading for each member, create the member stiffness and fixed-end-force vectors
# and then assemble {Pf} and [S]. Note that no units are shown, you just need to
# be consistent.

MemNum=1
Qf1 = np.zeros((NEdof,1),dtype=float)
A = 40
I = 300 #in^4
LX = 0    # sign (+/-) important
LY = 120  # sign (+/-) important
 
print 'Qf1 ='
print Qf1
print 
K1,Ff1,T1 = Create_TKF(E,I,A,LX,LY,Qf1) 
print 'K1 ='
print K1
print
print 'Ff1 ='
print Ff1
print 
print 'T1 ='
print T1
print
S, Pf = Assemble_S_Pf(S,Pf,K1,Ff1,MemNum,NEdof,NSdof,CodeNum)

MemNum=2
#Qf2 = np.zeros((NEdof,1),dtype=float)
A = 60
I = 600
LX = 144  # sign (+/-) important
LY = 60  # sign (+/-) important

L = np.sqrt(LX**2+LY**2)

W=.0417
WX = (W*LY/L)
WY = (W*LX/L)

FAb=WX*L/2
FSb=WY*L/2
FMb=WY*L*L/12
FAe=WX*L/2
FSe=WY*L/2
FMe=-WY*L*L/12

Qf2 = np.array([[FAb],
                [FSb],
                [FMb],
                [FAe],
                [FSe],
                [FMe]])

 
print 'Qf2 ='
print Qf2
print 
K2,Ff2,T2 = Create_TKF(E,I,A,LX,LY,Qf2)
print 'K2 ='
print K2
print
print 'Ff2 ='
print Ff2
print 
print 'T2 ='
print T2
print  
S, Pf = Assemble_S_Pf(S,Pf,K2,Ff2,MemNum,NEdof,NSdof,CodeNum)

MemNum=3
Qf3 = np.zeros((NEdof,1),dtype=float)
A = 60
I = 600
LX = 96  # sign (+/-) important
LY = 0 # sign (+/-) important

L = np.sqrt(LX**2+LY**2)

W=.0417

FAb=0
FSb=W*L/2
FMb=W*L*L/12
FAe=0
FSe=W*L/2
FMe=-W*L*L/12

Qf3 = np.array([[FAb],
                [FSb],
                [FMb],
                [FAe],
                [FSe],
                [FMe]])

 
print 'Qf3 ='
print Qf3
print 
K3,Ff3,T3 = Create_TKF(E,I,A,LX,LY,Qf3)
print 'K3 ='
print K3
print
print 'Ff3 ='
print Ff3
print 
print 'T3 ='
print T3
print  
S, Pf = Assemble_S_Pf(S,Pf,K3,Ff3,MemNum,NEdof,NSdof,CodeNum)


######## Solve ########
print 'P ='
print P
print
print 'Pf ='
print Pf
print
print 'S ='
print S
print


d = solve(S,(P-Pf))   
print 'd ='
print d
print

######## POST-PROCESS ########

# Use compatibility to create the {v} vectors for each element
v1 = np.array([[0],
               [0],
               d[0],
               d[1],
               d[2],
               d[3]])   # Indexing starts with zero (0) in Python

v2 = np.array([d[1],
               d[2],
               d[3],
               d[4],
               d[5],
               d[6]])

v3 = np.array([d[4],
               d[5],
               d[6],
                [0],
                [0],
                [0]])

# Calculate the member end forces (global)
F1=Ff1+K1.dot(v1)
F2=Ff2+K2.dot(v2)
F3=Ff3+K3.dot(v3)

# Calculate the member end forces (local)
Q1 = T1.dot(F1)            
Q2 = T2.dot(F2)
Q3 = T3.dot(F3)

print 'Q1 ='
print Q1
print
print 'Q2 ='
print Q2
print
print 'Q3 ='
print Q3
print

# Calculate the support reactions
R8 = F1[0]
R9 = F1[1]
R10 = F3[3]
R11 = F3[4]
R12 = F3[5]

print 'R8 =', R8[0]
print 'R9 =', R9[0]
print 'R10 =', R10[0]
print 'R11 =', R11[0]
print 'R12 =', R12[0]
