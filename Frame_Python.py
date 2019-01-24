'''
2D Frame Python Code with Member Loads
Units: kN and mm
'''
from __future__ import division #fixes floating point and integer division 
######## IMPORT MODULES ########
import numpy as np # import matrix operators
######## FORMAT FLOATS ########
float_formatter = lambda x: "%.5f" % x
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
NEdof = 6

# Define the number of Structural dofs
NSdof = 2

# Initialize {P} and {S} 
P = np.array([[1200],
             [ 0 ]])

S = np.zeros((NSdof,NSdof),dtype=float)
Pf = np.zeros((NSdof,1),dtype=float)

E = 4500 #ksi

print 'P ='
print P
print

# Define the code numbers
CodeNum=np.array([[3,4,5,6,7,1],
                  [6,7,1,8,9,2],
                  [8,9,2,10,11,12]]) # Here we write out each row explicitly
print 'CodeNum ='
print CodeNum
print 

# Define geometry and loading for each member, create the member stiffness and fixed-end-force vectors
# and then assemble {Pf} and [S]. Note that no units are shown, you just need to
# be consistent.

#==============================================================================
# MEMBER [1]
#==============================================================================

MemNum = 1
M = 3600    #kip-in
a = b = 120 #in

A = 80      #in^2
I = 600     #in^4
LX = 0      # sign (+/-) important
LY = -240   #in  # sign (+/-) important
L = 240     #in

FSb =  (6*M*a*b/L**3)
FMb = -(M*b/L**2*(b-2*a))
FSe = -(6*M*a*b/L**3)
FMe = -(M*a/L**2*(a-2*b))


Qf1 = np.array([[0],
                [FSb],
                [FMb],
                [0],
                [FSe],
                [FMe]])

 
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

#==============================================================================
# MEMBER [2]
#==============================================================================

MemNum = 2
A = 108    #in^2
I = 1300   #in^4
LX = -300  #in sign (+/-) important
LY = 0     #in sign (+/-) important

L = np.sqrt(LX**2+LY**2)
w = 2/12 #kip/in

FAb = 0
FSb = -w*L/2
FMb = -w*L**2/12
FAe = 0
FSe = -w*L/2
FMe = w*L**2/12

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

#==============================================================================
# MEMBER [3]
#==============================================================================

MemNum = 3
W = 20    #kip
a = 144   #in
b = 96    #in

A = 80    #in^2
I = 600   #in^4
LX = 0    # sign (+/-) important
LY = -240 #in  # sign (+/-) important
L = 240   #in

FAb = 0
FSb = (W*b**2/L**3)*(3*a+b)
FMb = W*a*b**2/L**2
FAe = 0
FSe = (W*a**2/L**3)*(3*b+a)
FMe = -W*b*a**2/L**2

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
               [0],
               [0],
               [0],
              d[0]])   # Indexing starts with zero (0) in Python

v2 = np.array([[0],
               [0],
              d[0],
               [0],
               [0],
              d[1]])

v3 = np.array([[0],
               [0],
              d[1],
               [0],
               [0],
               [0]])

# Calculate the member end forces (global)
F1 = Ff1 + K1.dot(v1)
F2 = Ff2 + K2.dot(v2)
F3 = Ff3 + K3.dot(v3)

print 'F1 ='
print F1
print
print 'F2 ='
print F2
print
print 'F3 ='
print F3
print
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
R3 =  F1[0]
R4 =  F1[1]
R5 =  F1[2]
R6 =  F1[3] + F2[0]
R7 =  F1[4] + F2[1]
R8 =  F2[3] + F3[0]
R9 =  F2[4] + F3[1]
R10 = F3[3]
R11 = F3[4]
R12 = F3[5]

print 'R3 =', R3[0], ' kips'
print 'R4 =', R4[0], ' kips'
print 'R5 =', R5[0], ' kip-in'
print 'R6 =', R6[0], ' kips'
print 'R7 =', R7[0], ' kips'
print 'R8 =', R8[0], ' kips'
print 'R9 =', R9[0], ' kips'
print 'R10 =', R10[0], ' kips'
print 'R11 =', R11[0], ' kips'
print 'R12 =', R12[0], ' kip-in'
