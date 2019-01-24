'''
Beam Python Code with Member Loads
Units: Kips (k) and Inches (in)
'''
from __future__ import division #fixes floating point and integer division
######## IMPORT MODULES ########
import numpy as np # import matrix operators
######## FORMAT FLOATS ########
float_formatter = lambda x: "%.2f" % x
np.set_printoptions(formatter={'float_kind':float_formatter})


######## DEFINE FUNCTIONS ########

solve = np.linalg.solve               # set up a shorthand for factorization
###
def Create_k(E,I,L):
    k = (E*I/L**3)*np.array([[12    ,6*L   ,-12   ,6*L   ],
                             [6*L   ,4*L**2,-6*L  ,2*L**2],
                             [-12   ,-6*L  ,12    ,-6*L  ],
                             [6*L   ,2*L**2,-6*L  ,4*L**2]])
    return k
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
NEdof=4

# Define the number of Structural dofs
NSdof=3

# Define {P} vector
P = np.array([[0],
              [0],
              [0]])
print 'P ='
print P
print

# Define the code numbers
CodeNum=np.array([[4,5,1,2],
                  [1,2,6,3]]) # Here we write out each row explicitly
print 'CodeNum ='
print CodeNum
print 

# Initialize S
S = np.zeros((NSdof,NSdof),dtype=float)
Pf = np.zeros((NSdof,1),dtype=float)
print "Initialize S:"
print S
# Define geometry and loading for each member. Create the member stiffness
# and then assemble [S]. Be consistent with units.

# Modulus of Elasticity
E = 4000 #ksi

MemNum=1
I=1500 
L=16*12 
W = 36 
Qf1 = np.array([[W/2],
                [W*L/8],
                [W/2],
                [-W*L/8]])
k1 = Create_k(E,I,L)
S, Pf = Assemble_S_Pf(S,Pf,k1,Qf1,MemNum,NEdof,NSdof,CodeNum)
print 'Qf1 ='
print Qf1
print
print 'k1 ='
print k1
print
print "S = "
print S
MemNum=2
I=1500 
L=8*12 
a = 2*12 
b = 6*12 
M = 8*12 
Qf2 = np.array([[-6*M*a*b/L**3],
                [M*b*(b-2*a)/L**2],
                [6*M*a*b/L**3],
                [M*a*(a-2*b)/L**2]])
k2 = Create_k(E,I,L)
S, Pf = Assemble_S_Pf(S,Pf,k2,Qf2,MemNum,NEdof,NSdof,CodeNum)
print 'Qf2 ='
print Qf2
print
print 'k2 ='
print k2
print
print "S = "
print S

print 'Pf='
print Pf
print
print 'S ='
print S
print


######## Solve ########

d = solve(S,(P-Pf))   
print 'd ='
print d
print

######## POST-PROCESS ########

# Use compatibility to create the {u} vectors for each element
u1 = np.array([[0],
               [0],
               d[0],
               d[1]])   # Indexing starts with zero (0) in Python

u2 = np.array([d[0],
               d[1],
               [0],
               d[2]])


# Calculate the member end forces
Q1 = Qf1 + k1.dot(u1)             
Q2 = Qf2 + k2.dot(u2)

print 'Q1 ='
print Q1
print
print 'Q2 ='
print Q2
print

# Calculate the support reactions
R4 = Q1[0]
R5 = Q1[1]
R6 = Q2[2]

print 'R4 =', R4
print 'R5 =', R5
print 'R6 =', R6


