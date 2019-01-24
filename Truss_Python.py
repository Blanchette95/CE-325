'''
2D Truss MATLAB Code
units: kN and mm
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
def Create_TK(LX, LY, E, A):
    L = np.sqrt(LX**2 + LY**2)
    c = LX/L #cosine
    s = LY/L #sine
    k = (E*A/L)*np.array([[1, 0,-1, 0],
                          [0, 0, 0, 0],
                          [-1,0, 1, 0],
                          [0, 0, 0, 0]])
    
    T = np.array([[c, s, 0, 0],
                  [-s,c, 0, 0],
                  [0, 0, c, s],
                  [0, 0,-s, c]])
    
    K = np.transpose(T).dot(k).dot(T)
    return T, K
###
def Assemble_S(S, k, MemNum, NEdof, NSdof, CodeNum):
    for KRow in range(NEdof):
        SRow = CodeNum[MemNum-1, KRow]
        if SRow <= NSdof:
            for KCol in range(NEdof):
                SCol = CodeNum[MemNum-1, KCol]
                if SCol <= NSdof:
                    S[SRow-1, SCol-1] = S[SRow-1, SCol-1] + k[KRow, KCol]
                    
    return S
###

######## DEFINE SYSTEM FOR ANALYSIS ########

# Define the number of element dofs (always 4 for the 2D truss element)
NEdof=4

# Define the number of Structural dofs
NSdof=2

# Initialize S
S = np.zeros((NSdof, NSdof),dtype=float)

#Modulus of Elasticity
E = 29000 # ksi

# Define {P} vector
P = np.array([[-100*np.cos(np.deg2rad(60))], [-100*np.sin(np.deg2rad(60))]])
print 'P = '
print P
print
    
# Define the code numbers
CodeNum=np.array([[3,4,1,2],
                  [5,6,1,2],
                  [7,8,1,2]]) # Here we write out each row explicitly
print 'CodeNum = ' 
print CodeNum
print




# Define geometry and loading for each member. Create the member stiffness
# and then assemble [S]. Note that no units are shown, you just need to be
# be consistent.

MemNum=1
A=4.0 #in^2
LX=120.0 # sign (+/-) important
LY=0   # sign (+/-) important 
T1, K1 = Create_TK(LX, LY, E, A)
S = Assemble_S(S, K1, MemNum, NEdof, NSdof, CodeNum)
print 'T1 ='
print T1
print
print 'K1 ='
print K1
print

MemNum=2
A=2.0
LX=0 # sign (+/-) important
LY=144.0  # sign (+/-) important
T2, K2 = Create_TK(LX, LY, E, A)
S = Assemble_S(S, K2, MemNum, NEdof, NSdof, CodeNum)
print 'T2 ='
print T2
print
print 'K2 ='
print K2
print

MemNum=3
A=6.0
LX=-9.0*12 # sign (+/-) important
LY=144.0 # sign (+/-) important
T3, K3 = Create_TK(LX, LY, E, A)
S = Assemble_S(S, K3, MemNum, NEdof, NSdof, CodeNum)
print 'T3 ='
print T3
print
print 'K3 ='
print K3
print

print 'S ='
print S
print


######## Solve ########

d = solve(S, P)

print 'd ='
print d
print

######## POST-PROCESS ########

# Use compatibility to create the {v} vectors for each element
v1 = np.array([ [0], [0], d[0], d[1]])
v2 = np.array([ [0], [0], d[0], d[1]])
v3 = np.array([ [0], [0], d[0], d[1]])

# Calculate the member end forces
F1 = K1.dot(v1)       # GLOBAL {F} = [K]{v}      
F2 = K2.dot(v2)  
F3 = K3.dot(v3)

print 'F1 ='
print F1
print
print 'F2 ='
print F2
print
print 'F3 ='
print  F3
print

Q1 = T1.dot(F1)       # LOCAL {Q} = [T]{F}      
Q2 = T2.dot(F2)  
Q3 = T3.dot(F3)

print 'Q1 ='
print Q1
print
print 'Q2 ='
print Q2
print
print 'Q3 ='
print  Q3
print

#Calculate the bar (axial) forces
BarForces = np.array([[Q1[2]],[Q2[2]],Q3[2]])

# Calculate the support reactions
R3 = F1[0]
R4 = F1[1]
R5 = F2[0] 
R6 = F2[1]
R7 = F3[0]
R8 = F3[1]

print "R3 = ",R3[0]
print
print 'R4 =', R4[0]
print
print 'R5 =', R5[0]
print
print 'R6 =', R6[0]
print
print "R7 =", R7[0]
print 
print "R8 =", R8[0]




