'''
Beam Python Code with Member Loads
Units: Kips (k) and Inches (in)
'''
from __future__ import division #fixes floating point and integer division 
######## CLEAR VARS ########
from IPython import get_ipython
get_ipython().magic('reset -sf') 

######## IMPORT MODULES ########
import numpy as np # import matrix operators

######## FORMAT FLOATS ########
float_formatter = lambda x: "%8f" % x
np.set_printoptions(precision = 3)#,formatter={'float_kind':float_formatter})
fmt = "[{0:>5.5}] {1:<5}"

######## FORMATTING FUNCTIONS ########
def print_Q(arr):
    units = [' kN\n', ' kN.m\n']
    lst = arr.tolist()
    string = ''
    count = 0
    for num in lst:
        for n in num:
            string += fmt.format(n, units[count%2])
        count += 1
    return string

def print_d(arr):
    units = ' rad\n'
    string = ''
    lst = arr.tolist()
    
    for num in lst:
        for n in num:
            string += fmt.format(n, units)
    return string
######################################
######## DEFINE FUNCTIONS ########
solve = np.linalg.solve               # set up a shorthand for factorization
###
def Create_k(E, I, L):
    k = (E*I/L**3)*np.array([[12   , 6*L   ,-12   , 6*L   ],
                             [6*L  , 4*L**2,-6*L  , 2*L**2],
                             [-12  ,-6*L   , 12   ,-6*L   ],
                             [6*L  , 2*L**2,-6*L  , 4*L**2]])
    return k
###
def Assemble_S_Pf(S, Pf, k, Qf, MemNum, NEdof, NSdof, CodeNum):
    for KRow in range(NEdof):
        SRow = CodeNum[MemNum-1, KRow]
        if SRow <= NSdof:
            Pf[SRow-1] = Pf[SRow-1] + Qf[KRow]
            for KCol in range(NEdof):
                SCol = CodeNum[MemNum-1, KCol]
                if SCol <= NSdof:
                    S[SRow-1,SCol-1] = S[SRow-1, SCol-1] + k[KRow, KCol]
    return S, Pf
###

######## DEFINE SYSTEM FOR ANALYSIS ########
# Define the number of element dofs (always 2 for the uniaxial element)
NEdof = 4

# Define the number of Structural dofs
NSdof = 4

# Define {P} vector
P = np.array([[ 0 ],
              [ 0 ],
              [ 0 ],
              [120.0]])#kN-m
print 'P ='
print P
print

# Define the code numbers
CodeNum = np.array([[5,1,6,2],
                    [6,2,7,3],
                    [7,3,8,4]]) # Here we write out each row explicitly
print 'CodeNum ='
print CodeNum
print 

# Initialize S
S = np.zeros((NSdof, NSdof), dtype = float)
Pf = np.zeros((NSdof, 1), dtype = float)

# Define geometry and loading for each member. Create the member stiffness
# and then assemble [S]. Be consistent with units.

# Modulus of Elasticity
E = 200000000.0 #KPa
I = 0.0004 #m^4 
#### MEMBER [1] ####
print "Assembling matrices for Member [1]"
MemNum = 1
L = 15.0 #m
w = 18.0 #kN/m
Qf1 = np.array([[   w*L/2  ],
                [ w*L**2/12],
                [   w*L/2  ],
                [-w*L**2/12]])

k1 = Create_k(E, I, L)
S, Pf = Assemble_S_Pf(S, Pf, k1, Qf1, MemNum, NEdof, NSdof, CodeNum)
print 'Qf1 ='
print print_Q(Qf1)
print 'k1 ='
print k1
print '\n'
####################
#### MEMBER [2] ####
print "Assembling matrices for Member [2]"
MemNum = 2
W = 90.0 #kN
L = 15.0 #m
l1 = 5.0 #m
l2 = 10.0#m 
 
Qf2 = np.array([[(W/L**3)*(l2**2*(3*l1 + l2) + l1**2*(l1 + 3*l2))],
                [             (W*l1*l2/L**2)*(l1+l2)             ],
                [(W/L**3)*(l1**2*(l1 + 3*l2) + l2**2*(3*l1 + l2))],
                [            -(W*l1*l2/L**2)*(l1+l2)             ]])

k2 = Create_k(E, I, L)
S, Pf = Assemble_S_Pf(S, Pf, k2, Qf2, MemNum, NEdof, NSdof, CodeNum)
print 'Qf2 ='
print print_Q(Qf2)
print 'k2 ='
print k2
print '\n'
####################
#### MEMBER [3] ####
print "Assembling matrices for Member [3]"
MemNum = 3
L = 15.0 #m 
w = 25.0 #kN
Qf3 = np.array([[ 7*w*L/20 ],
                [ w*L**2/20],
                [ 3*w*L/20 ],
                [-w*L**2/30]])

k3 = Create_k(E, I, L)
S, Pf = Assemble_S_Pf(S, Pf, k3, Qf3, MemNum, NEdof, NSdof, CodeNum)
print 'Qf3 ='
print print_Q(Qf3)
print 'k3 ='
print k3
print '\n'
print 'Pf='
print Pf
print
print 'S ='
print S
print

######## Solve ########
print "Finding Displacements:"
d = solve(S, (P-Pf))   
print 'd ='
print print_d(d)
print

######## POST-PROCESS ########

# Use compatibility to create the {u} vectors for each element
u1 = np.array([ [0],
               d[0],
                [0],
               d[1]])   # Indexing starts with zero (0) in Python

u2 = np.array([ [0],
               d[1],
                [0],
               d[2]])

u3 = np.array([ [0],
               d[2],
                [0],
               d[3]])

# Calculate the member end forces
Q1 = Qf1 + k1.dot(u1)             
Q2 = Qf2 + k2.dot(u2)
Q3 = Qf3 + k3.dot(u3)

print "Finding Q vectors:"
print 'Q1 ='
print print_Q(Q1)
print 'Q2 ='
print print_Q(Q2)
print 'Q3 ='
print print_Q(Q3)

# Calculate the support reactions
R5 = Q1[0]
R6 = Q1[2] + Q2[0]
R7 = Q2[2] + Q3[0]
R8 = Q3[2]

print
print "Reactions: "
print 'R5 =', R5, "kN"
print 'R6 =', R6, "kN"
print 'R7 =', R7, "kN"
print 'R8 =', R8, " kN"