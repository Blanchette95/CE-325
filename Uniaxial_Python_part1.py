'''
Uniaxial Python Code with Member Loads
Units: kN and mm
'''
######## IMPORT MODULES ########
 
import numpy as np # import matrix operators

solve = np.linalg.solve # set up a shorthand for factorization

######## DEFINE FUNCTIONS ########
###
def Create_k(E, A, L):
    """This function creates the k matrix
    Args:
        E (int): Modulus of Elasticity (in GPa)
        A (int): Area of Member (in mm^2)
        L (int): Length of Member (in mm)
    Returns:
        k matrix (array)
    """
    k = (E * A/L) * np.array([[+1, -1],
                              [-1, +1] ])
    return k


def Assemble_S(S, k, MemNum, NEdof, NSdof, CodeNum):
    """Assembles the S matrix
    Args:
        S (array): array of 0's
        k (array): the k matrix from Create_k function
        MemNum (int): Member Number
        NEdof (int): Element degrees of freedom
        NSdof (int): Structural degrees of freedom
        CodeNum (array): Code Numbers
    Returns:
        S (array): The S matrix
    """
    for KRow in range(0, NEdof, 1):
        SRow = CodeNum[MemNum - 1, KRow]
        if SRow <= NSdof:
            for KCol in range(0, NEdof, 1):
                SCol = CodeNum[MemNum - 1, KCol]
                if SCol <= NSdof:
                    S[SRow-1, SCol - 1] = S[SRow - 1, SCol - 1] + k[KRow, KCol]
    return S
###

######## DEFINE SYSTEM FOR ANALYSIS ########

# Define the number of element dofs (always 2 for the uniaxial element)
NEdof = 2

# Define the number of Structural dofs
NSdof = 2

# Define the code numbers
CodeNum = np.array([[3, 1],
                    [1, 2],
                    [2, 4]]) # Here we write out each row explicitly

# Define {P} vector
P = np.array([[300],
              [-800]])
# Initialize S
S = np.zeros((NSdof, NSdof), dtype=float)

# Define geometry and loading for each member. Create the member stiffness
# and then assemble [S]. Note that no units are shown, you just need to be
# be consistent.

# Modulus of Elasticity
E = 200

MemNum = 1
A = 7000
L = 4000
k1 = Create_k(E, A, L)
S = Assemble_S(S, k1, MemNum, NEdof, NSdof, CodeNum)

MemNum = 2
A = 3000
L = 3000
k2 = Create_k(E, A, L)
S = Assemble_S(S, k2, MemNum, NEdof, NSdof, CodeNum)


MemNum = 3
A = 9000
L = 3500
k3 = Create_k(E, A, L)
S = Assemble_S(S, k3, MemNum, NEdof, NSdof, CodeNum)

print 'S =', S
print

######## Solve ########

d = solve(S, P)   #d = linalg.inv(S).dot(P)  alternative solve via invert [S]

print 'd =', d.round(3), '\n'

######## POST-PROCESS ########

# Use compatibility to create the {u} vectors for each element
u1 = np.array([ [0], d[0]])   # Indexing starts with zero (0) in Python
u2 = np.array([d[0], d[1]])
u3 = np.array([d[1],  [0]])

# Calculate the member end forces
Q1 = k1.dot(u1)             
Q2 = k2.dot(u2)
Q3 = k3.dot(u3)

print "Calculating the Member End Forces:"
print 'Q1 =\n', Q1.round(2)
print 'Q2 =\n', Q2.round(2)
print 'Q3 =\n', Q3.round(2), '\n'

# Calculate the member stress
S1 = Q1[1, 0]/7000             
S2 = Q2[1, 0]/3000             
S3 = Q3[1, 0]/9000             

print 'S1 =', S1.round(5)
print 'S2 =', S2.round(5)
print 'S3 =', S3.round(5), '\n'

# Calculate the support reactions
print "Calculating the Reactions:"
R3 = Q1[0, 0]
R4 = Q3[1, 0]

print 'R3 =', R3.round(2)
print 'R4 =', R4.round(2)

