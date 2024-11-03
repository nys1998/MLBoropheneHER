from SymmetryFunction import Sym_Function
from QE_Class import QE_Input
import numpy as np
import copy

def Generate_Feature(qe,atom):
    sym = Sym_Function(qe)
    temp = []
    R = [0,2.8,4.2,5.3] # In angstrom unit
    for i in range(3):
        # Here it generates A, B, C
        temp.append(sym.G1_mod(atom,Rb=R[i]/0.529,Rc=R[i+1]/0.529))
    # Here it generates alpha
    temp.append(sym.G1_neigh(atom,Rb=R[0]/0.529,Rc=R[1]/0.529))
    # Here it generates beta
    temp.append(sym.G1_neigh(atom,Rb=R[1]/0.529,Rc=R[2]/0.529))
    return temp

# True to generate testing data only
train = False

'''
qe: QE input file, converted to QE input type
dict: contains atom number of each defined site
CN: coordination number
'''
if train == True:
    training_data = []
    fil = open('training_name.dat', 'w')

    #alpha1
    dict = {'S1': 2, 'S2': 50,'S3': 22}
    CN = {'S1':5,'S2':6,'S3':5}
    for S in dict:
        # Load file containing info on strain and dG
        GH = np.loadtxt(f'./alpha1/{S}.txt')
        ind = np.lexsort((GH[:,1],GH[:,0]))
        GH =GH[ind]
        for line in GH:
            '''
            Load each input file (according to strain)
            Generate sym function and store along with dG
            '''
            temp = []
            file = f"./alpha1/a{int(line[0])}b{int(line[1])}.in"
            qe = QE_Input(file)
            temp = Generate_Feature(qe,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'alpha1 {S} {int(line[0])} {int(line[1])} {CN[S]}\n')

    #Beta1
    f = QE_Input('./beta1/beta1.in')
    dict =   {'S1': 34, 'S2': 22, 'S3': 26, 'S4':38,'S5':18,'S6':2,'S7':42} 
    CN =   {'S1': 5, 'S2': 5, 'S3': 4, 'S4':5,'S5':5,'S6':6,'S7':6} 
    for S in dict:
        GH = np.loadtxt(f'./beta1/{S}.txt')
        for line in GH:
            temp = []
            g =copy.deepcopy(f)
            temp = Generate_Feature(g,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'B1 {S} 0 0 {CN[S]}\n')
            break

    # Beta12
    dict = {'S1':28,'S2':27,'S3':30}
    CN = {'S1':4,'S2':6,'S3':5}
    for S in dict:
        GH = np.loadtxt(f'./b12/{S}train.txt')
        ind = np.lexsort((GH[:,1],GH[:,0]))
        GH =GH[ind]
        for line in GH:
            temp = []
            file = f"./b12/a{int(line[0])}b{int(line[1])}.in"
            qe = QE_Input(file)
            temp = Generate_Feature(qe,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'B12 {S} {int(line[0])} {int(line[1])} {CN[S]}\n')

    # trigonal
    dict = {'S1':30}
    CN={'S1':'7'}
    f = QE_Input('./trigonal/trigonal.in')
    for S in dict:
        GH = np.loadtxt(f'./trigonal/{S}.txt')
        ind = np.lexsort((GH[:,1],GH[:,0]))
        GH =GH[ind]
        for line in GH:
            file = f"./trigonal/a{int(line[0])}b{int(line[1])}.in"
            f = QE_Input(file)
            temp = Generate_Feature(f,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'trigonal {S} {int(line[0])} {int(line[1])} {CN[S]}\n')

    # x3
    dict = {'S1':48 ,'S2':47}
    CN = {'S1':4,'S2':5}
    f = QE_Input('./x3/x3.in')
    for S in dict:
        GH = np.loadtxt(f'./x3/{S}.txt')
        ind = np.lexsort((GH[:,1],GH[:,0]))
        GH =GH[ind]
        for line in GH:
            file = f"./x3/a{int(line[0])}b{int(line[1])}.in"
            f = QE_Input(file)
            temp = Generate_Feature(f,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'x3 {S} {int(line[0])} {int(line[1])} {CN[S]}\n')

    #alpha
    dict = {'S1':31,'S2':32}
    CN = {'S1':5,'S2':6}
    for S in dict:
        GH = np.loadtxt(f'./alpha/{S}.txt')
        ind = np.lexsort((GH[:,1],GH[:,0]))
        GH =GH[ind]
        for line in GH:
            temp = []
            file = f"./alpha/a{int(line[0])}b{int(line[1])}.in"
            f = QE_Input(file)
            temp = Generate_Feature(f,dict[S])
            temp.append(line[2])
            training_data.append(temp)
            fil.write(f'alpha {S} {int(line[0])} {int(line[1])} {CN[S]}\n')
    fil.close()
    with open('training9Sep.dat','w') as f:
        for line in training_data:
            l = ''
            for i in line:
                l += f'{i}\t'
            f.write(l)
            f.write('\n')

# Generate testing file
if train == False:
    temp = []
    for i in [-2,-1,0,1,2]:
        #Beta1
        f = QE_Input(f'./beta1test/a{i}b{i}.scf.in')
        dict =   {'S1': 34, 'S2': 22} 
        CN =   {'S1': 5, 'S2': 5, 'S3': 4, 'S4':5,'S5':5,'S6':6,'S7':6} 
        for S in dict:
            g =copy.deepcopy(f)
            temp.append(Generate_Feature(g,dict[S]))
            # fil.write(f'B1 {S} 0 0 {CN[S]}\n')
    with open('beta1test.dat','w') as f:
        for line in temp:
            l = ''
            for i in line:
                l += f'{i}\t'
            f.write(l)
            f.write('\n')
