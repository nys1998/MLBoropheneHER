import numpy as np
import re

class QE_Output: 
    ''' 
    Class for QE Input File
    Input:
        path : type str : path to desired input file
    '''
    def __init__(self,path:str):
        self.path = path
        self.finalenergy = 0
        self.converged = False
        self.finalforce = 0
        self.coordinate = []
        self.cell = []
        self.alat = 0
        self.load_file()
        self.celltype = None
        self.atomtype = None
    
    def load_file(self):
        with open(self.path,'r') as f:
            data = f.read()
            line = []
            self.alat = float(data[data.find('lattice parameter (alat)' ):data.find('a.u.')].split()[-1])
            if data.find('final') != -1:
                self.converged = True
            iteration = ([m.start() for m in re.finditer('!', data)])
            if len(iteration) == 1:
                data = data[iteration[0]:].split('\n')
            else:
                data = data[iteration[len(iteration)-2]:iteration[len(iteration)-1]].split('\n')
            self.finalenergy = float(re.findall(r'-[0-9]+\.[0-9]+',data[0])[0])

            read_atom = False
            for i in range(len(data)):
                line = data[i]
                if 'Total force' in line:
                    self.finalforce=float(line.split()[3])
                if 'CELL_PARAMETERS' in line:
                    if 'alat' in line: 
                        self.celltype = 'alat'
                    if 'angstrom' in line:
                        self.celltype = 'angstrom'
                    self.cell.append([float(m) for m in data[i+1].strip().split()])
                    self.cell.append([float(m) for m in data[i+2].strip().split()])
                    self.cell.append([float(m) for m in data[i+3].strip().split()])
                if 'ATOMIC_POSITIONS' in line:
                    if 'alat' in line: 
                        self.atomtype = 'alat'
                    if 'angstrom' in line:
                        self.atomtype = 'angstrom'  
                    read_atom = True
                    continue
                if read_atom:
                    if len(line.strip().split()) !=4:
                        break
                    self.coordinate.append(line.split())

class QE_Input:
    ''' 
    Class for QE Input File
    Input:
        path : type str : path to desired input file
    '''
    def __init__(self, path:str):
        self.path = path
        self.dict = {}
        self.EOF = 0
        self.atom_isalat = True
        self.prefix = ''
        self.ntyp = 0
        self.nat = 0
        self.ibrav = 0
        self.typorder = []
        self.load_file()
        pass

    def load_file(self):
        with open(self.path) as f:
            line_vec= []
            category = []
            data = f.read().split('\n')
            temp = []
            for l in range(len(data)): 
                if data[l].strip() == '':
                    continue
                temp.append(data[l])
            data = temp
            self.EOF = len(data)
            for count, line in enumerate(data):
                line = line.strip()
                if line[0] == '&':
                    line_vec.append(count)
                    category.append(line.upper())
                    temp = []
                if line.upper() == 'ATOMIC_SPECIES':
                    line_vec.append(count)
                    category.append(line.upper())
                if line.upper()[:15] == 'CELL_PARAMETERS':
                    line_vec.append(count)
                    category.append(line.upper()[:15])                    
                if line.upper()[:16] == 'ATOMIC_POSITIONS':
                    line_vec.append(count)
                    category.append(line.upper()[:16])
                    if line.upper().find('ANGSTROM') != -1:
                        self.atom_isalat = False 
                if line.upper()[:8] == 'K_POINTS':
                    line_vec.append(count)
                    category.append(line.upper())
        line_vec.append(self.EOF)
        for c,i in enumerate(line_vec):
            temp = []
            temp_dict = {}
            if c > len(category)-1: 
                break
            if category[c][0] == '&':
                for j in range(i+1,len(data)):
                    line = data[j].strip()
                    if j in line_vec or line == '/':
                        self.dict.update({category[c]:temp_dict})
                        break
                    key = line.split('=')[0].strip().lower()
                    value = line.split('=')[1].strip()
                    if key == 'prefix': self.prefix = value.strip(',')
                    if key == 'ntyp': self.ntyp = int(value.strip(','))
                    if key =='nat':self.nat = int(value.strip(','))
                    if key == 'ibrav':self.ibrav = value.strip(',')
                    temp_dict.update({key:value})
            else:
                for j in range(i+1,len(data)):
                    line = data[j].strip()
                    if j in line_vec or line== '/':
                        self.dict.update({category[c]:temp})
                        break
                    if j+1 == len(data):
                        temp.append(line)
                        self.dict.update({category[c]:temp})
                        break
                    try:
                        temp.append([float(k) for k in line.split()])
                    except:
                        temp.append([k for k in line.split()])

    def update(self): 
        self.dict['&CONTROL']['prefix'] = f"{self.prefix}"
        self.dict['&SYSTEM']['nat'] = self.nat
        self.dict['&SYSTEM']['ntyp'] = self.ntyp

    def print(self,f=''):
        # Print input file
        self.update()
        if f == '': #If file name is blank, print to console
            for k in self.dict:
                print(k)
                for i in self.dict[k]:
                    if '=' in i:
                        print(f'\t{i}')
                    else:
                        print(f'{i}')
        else:
            if f.find('.in') > -1:
                f = f.strip('.in')
            with open(f'{f}.in','w') as g:
                for k in self.dict:
                    if k == 'ATOMIC_POSITIONS':
                        if self.atom_isalat == True:
                            kp = k + ' {alat}'
                            g.write(f'{kp}\n')
                        else:
                            kp = k + ' {angstrom}'
                            g.write(f'{kp}\n')
                    elif k == 'CELL_PARAMETERS':
                        kp = k + ' (alat)'
                        g.write(f'{kp}\n')
                    else:
                        g.write(f'{k}\n')
                    for i in self.dict[k]:
                        if k[0] == '&':
                            # if i.strip() == 'nat': 
                            #     self.dict[k][i] = self.nat
                            # if i.strip() == 'ntyp': 
                            #     self.dict[k][i] = self.ntyp
                            g.write(f'\t{i}={self.dict[k][i]}\n')
                            continue
                        if k == 'ATOMIC_POSITIONS' or k == 'ATOMIC_SPECIES' or k == 'CELL_PARAMETERS':
                            line = ''
                            for j in i:
                                try:
                                    j=float(j)
                                    if k =='ATOMIC_POSITIONS'or k == 'CELL_PARAMETERS':
                                        line += f'{j:.8f}\t'
                                    else:
                                        line+= f'{j}\t'
                                except:
                                    line+= f'{j}\t'
                            line += '\n'
                            g.write(line)
                        else:
                            g.write(f'{i}\n')
                    if k.find('&') == 0:
                        g.write('/ \n')

    def get_crystal_info(self):
        sys = self.dict['&SYSTEM']
        crystal={'ibrav':0,'BohrAng':'Bohr','a1':0,'a2':0,'a3':0}
        for i in sys:
            if i == '':
                continue
            i = i.split('=')
            if i[0].strip().lower().find('celldm(1)') != -1:
                return np.round(float(i[1].strip().strip(','))*0.529,decimals=4)
            if i[0].strip().find('A') != -1:
                return float(i[1].strip())    
               
    def get_alat(self):
        # Always return in angstrom
        sys = self.dict['&SYSTEM']
        for i in sys:
            if i == '':
                continue
            i = i.split('=')
            if i[0].strip().lower() == 'celldm(1)':
                return np.round(float(sys[i[0]].strip().strip(','))*0.529,decimals=4)
            if i[0].strip() == 'a':
                return float(sys[i[0]].strip(','))
            
    def supercell(self,m,n,z=1):
        # First make sure the coordinates are in alat form, angstrom wouldn't work
        self.angstrom_to_alat()

        # Stash away the original atoms
        old_coordinate = self.dict['ATOMIC_POSITIONS']
        cell = np.array(self.dict['CELL_PARAMETERS'])

        # Add new atoms
        temp =[]
        for i in range(abs(m)):
            for j in range(abs(n)):
                if i+j==0: continue
                add_vector = np.add(np.sign(m)*i*cell[0],np.sign(n)*j*cell[1])
                for atom in old_coordinate:
                    num_coord = np.array(atom[1:],dtype=float)+ add_vector     
                    temp.append(np.hstack(([atom[0]],num_coord)))
        
        # Expands the cell and store
        cell = np.multiply(cell,[m,n,z])
        self.dict['ATOMIC_POSITIONS'] = np.vstack((old_coordinate,temp))
        self.dict['CELL_PARAMETERS'] = cell
        self.nat = len(self.dict['ATOMIC_POSITIONS'])
    
    def fold_to_1x1(self,m,n):
        m = int(m)
        n = int(n)

        # First make sure the coordinates are in alat form, angstrom wouldn't work
        self.angstrom_to_alat() 
        # Stash away the original atoms
        old_coordinate = self.dict['ATOMIC_POSITIONS']
        cell = np.array(self.dict['CELL_PARAMETERS'])
        divisor = np.round(cell[0])
        new_cell = np.copy(cell)
        new_cell[0] /= m
        new_cell[1] /= n

        # New atoms
        temp =[]

        for atom in old_coordinate:
            if float(atom[1]) >= new_cell[0][0] or float(atom[2]) >= new_cell[1][1]:
                continue
            else:
                temp.append(atom)

        # Store
        self.dict['ATOMIC_POSITIONS'] = temp
        self.dict['CELL_PARAMETERS'] = new_cell
        self.nat = len(self.dict['ATOMIC_POSITIONS'])

    def angstrom_to_alat(self):
        if self.atom_isalat == True:
            return None
        divisor = self.get_alat()
        self.atom_isalat = True
        old_coordinate = self.dict['ATOMIC_POSITIONS']
        temp =[]
        for atom in old_coordinate:
            num_coord = np.array(atom[1:],dtype=float)/divisor   
            temp.append(np.hstack(([atom[0]],num_coord)))     
        self.dict['ATOMIC_POSITIONS'] = temp

    def strain(self,m,n):
        m = float(m)/100+1
        n = float(n)/100+1
        self.dict['CELL_PARAMETERS'][0] = np.multiply(self.dict['CELL_PARAMETERS'][0],m)
        self.dict['CELL_PARAMETERS'][1] = np.multiply(self.dict['CELL_PARAMETERS'][1],n)
        for j in range(len(self.dict['ATOMIC_POSITIONS'])):
            self.dict['ATOMIC_POSITIONS'][j][1] = np.multiply(float(self.dict['ATOMIC_POSITIONS'][j][1]), m)
            self.dict['ATOMIC_POSITIONS'][j][2] = np.multiply(float(self.dict['ATOMIC_POSITIONS'][j][2]), n)
        
    def get_output(self,output:QE_Output):
        self.dict['ATOMIC_POSITIONS'] = output.coordinate
        if output.cell != []:
            self.dict['CELL_PARAMETERS'] = output.cell

    def add_adsorbate(self,species: str, atom_no: int, distance = 1, pseudo_name = '', elem_mass = 1):
        new_coord = [species]
        for i in range(3):
            new_coord.append(self.dict['ATOMIC_POSITIONS'][atom_no-1][i+1])
        new_coord[3] = str(float(new_coord[3])+np.round(1/(self.get_alat()),7))
        self.dict['ATOMIC_POSITIONS'].append(new_coord)
        if pseudo_name != '':
            self.dict['ATOMIC_SPECIES'].append([species,str(elem_mass),pseudo_name])
        self.ntyp += 1
        self.nat += 1