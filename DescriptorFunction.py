import numpy as np
from Input_Parser import QE_Input

class Sym_Function:
    '''
    Sym_Function:
        Symmetry function based on modified ACSF
    '''
    def __init__(self, input: QE_Input): 
        '''
        Parse input file into required data
        Input: QE_Input from QE_Class
        '''
        self.alat = input.get_alat()/0.529 #in Bohr
        self.input = input
        self.GenerateImage()
        self.coordinates = []
        self.atom_species = []
        for atom in self.input.dict['ATOMIC_POSITIONS']:
            self.coordinates.append(np.array(atom[1:],dtype=float)*self.alat)
            self.atom_species.append(atom[0])
        self.cell = self.input.dict['CELL_PARAMETERS']
    
    def cutoff_Func(self,Rc:float,Rij: float):
        '''
        Cutoff function as defined in ACSF paper
        '''
        if Rij > Rc:
            # When Rij exceed cutoff. set to 0
            return 0
        return 0.5*(np.cos(np.pi*(Rij/Rc))+1)

    def Compute(self,atom, Rc: float, eta = 0, Rs=0, kappa=0):
        '''
        The main computation function, this is exactly 
        the same as the one in ACSF paper
        '''
        sum = 0
        atom = atom -1

        coord = np.copy(self.coordinates)
        Ri=self.coordinates[atom]
        Rj=np.delete(coord,atom,0)
        Rij = np.sqrt(np.sum((Ri-Rj)*(Ri-Rj),axis=1))
        Rij = Rij[Rij < Rc]
        sum =  np.sum(np.cos(kappa*Rij)*np.exp(-eta*(Rij-Rs)**2)*0.5*(np.cos(np.pi*(Rij/Rc))+1))
        return sum
    
    def Compute_Trunc(self,atom,Rb:float, Rc: float, eta = 0, Rs=0, kappa=0):
        '''
        This is the modified version of 
        compute. We also exclude atom within 
        radius Rb. This is to get descriptor
        B and C.
        '''
        sum = 0
        atom = atom -1
        Rs = Rb
        coord = np.copy(self.coordinates)
        Ri=self.coordinates[atom]
        Rj=np.delete(coord,atom,0)
        Rij = np.sqrt(np.sum((Ri-Rj)*(Ri-Rj),axis=1))
        Rij = Rij[Rij > Rb]
        Rij = Rij[Rij < Rc]
        sum =  np.sum(np.cos(kappa*Rij)*np.exp(-eta*(Rij-Rs)**2)*0.5*(np.cos(np.pi*(Rij/Rc))+1))
        return sum
        
    def G1(self,atoms,Rc):
        '''
        As defined in ACSF paper.
        '''
        if type(atoms) is not list:
            atoms = [atoms]
        sum = 0 
        for atom in atoms:
            sum += self.Compute(atom,Rc=Rc)
        return sum
    
    def G1_mod(self,atoms,Rb,Rc):
        '''
        Modified G1, which do not consider
        atom within range Rb. Uses 
        Compute_trunc
        '''
        if type(atoms) is not list:
            atoms = [atoms]
        sum = 0 
        for atom in atoms:
            sum += self.Compute_Trunc(atom,Rb=Rb,Rc=Rc)
        return sum
    
    def G1_neigh(self,atom,Rb,Rc):
        '''
        This calculates G1 for each 
        neigbor found between Rb and Rc.
        '''
        sum = 0
        atom = atom -1

        coord = np.copy(self.coordinates)
        Ri=self.coordinates[atom]
        Rj=coord
        Rij = np.sqrt(np.sum((Ri-Rj)*(Ri-Rj),axis=1))
        for i in np.where((Rij < Rc) * Rij > Rb )[0]:
            sum += self.G1(i+1,2.8/0.529)
        sum /= len(np.where((Rij < Rc) * Rij > Rb )[0])
        print( len(np.where((Rij < Rc) * Rij > Rb )[0]))
        return sum

   
    def GenerateImage(self):
        '''
        To make cell larger so that the neighbouring 
        image atoms are captured. Ensure all local
        details are included.
        '''
        self.input.supercell(2,2)
        
        
