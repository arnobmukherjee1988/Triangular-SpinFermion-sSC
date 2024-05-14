from utils import *

# Lattice Type and lattice parameters
Lx = 30 ; Ly = 30
Nsites = Lx * Ly
BoundaryConditionList = ['PBCx + PBCy',
                         'OBCx + PBCy',
                         'PBCx + OBCy',
                         'OBCx + OBCy']
BoundaryCondition = BoundaryConditionList[3]
BraviasLattice = "TriangularLattice"

# Hamiltonian coupling parameters
# hopping parameters
tx = 1.0 ; ty = 1.0 ; td = 1.0
# superconducting order parameters
Delta_s = 0.0
Delta_px = 0.0 ; Delta_py = 0.0
Delta_dx = 0.0 ; Delta_dy = 0.0
# Hunds coupling parameters
J_Hundx = 0.0 ; J_Hundy = 0.0 ; J_Hundz = 0.0
# chemical potential & crystal field splitting
mu = 0.0 ; epsilon0 = 0.0
# the SOC strength
Lambda_SOC_x = 0.0 ; Lambda_SOC_y = 0.0 ; Lambda_SOC_d = 0.0
Lambda_WD_x = 0.0 ; Lambda_WD_y = 0.0 ; Lambda_WD_d = 0.0

# Hamiltonian size parameters
# degrees of freedom, spin (SDOF), particle-hole (PHDOF), orbital (ODOF)
SDOF = 1 ; ODOF = 1 ; PHDOF = 1 ; TDOF = PHDOF * ODOF * SDOF
PHFactor = 1.0 if PHDOF == 2 else 1.0
n = Nsites * TDOF

# Spin Configuration Parameters
SpinConfigTypeList = ["NeelSkyrmion",      # 0
                      "BlochSkyrmion",     # 1
                      "AntiSkyrmion",      # 2
                      "SpinSpiral",        # 3
                      "Ferro",             # 4
                      "AntiFerro",         # 5
                      "CustomSkyrmion",    # 6
                      "PRB_109_L041409",   # 7
                      "PRR_4_013225_SWCB4"] # 8
SpinConfig_Type = SpinConfigTypeList[7]
SkyrDiameter = 8
vorticity = 1.0
helicity = 0.0
polarity = -1.0

qx = PI/2.0 ; qy = PI/2.0 # spiral wave vector components

# Skyrmion optimization paramters
no_Beta = 1000
Beta_min = 0.05
Beta_max = 0.3

# generic offset
offset = 0.00001

BC_Error = 'Invalid BoundaryCondition! You have entered ' + BoundaryCondition + '\n' + 'Valid options are: PBCx + PBCy, OBCx + PBCy, PBCx + OBCy, OBCx + OBCy'
line1 = 'Invalid SpinConfig_Type! You have entered ' + SpinConfig_Type + '\n'
line2 = 'Valid options are: "NeelSkyrmion", "BlochSkyrmion", "AntiSkyrmion", "SpinSpiral", "Ferro", "AntiFerro"'
SpinConfig_Error = line1 + line2

def print_parameters(print_info=True):
  if not print_info:
    return
  print("-" * 78)
  
  print("\nParameters:")
  print(f"  Bravais Lattice: {BraviasLattice}")
  print(f"  Lattice Size: Lx = {Lx}, Ly = {Ly}")
  # print(f"  Number of Sites: {Nsites}")
  # print(f"  Degrees of Freedom (SDOF, PHDOF): {SDOF}, {PHDOF}")
  # print(f"  Total Degrees of Freedom (n): {n}")
  # print(f"  Shift: {Shift}")
  print(f"  Boundary Condition: {BoundaryCondition}")
  print(f"  Tight-Binding Parameters: tx = {tx}, ty = {ty}")
  # print(f"  Hund's Coupling: J_Hund = {J_Hund}")
  print(f"  Chemical Potential: mu = {mu}")
  print(f"  S-wave pairing: Delta_s = {Delta_s}")
  print(f"  P-wave pairing: Delta_px = {Delta_px}, Delta_py = {Delta_py}")
  print(f"  D-Wave pairing: Delta_dx = {Delta_dx}, Delta_dy = {Delta_dy}")
  print(f"  Spin Configuration Type: {SpinConfig_Type}")
  if SpinConfig_Type == 'SpinSpiral':
    print(f"  Wave vector: qx = {qx}, qy = {qy}")
  elif SpinConfig_Type == 'NeelSkyrmion' or SpinConfig_Type == 'BlochSkyrmion' or SpinConfig_Type == 'AntiSkyrmion':
    print(f"  Skyrmion Parameters: SkyrDiameter = {SkyrDiameter}, SkyrRadius = {int(SkyrDiameter/2.0)}")
    print(f"  Number of Skyrmions Along X, Y: {int(Lx/SkyrDiameter)}, {int(Ly/SkyrDiameter)}")
    print(f"  Total Skyrmions = {int(Lx/SkyrDiameter) * int(Ly/SkyrDiameter)}")
  
  print("-" * 78)



'''
Put errors for wrong boundary condition and spin config type
Fix the if conditions for 'Create_SpinConfig' in the GenerateSpinConfig file 

'''

energy_of_interest = 0.0
site_of_interest = 1
eta = 0.01#1e-2
DOS_NE = 200 ; LDOS_NE = 200

filename_spin = 'spin_config.txt'
filename_eigenval = 'eigenvals.txt'
filename_dos = 'dos.txt'
filename_ldos_E_fix = 'ldos_EnergyFixed.txt'
filename_ldos_Site_fix = 'ldos_SiteFixed.txt'