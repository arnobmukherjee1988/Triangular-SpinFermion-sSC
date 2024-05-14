from utils import *
from config import *
import GenerateSpinConfig as spin
import GenerateHamiltonian as ham
import Observables as obs
import Plotting as plot

def main():
    # Print simulation parameters
    print_parameters(print_info = True)
    
    # J11, J22, J33 = ham.neighbor_array_triangular_lattice()
    # for iy in range(Ly):
    #   for ix in range(Lx):
    #     i = (iy * Lx) + ix
    #     print(ix, iy, i, J11[i], J22[i], J33[i])
    
    # Getting Skyrmion spin configuration and calculating skyrmion number
    Sx, Sy, Sz, beta_opt = spin.get_spin()
    Skyrmion_Number = spin.CalculateSkyrmionNumber(Sx, Sy, Sz)
    print("Skyrmion_Number = ", Skyrmion_Number)
    
    # saving the spin configuration in datafile
    datafile = SpinConfig_Type+'_SpinConfig.txt'    
    save_matrices_to_file (datafile, Sx, Sy, Sz)

    # Plot the spin configuration
    plot.SpinPlot(Mat1=Sx, Mat2=Sy, Mat3=Sz, arrow_length=1.5, SaveAs=SpinConfig_Type+'.pdf')

    # Creating the Hamiltonian matrix
    Ham_Matrix = ham.BHZ_hamiltonian(Sx, Sy, Sz)
    
    # diagonalization of the Ham_Matrix to get the eigenvalue array W(n) and eigenfunction matrix Z(n,n)
    W = sp.linalg.eigvalsh(Ham_Matrix)
    
    # QM = obs.QuadrupoleMoment (Z)
    # print(QM)
    
    # save the eigenvalues W in datafile
    datafile = SpinConfig_Type + '_eigenvals.txt'
    save_matrices_to_file (datafile, W)

    # plot the eigen spectrum (every 20 points)
    # X, Y = np.arange(len(W))[::20], W[::20]
    X, Y = np.arange(len(W))[:], W[:]
    plot.plot_data([[X,Y]], xlim=None, ylim=None, linetype='', pointtype='.', 
                    pointsize=20, color='red', xlabel='Eigen Index', ylabel='Eigen Values',
                    title='Eigen Spectrum', SaveAs=SpinConfig_Type+'_eigenvals.pdf')
    
    Energy, dos = obs.dos_cal (W, filename_dos)
    plot.plot_data([[Energy,dos]], xlim=None, ylim=None, linetype='-', pointtype='', 
                    pointsize=20, color='red', xlabel='Energy', ylabel='DOS',
                    title='Total DOS', SaveAs=SpinConfig_Type+'_dos.pdf')

if __name__ == "__main__":
    main()