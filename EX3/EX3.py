import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.signal import welch


class EX3:
    def __init__(self, dataset: int = 12, sample_rate: float = 50000) -> None:
        '''
        Create an instance of the EX3 class

        Args:
            dataset (int): The dataset to use: 1-12 (default 12)
            sample_rate (float): The sample rate of the data [Hz] (default 50000 Hz)
        '''
        self.dataset = dataset-1
        self.fs = sample_rate # Hz
        self.load_data()
        self.dt = 1 / self.fs # s


    def load_data(self) -> None:
        '''
        Load the data from the .mat file
        '''
        data = sp.loadmat('EX3/Exercise3.mat')
        data = data['Jet']
        no_meas = 500000
        
        self.L = float(data['L'][0][self.dataset][0][0])
        self.D = float(data['D'][0][self.dataset][0][0])
        self.nu = float(data['nu'][0][self.dataset][0][0])
        self.V = float(data['V'][0][self.dataset][0][0])

        
        self.t = np.zeros(no_meas)
        self.u = np.zeros(no_meas)

        for i in range(no_meas):
            self.t[i] = float(data['t'][0][self.dataset][i][0])
            self.u[i] = float(data['u'][0][self.dataset][i][0])


    def mean(self, printbol: bool = False) -> None:
        '''
        Calculate the mean velocity ubar

        Args:
            printbol (bool): Print the mean velocity (default False)
        '''
        self.ubar = np.mean(self.u)

        if printbol:
            print(f'The mean velocity is {self.ubar:.4g} m/s')


    def standard_deviation(self, printbol: bool = False) -> None:
        '''
        Calculate the standard deviation sigma of the velocity

        Args:
            printbol (bool): Print the standard deviation (default False)
        '''
        self.sigma = np.std(self.u)

        if printbol:
            print(f'The standard deviation of the velocity is {self.sigma:.4g} m/s')


    def variance(self, printbol: bool = False) -> None:
        '''
        Calculate the variance sigma2 of the velocity

        Args:
            printbol (bool): Print the variance (default False)
        '''
        self.sigma2 = np.var(self.u)

        if printbol:
            print(f'The variance of the velocity is {self.sigma2:.4g} m^2/s^2')
    

    def skewness(self, printbol: bool = False) -> None:
        '''
        Calculate the skewness Su of the velocity

        Args:
            printbol (bool): Print the skewness (default False)
        '''
        self.Su = np.mean((self.u - self.ubar)**3) / self.sigma2**(3/2) # (4.11)

        if printbol:
            print(f'The skewness of the velocity is {self.Su:.4g}')
    

    def kurtosis(self, printbol: bool = False) -> None:
        '''
        Calculate the kurtosis Fu of the velocity

        Args:
            printbol (bool): Print the kurtosis (default False)
        '''
        self.Fu = np.mean((self.u - self.ubar)**4) / self.sigma2**2 # (4.12)

        if printbol:
            print(f'The kurtosis of the velocity is {self.Fu:.4g}')

    
    def turbulence_intensity(self, printbol: bool = False) -> None:
        '''
        Calculate the turbulence intensity I

        Args:
            printbol (bool): Print the turbulence intensity (default False)
        '''
        self.I = np.sqrt(self.sigma2) / self.ubar # (4.53)

        if printbol:
            print(f'The turbulence intensity is {self.I:.4g}')
        

    def fluc_vel(self, plot: bool = False, save: bool = False) -> None:
        '''
        Calculate the fluctuating velocity u'

        Args:
            plot (bool): Plot the fluctuating velocity (default False)
            save (bool): Save the plot in .txt file (default False)
        '''
        self.uprime = self.u - self.ubar
        
        if plot:
            plt.plot(self.t, self.uprime)
            plt.xlabel('Time [s]')
            plt.ylabel('Fluctuating velocity [m/s]')
            plt.title('Fluctuating velocity')
            plt.grid()

        if save:
            np.savetxt('EX3/data/fluctuating_velocity.txt', np.vstack((self.uprime, self.t)).T)


    def probability_density(self, plot: bool = False) -> None:
        '''
        Calculate the probability density function of the velocity

        Args:
            plot: bool - Plot the PDF (default False)
        '''
        x = np.linspace(min(self.u), max(self.u), 1000)
        y = norm.pdf(x, self.ubar, self.sigma)

        gauss = 1 / (self.sigma * np.sqrt(2*np.pi)) * np.exp(- (self.u - self.ubar)**2 / (2 * self.sigma**2)) # (4.3)

        if plot:
            plt.scatter(self.u, gauss, label='Gaussian')
            plt.plot(x, y, label='pdf')
            plt.hist(self.u, bins=100, density=True, alpha=0.5, label='Histogram')
            plt.xlabel('Velocity [m/s]')
            plt.ylabel('Probability density')
            plt.title('Probability density function')
            plt.legend()
            plt.grid()


    def time_correlation(self, plot: bool = False, printbool: bool = False) -> None:
        '''
        Calculate the time correlation function and the Eulerian macro and micro scales

        Args:
            plot (bool): Plot the time correlation function (default False)
            printbool (bool): Print the Eulerian macro and micro scales (default False)
        '''
        max_tau = 1000
        R_E = np.zeros(max_tau)
        R_E[0] = np.mean(self.uprime**2) / self.sigma2 # p. 183
        for tau in range(1, max_tau):
            R_E[tau] = np.mean(self.uprime[:len(self.u) - tau] * self.uprime[tau:]) / self.sigma2 # (4.38)

            if R_E[tau] < 0:
                tau_cross = tau+1
                break

            if tau == max_tau-1:
                raise ValueError('The time correlation function does not cross zero within the range of tau')
        
        self.R_E = R_E[:tau_cross] # First negative value is included

        self.T_E = np.trapezoid(self.R_E, self.t[:tau_cross]) # Eulerian macro scale (4.45)

        self.tau_E = np.sqrt(2 * self.sigma2 / np.mean((np.diff(self.uprime)/self.dt)**2))

        if printbool:
            print(f'Zero crossing at tau =  {tau_cross} with R_E = {R_E[tau_cross]:.4g}')
            print(f'The Eulerian macro scale is {self.T_E:.4g} s')
            print(f'The Eulerian micro scale is {self.tau_E:.4g} s')


        if plot:
            np.savetxt('EX3/data/time_correlation.txt', np.vstack((self.t[:tau_cross], self.R_E)).T)
            plt.plot(self.t[:tau_cross], self.R_E)
            plt.xlabel('Time (tau) [s]')
            plt.ylabel('Time correlation function (R_E(tau))')
            plt.title('Time correlation function')
            plt.grid()


    def micro_macro(self, printbool:bool=False) -> None:
        '''
        Calculates the approximate micro and macro turbulent scales with Taylor's frozen turbulence hypothesis.

        Args:
            printbool (bool): Print the micro and macro turbulent scales (default False)
        
        Attributes:
            lambda_f (float): The micro turbulent scale
            Lambda_f (float): The macro turbulent scale
        '''
        self.lambda_f = self.tau_E * self.ubar # (4.55)
        
        self.Lambda_f = self.T_E * self.ubar # (4.55)

        if printbool:
            print(f'The micro turbulent scale is {self.lambda_f:.4g} m')
            print(f'The macro turbulent scale is {self.Lambda_f:.4g} m')


    def dissipation_rate(self, printbool: bool = False) -> None:
        '''
        Calculates the dissipation rate of the turbulent kinetic energy

        Args:
            printbool (bool): Print the dissipation rate (default False)
        
        Attributes:
            epsilon (float): The dissipation rate
        '''
        self.epsilon = 30 * self.nu * self.sigma2 / self.lambda_f**2 # (4.83)

        if printbool:
            print(f'The dissipation rate is {self.epsilon:.4g} m^2/s^3')


    def kolomogorov_length(self, printbool: bool = False) -> None:
        '''
        Calculates the Kolmogorov length scale

        Args:
            printbool (bool): Print the Kolmogorov length scale (default False)

        Attributes:
            eta_K (float): The Kolmogorov length scale
        '''
        self.eta_K = (self.nu**3 / self.epsilon)**(1/4) # (4.86)

        if printbool:
            print(f'The Kolmogorov length scale is {self.eta_K:.4g} m')


    def kolomogorov_time(self, printbool: bool = False) -> None:
        '''
        Calculates the Kolmogorov time scale

        Args:
            printbool (bool): Print the Kolmogorov time scale (default False)

        Attributes:
            tau_K (float): The Kolmogorov time scale
        '''
        self.tau_K = (self.nu / self.epsilon)**(1/2) # (4.87)

        if printbool:
            print(f'The Kolmogorov time scale is {self.tau_K:.4g} s')


    def energy_spectra(self, plot: bool = False, printbool: bool = False) -> None:
        '''
        Computes the energy spectra of the fluctuating velocity and checks if the integral matches the variance of the velocity

        Args:
            plot (bool): Plot the energy spectra (default False)
            printbool (bool): Print the integral of the energy spectra (default False)
        
        Raises:
            ValueError: If the integral of the energy spectra does not match the variance of the velocity
        '''
        # Er nappet direkte fra side 216
        N = len(self.uprime)
        dt = 1/self.fs
        f_N = 1 / (2*dt)
        df = f_N / (N / 2)
        self.f = np.arange(0, f_N, df)
        U = np.fft.fft(self.uprime)/N
        A = np.abs(2*U[0:N//2])
        self.S = 0.5*A**2/df

        # Smoothing
        nw = 100
        windowlenght = round(N / nw)
        freqs, psd = welch(self.uprime, fs=self.fs, nperseg=windowlenght) 

        area = np.trapezoid(self.S, self.f)
        if not np.isclose(area, self.sigma2, rtol=1e-6):
            raise ValueError(f'The integral of the energy spectra does not match the variance of the velocity. The integral is {area} m^2/s^2 and the variance is {self.sigma2} m^2/s^2')

        
        if printbool:
            print(f'The integral of the energy spectra is {area:.4g} m^2/s')

        if plot:
            np.savetxt('EX3/data/energy_spectra.txt', np.vstack((self.f, self.S)).T)
            np.savetxt('EX3/data/energy_spectra_welch.txt', np.vstack((freqs, psd)).T)

            plt.loglog(self.f, self.S, label='Energy spectra')
            plt.loglog(freqs, psd, label='Smoothed energy spectra')
            plt.legend()
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Energy spectra m^2/s')
            plt.title('Energy spectra of the fluctuating velocity')
            plt.grid()


    def wave_spectrum(self, plot:bool = False, printbool:bool=False) -> None:
        '''
        Computes the wave spectrum of the fluctuating velocity

        Args:
            printbool (bool): Print the wave spectrum integral (default False)
            plot (bool): Plot the wave spectrum and compared to the von Karmen model spectrum (default False)
            
        Attributes:
            F (np.array): The wave spectrum
            k (np.array): The wave number
        '''
        self.F = self.ubar/(4*np.pi) * self.S  # (4.125)
        self.k = 2*np.pi*self.f/self.ubar # P. 707

        area = np.trapezoid(self.F, self.k)
        if not np.isclose(area, self.sigma2/2, rtol=1e-6):
            raise ValueError(f'The integral of the wave spectrum does not match hald of the variance of the velocity. The integral is {area} m^2/s^2 and the variance is {self.sigma2} m^2/s^2')

        VK = self.Lambda_f * self.sigma2 / np.pi * 1/(1 + 70.78 * (self.k * self.Lambda_f/(2*np.pi))**2)**(5/6) # (4.131)

        if printbool:
            print(f'The integral of the wave spectrum is {area:.4g} m^2/s^2')

        if plot:
            np.savetxt('EX3/data/wave_spectrum.txt', np.vstack((self.k, self.F)).T)
            np.savetxt('EX3/data/von_karman_spectrum.txt', np.vstack((self.k, VK)).T)

            plt.loglog(self.k, self.F, label='Wave spectrum')
            plt.loglog(self.k, VK, label='von Karmen model')
            plt.legend()
            plt.xlabel('Wave number [1/m]')
            plt.ylabel('Wave spectrum m^2/s')
            plt.title('Wave spectrum of the fluctuating velocity')
            plt.grid()


if __name__ == '__main__':
    if False: # Test
        a = np.array([1, 2, 3])
        b = np.array([4, 5, 6])
        c = np.vstack((a, b)).T
        print(c)


    if False: # Part AI1
        print('Part AI1')
        ex1 = EX3()
        ex1.mean(printbol=True)
        ex1.fluc_vel(plot=True, save=True)
        plt.show()
    

    if False: # Part AI2
        print('Part AI2')
        ex1 = EX3()
        ex2.mean(printbol=True)
        ex2.variance(printbol=True)
        ex2.skewness(printbol=True)
        ex2.kurtosis(printbol=True)
        ex2.turbulence_intensity(printbol=True)


    if False: # Part AI3
        print('Part AI3')
        ex3 = EX3()
        ex3.mean()
        ex3.standard_deviation()
        ex3.probability_density(plot=True)
        plt.show()


    if False: # Part AI4
        print('Part AI4')
        ex4 = EX3()
        ex4.mean()
        ex4.variance(printbol=False)
        ex4.fluc_vel()
        ex4.time_correlation(plot=True, printbool=True)
        plt.show()


    if False: # Part AI5
        print('Part AI5')
        ex5 = EX3()
        ex5.mean()
        ex5.variance()
        ex5.fluc_vel()
        ex5.time_correlation()
        ex5.micro_macro(printbool=True)
        
    
    if False: # Part AI6
        print('Part AI6')
        ex6 = EX3()
        ex6.mean()
        ex6.variance()
        ex6.fluc_vel()
        ex6.time_correlation()
        ex6.micro_macro()
        ex6.dissipation_rate(printbool=True)
        ex6.kolomogorov_length(printbool=True)
        ex6.kolomogorov_time(printbool=True)


    if False: # Part AI7
        print('Part AI7')
        ex3 = EX3()
        ex3.mean()
        ex3.variance(printbol=True)
        ex3.fluc_vel()
        ex3.energy_spectra(plot=True, printbool=True)
        plt.show()


    if True: # Part AI8
        print('Part AI8')
        ex8 = EX3()
        ex8.mean()
        ex8.variance()
        ex8.fluc_vel()
        ex8.time_correlation()
        ex8.micro_macro()
        ex8.energy_spectra()
        ex8.wave_spectrum(plot=True)
        plt.show()
        

