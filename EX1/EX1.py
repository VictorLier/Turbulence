import scipy.io as sp
from scipy.integrate import cumulative_trapezoid as cumtrapz
import numpy as np
import matplotlib.pyplot as plt

def save_as_txt(filename: str, x_data, y_data) -> None:
    '''
    Save the data as a .txt file
    '''
    np.savetxt('EX1/Plot/'+filename, np.column_stack((x_data, y_data)))
    


class EX:
    def __init__(self):
        '''
        Starts with loading the data from the .mat file
        '''
        self.load_data()
        self.b = 0.3    # Width of the channel

    def load_data(self):
        '''
        Loading the data from the .mat file
        '''
        # Load the data
        data = sp.loadmat('EX1/Exercise1.mat')
        data = data['Channel']
        data = data[0]

        data_n = 23 # Number of data points

        self.nu = data[0][6][0][0]

        self.h = data[0][5][0][0]

        self.y = np.zeros(data_n+2)
        self.y[0] = 0
        self.y[-1] = self.h
        for i in range(data_n):
            self.y[i+1] = data[i][4][0][0]

        self.v = [0]*(data_n+1)
        self.v[0] = 0
        for i in range(data_n):
            j_in = len(data[i][3])
            v_temp = np.zeros(j_in)
            for j in range(j_in):
                v_temp[j] = data[i][3][j][0]
            self.v[i+1] = v_temp

        self.u = [0]*(data_n+2)
        self.u[0] = 0
        self.u[-1] = 0.3
        for i in range(data_n):
            j_in = len(data[i][2])
            u_temp = np.zeros(j_in)
            for j in range(j_in):
                u_temp[j] = data[i][2][j][0]
            self.u[i+1] = u_temp

        self.tt = [0]*data_n
        for i in range(data_n):
            j_in = len(data[i][1])
            tt_temp = np.zeros(j_in)
            for j in range(j_in):
                tt_temp[j] = data[i][1][j][0]
            self.tt[i] = tt_temp

        self.t = [0]*data_n
        for i in range(data_n):
            j_in = len(data[i][0])
            t_temp = np.zeros(j_in)
            for j in range(j_in):
                t_temp[j] = data[i][0][j][0]
            self.t[i] = t_temp

    def mean_velocity(self):
        '''
        Calculate the mean velocity for ubar, vbar and uvbar
        '''
        self.ubar = np.zeros(23+2)
        self.ubar[0] = 0
        self.ubar[-1] = 0.3
        for i in range(23):
            self.ubar[i+1] = (np.sum(self.u[i+1] * self.tt[i])) / np.sum(self.tt[i])    # eq. 10.1
        
        self.vbar = np.zeros(23)
        for i in range(23):
            self.vbar[i] = (np.sum(self.v[i+1] * self.tt[i])) / np.sum(self.tt[i])    # eq. 10.1
        
        self.uvbar = np.zeros(23)
        for i in range(23):
            self.uvbar[i] = (np.sum(self.u[i+1] * self.v[i+1] * self.tt[i])) / np.sum(self.tt[i])    # eq. 10.1

    def rms(self):
        '''
        Calculate the RMS velocity u_rms, v_rms and uv_rms
        '''
        u_prime = self.u.copy()
        v_prime = self.v.copy()

        for i in range(23):
            u_prime[i+1] = self.u[i+1] - self.ubar[i+1]
            v_prime[i+1] = self.v[i+1] - self.vbar[i]

        self.u_rms = np.zeros(23)
        for i in range(23):
            self.u_rms[i] = ((np.sum(u_prime[i+1]**2 * self.tt[i])) / np.sum(self.tt[i]))**0.5    # eq. 10.2
        
        self.v_rms = np.zeros(23)
        for i in range(23):
            self.v_rms[i] = ((np.sum(v_prime[i+1]**2 * self.tt[i])) / np.sum(self.tt[i]))**0.5   # eq. 10.2 

        self.uv_rms = np.zeros(23)
        for i in range(23):
            self.uv_rms[i] = ((np.sum(-v_prime[i+1]*u_prime[i+1] * self.tt[i])) / np.sum(self.tt[i]))**0.5    # eq. 10.2
           
    def depth_averaged_velocity(self):
        '''
        Calculate the depth-averaged velocity
        '''
        self.V = 1 / self.h * np.trapezoid(self.ubar, self.y)    # eq. 10.3

    def friction_velocity(self):
        '''
        Calculate the friction velocity
        '''
        A = self.h * self.b
        P = 2 * self.h + self.b
        r_h = A / P

        Re = r_h * self.V / self.nu
        f = 0.0557 / (Re**0.25)    # eq. 10.5
        
        self.u_f = np.sqrt(f/2) * self.V    # eq. 10.4       

    def bounds(self):
        '''
        Calculate the bounds with the friction velocity
        '''
        self.y_plus = self.y * self.u_f / self.nu
        self.Re_tau = self.h * self.u_f / self.nu

        self.upper_bound = int(np.where(self.y_plus <= 0.1*self.Re_tau)[0][-1])
        self.lower_bound = int(np.where(self.y_plus >= 30)[0][0])

    def friction_velocity_calc(self):
        '''
        Calculate the friction velocity with the new bounds
        '''
        A = 2.5
        slope, intercept = np.polyfit(np.log(self.y[self.lower_bound:self.upper_bound]), self.ubar[self.lower_bound:self.upper_bound], 1)
        self.u_f = slope / A    # eq. 3.42

    def van_driest(self):
        '''
        Calculate the van Driest velocity profile
        '''
        kappa = 0.4 # Side 99
        Ad = 25 # Side 99

        den = 1 + (1 + 4 * kappa**2 * self.y_plus**2 * (1 - np.exp(- self.y_plus / Ad))**2)**0.5

        self.u_vd = 2*self.u_f * cumtrapz(1/den, self.y_plus, initial=0)    # eq. 3.108

    def turbulent_kinetic_energy(self):
        '''
        Calculate the turbulent kinetic energy
        '''
        w = 1.8 * self.v_rms**2

        self.k = 1/2 * (self.u_rms**2 + self.v_rms**2 + w)    # eq. 10.6

    def reynolds_stress(self):
        '''
        Calculate the Reynolds stress
        '''
        self.rho = 998 # Density of water [kg/m^3]
        
        tau_bar = self.rho * self.u_f**2 * (1- self.y/self.h)    # eq. 10.8
        self.mu = self.rho * self.nu
        self.rey_stress = tau_bar - self.mu * np.gradient(self.ubar, self.y)    # eq. 10.7

    def energy_production(self):
        '''
        Calculate the energy production
        '''
        self.P = self.rey_stress * np.gradient(self.ubar, self.y)    # p. 701

if __name__ == '__main__':
    if False: # EX1 - Plot the velocity profile
        EX1 = EX()
        EX1.mean_velocity()


        save_as_txt('velocity_profile.txt', EX1.y, EX1.ubar)
        plt.figure()
        plt.plot(EX1.y, EX1.ubar, label='ubar')
        plt.legend()
        plt.xlabel('y [m]')
        plt.ylabel('ubar [m/s]')
        plt.title('Velocity profile')
        plt.show()

    if False: # EX2 - Calculate the depth-averaged velocity
        EX2 = EX()
        EX2.mean_velocity()
        EX2.depth_averaged_velocity()
        print("The depth-averaged velocity is: ", EX2.V, "m/s")

    if False: # EX3 - Friction velocity
        EX3 = EX()
        EX3.mean_velocity()
        EX3.depth_averaged_velocity()
        EX3.friction_velocity()
        print("The friction velocity is: ", EX3.u_f, "m/s")

    if False: # EX4 - Compare new U_f
        EX4 = EX()
        EX4.mean_velocity()
        EX4.depth_averaged_velocity()
        EX4.friction_velocity()
        EX4.bounds()

        save_as_txt('bound.txt', EX4.y, EX4.ubar)

        plt.figure()
        plt.semilogx(EX4.y, EX4.ubar, label='ubar')
        plt.axvline((0.1*EX4.Re_tau) * EX4.nu / EX4.u_f, color='r', linestyle='--', label='Upper bound')
        plt.axvline(30 * EX4.nu / EX4.u_f, color='g', linestyle='--', label='Lower bound')
        plt.legend()
        plt.ylabel('ubar [m/s]')
        plt.xlabel('y [m]')
        plt.title('Velocity profile - Semilog')
        plt.show()

        print("The old friction velocity is: ", EX4.u_f, "m/s")
        print("The old upper bound is: ", (0.1*EX4.Re_tau) * EX4.nu / EX4.u_f) # y [m]
        print("The old lower bound is: ", 30 * EX4.nu / EX4.u_f) # y [m]
        
        
        EX4.friction_velocity_calc()
        print("The new friction velocity is: ", EX4.u_f, "m/s")
        print("The new upper bound is: ", (0.1*EX4.Re_tau) * EX4.nu / EX4.u_f) # y [m]
        print("The new lower bound is: ", 30 * EX4.nu / EX4.u_f) # y [m]

    if False: # EX5 Dimensionless velocity profile
        EX5 = EX()
        EX5.mean_velocity()
        EX5.depth_averaged_velocity()
        EX5.friction_velocity()
        EX5.bounds()
        EX5.friction_velocity_calc()

        save_as_txt('yplusubar.txt', EX5.y_plus, EX5.ubar / EX5.u_f)
        print(EX5.y_plus[EX5.upper_bound])
        plt.figure()
        plt.semilogx(EX5.y_plus, EX5.ubar / EX5.u_f, label='ubar / u_f')
        plt.axvline(5, color='r', linestyle='--', label='y+ = 5')
        plt.axvline(30, color='g', linestyle='--', label='y+ = 30')
        plt.axvline(EX5.y_plus[EX5.upper_bound], color='b', linestyle='--', label='Upper bound')
        plt.legend()
        plt.xlabel('y+')
        plt.ylabel('ubar / u_f')
        plt.title('Dimensionless velocity profile - yplus')
        
        save_as_txt('yhubar.txt', EX5.y/EX5.h, EX5.ubar / EX5.u_f)
        plt.figure()
        plt.semilogx(EX5.y/EX5.h, EX5.ubar / EX5.u_f, label='ubar / u_f')
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('ubar / u_f')
        plt.title('Dimensionless velocity profile - y/h')
        plt.show()

    if False: # EX6 - van Driest velocity profile
        EX6 = EX()
        EX6.mean_velocity()
        EX6.depth_averaged_velocity()
        EX6.friction_velocity()
        EX6.bounds()
        EX6.friction_velocity_calc()
        EX6.van_driest()

        save_as_txt('van_driest.txt', EX6.y_plus, EX6.u_vd/EX6.u_f)
        plt.figure()
        plt.semilogx(EX6.y_plus, EX6.ubar/EX6.u_f, label='ubar')
        plt.semilogx(EX6.y_plus, EX6.u_vd/EX6.u_f, label='van Driest')
        plt.legend()
        plt.xlabel('y+')
        plt.ylabel('ubar/u_f and u_vd/u_f')
        plt.title('Dimensionless velocity profile - yplus')
        plt.show()

    if False: # EX7 - Turbulence data
        EX7 = EX()
        EX7.mean_velocity()
        EX7.rms()
        EX7.depth_averaged_velocity()
        EX7.friction_velocity()
        EX7.bounds()
        EX7.friction_velocity_calc()

        # print("The boundary for log-region is: ", EX7.y_plus[EX7.lower_bound])

        save_as_txt('urms.txt', EX7.y_plus[1:-1], EX7.u_rms / EX7.u_f)

        plt.figure()
        plt.plot(EX7.y_plus[1:-1], EX7.u_rms / EX7.u_f, label='u_rms / u_f')
        plt.xlim(0, 100)
        plt.legend()
        plt.xlabel('y+')
        plt.ylabel('Velocity')
        plt.title('u_rms')

        save_as_txt('vrms.txt', EX7.y_plus[1:-1], EX7.v_rms / EX7.u_f)

        plt.figure()
        plt.plot(EX7.y_plus[1:-1], EX7.v_rms / EX7.u_f, label='v_rms / u_f')
        plt.xlim(0, 100)
        plt.legend()
        plt.xlabel('y+')
        plt.ylabel('Velocity')
        plt.title('v_rms')

        save_as_txt('uvrms.txt', EX7.y_plus[1:-1], EX7.uv_rms / EX7.u_f)

        plt.figure()
        plt.plot(EX7.y_plus[1:-1], EX7.uv_rms / EX7.u_f, label='uv_rms / u_f')
        plt.xlim(0, 100)
        plt.legend()
        plt.xlabel('y+')
        plt.ylabel('Velocity')
        plt.title('uv_rms')
        plt.show()

    if False: # EX8 - Tubulence outer-flow
        EX8 = EX()
        EX8.mean_velocity()
        EX8.rms()
        EX8.depth_averaged_velocity()
        EX8.friction_velocity()
        EX8.bounds()
        EX8.friction_velocity_calc()

        save_as_txt('urms_yh.txt', EX8.y[1:-1]/EX8.h, EX8.u_rms / EX8.u_f)

        plt.figure()
        plt.plot(EX8.y[1:-1]/EX8.h, EX8.u_rms / EX8.u_f, label='u_rms / u_f')
        plt.xlim(0, 1)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Velocity')
        plt.title('u_rms')

        save_as_txt('vrms_yh.txt', EX8.y[1:-1]/EX8.h, EX8.v_rms / EX8.u_f)

        plt.figure()
        plt.plot(EX8.y[1:-1]/EX8.h, EX8.v_rms / EX8.u_f, label='v_rms / u_f')
        plt.xlim(0, 1)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Velocity')
        plt.title('v_rms')

        save_as_txt('uvrms_yh.txt', EX8.y[1:-1]/EX8.h, EX8.uv_rms / EX8.u_f)

        plt.figure()
        plt.plot(EX8.y[1:-1]/EX8.h, EX8.uv_rms / EX8.u_f, label='uv_rms / u_f')
        plt.xlim(0, 1)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Velocity')
        plt.title('uv_rms')
        plt.show()

    if False: # EX9 - Turbulent kinetic energy
        EX9 = EX()
        EX9.mean_velocity()
        EX9.rms()
        EX9.depth_averaged_velocity()
        EX9.friction_velocity()
        EX9.bounds()
        EX9.friction_velocity_calc()
        EX9.turbulent_kinetic_energy()

        save_as_txt('k_yh.txt', EX9.y[1:-1]/EX9.h, EX9.k / EX9.u_f**2)

        plt.figure()
        plt.plot(EX9.y[1:-1]/EX9.h, EX9.k / EX9.u_f**2, label='k / u_f^2')
        plt.xlim(0, 1)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Velocity [m/s]')
        plt.title('Turbulent kinetic energy')
        plt.show()

    if False: # EX10 - Reynolds stress
        EX10 = EX()
        EX10.mean_velocity()
        EX10.rms()
        EX10.depth_averaged_velocity()
        EX10.friction_velocity()
        EX10.bounds()
        EX10.friction_velocity_calc()
        EX10.reynolds_stress()

        save_as_txt('rey_yh.txt', EX10.y/EX10.h, EX10.rey_stress / EX10.u_f**2)
        save_as_txt('rey_compare.txt', EX10.y[1:-1]/EX10.h, EX10.rho * (EX10.uv_rms / EX10.u_f)**2)

        plt.figure()
        plt.plot(EX10.y/EX10.h, EX10.rey_stress / EX10.u_f**2, label='Reynolds stress / u_f^2')
        plt.plot(EX10.y[1:-1]/EX10.h, EX10.rho * (EX10.uv_rms / EX10.u_f)**2, label='uv_rms / u_f')
        plt.xlim(0, 1)
        plt.ylim(0, 1000)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Reynolds stress')
        plt.title('Reynolds stress')
        plt.show()

    if True: # EX11 - Energy production
        EX11 = EX()
        EX11.mean_velocity()
        EX11.rms()
        EX11.depth_averaged_velocity()
        EX11.friction_velocity()
        EX11.bounds()
        EX11.friction_velocity_calc()
        EX11.reynolds_stress()
        EX11.energy_production()

        save_as_txt('P_yh.txt', EX11.y/EX11.h, EX11.P)

        plt.figure()
        plt.plot(EX11.y/EX11.h, EX11.P, label='P')
        plt.xlim(0, 1)
        plt.legend()
        plt.xlabel('y/h')
        plt.ylabel('Energy production')
        plt.title('Energy production')
        plt.show()