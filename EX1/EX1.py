import scipy.io as sp
import numpy as np
import matplotlib.pyplot as plt

def save_as_txt(filename: str, x_data, y_data) -> None:
    '''
    Save the data as a .txt file
    '''
    np.savetxt(filename, np.column_stack((x_data, y_data)))
    


class EX1:
    def __init__(self):
        '''
        Starts with loading the data from the .mat file
        '''
        self.load_data()

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
        self.u_rms = np.zeros(23)
        for i in range(23):
            self.u_rms[i] = np.sqrt((np.sum(self.u[i+1]**2 * self.tt[i])) / np.sum(self.tt[i]))    # eq. 10.2
        
        self.v_rms = np.zeros(23)
        for i in range(23):
            self.v_rms[i] = np.sqrt((np.sum(self.v[i+1]**2 * self.tt[i])) / np.sum(self.tt[i]))    # eq. 10.2 

        self.uv_rms = np.zeros(23)
        for i in range(23):
            self.uv_rms[i] = np.sqrt((np.sum(self.u[i+1]*self.v[i+1] * self.tt[i])) / np.sum(self.tt[i]))    # eq. 10.2
           
        


if __name__ == '__main__':
    if False: # Test
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 2, 3, 4, 5])
        save_as_txt('test.txt', x, y)

    if True: # EX1
        EX1 = EX1()
        EX1.mean_velocity()

        plt.figure()
        plt.plot(EX1.y, EX1.ubar, label='ubar')
        plt.legend()
        plt.xlabel('y [m]')
        plt.ylabel('ubar [m/s]')
        plt.title('Velocity profile')
        plt.show()
    

    print("Stop")
