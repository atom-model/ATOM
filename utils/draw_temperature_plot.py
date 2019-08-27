import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def draw_temperature_plot(lon=180, data_dir='../benchmark/output/', output_dir='../benchmark/output/'):
    plt.figure(figsize=(20,10))
    plt.xlim(90, -90)
    plt.xticks(range(90,-91,-10))
    plt.title('Temperature along longitude 180 at different times')

    x=range(90,-91,-1)

    for time in range(0,101,20):
        adata = np.loadtxt(data_dir + '[{}Ma_smooth.xyz]_PlotData_Atm.xyz'.format(time), skiprows=1)
        plt.plot(x, adata[:, 6].reshape((361,181))[180],  label='time: {} Ma'.format(time), linewidth=2)
    
    plt.legend(title='Times', loc="lower center")
    plt.xlabel('Latitude')
    plt.ylabel('Temperature(Celsius)')
    plt.savefig(output_dir + 'temperature_at_lon_180.png')
    #plt.show()

if __name__ == "__main__":
    draw_temperature_plot()

