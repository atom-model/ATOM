#!/usr/bin/python
# -*- coding: utf-8 -*-
# draw quiver
import os, sys, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def down_sample(a):
    aa = a.reshape((361,181))
    tmp=[]
    for t in aa:
        tmp.append(t[3::4])
        #print t
    #print tmp
    return tmp[3::4]

def fix_wrong_lats(a):
    aa = a.reshape((361,181))
    tmp = []
    for r in aa:
        r0=np.fliplr(r[0:90].reshape((1,90)))
        #print r0
        r1=np.fliplr(r[90:181].reshape((1,91)))
        #print r0[0]
        tmp.append(np.concatenate((r1[0],r0[0])))
    return np.array(tmp)

def draw_velocity(time, output_dir, data_dir, topo_suffix, atm = 'Atm'):
    data = np.genfromtxt(data_dir + '/[{0}Ma_{1}.xyz]_PlotData_{2}.xyz'.format(time, topo_suffix, atm), skip_header=1)
    for d in data:
        d[1]=90-d[1]
    x = data[:,0]
    y = data[:,1]
    z = data[:,6]

    #vx=np.sin(data[:,3])*np.cos(data[:,4])
    #vy=np.sin(data[:,3])*np.sin(data[:,4])
    vx=-data[:,3]
    vy=data[:,4]
    magitude = data[:,5]
    topo = data[:,2]

    plt.figure(figsize=(15, 8))

    m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='cyl', lon_0=0)  

    xx = x
    x, topo = m.shiftdata(xx, datain = topo, lon_0=0)
    x, z = m.shiftdata(xx, datain = z, lon_0=0)
    x, vx = m.shiftdata(xx, datain = vx, lon_0=0)
    x, vy = m.shiftdata(xx, datain = vy, lon_0=0)
    x, magitude = m.shiftdata(xx, datain = magitude, lon_0=0)

    xi, yi = m(x, y)
    #cs = m.scatter(xi, yi, marker='.', c=z, alpha=0.5, lw=0)

    #m.quiver(xi[::39],yi[::39],vx[::39],vy[::39])
    #vx = fix_wrong_lats(vx)
    #vy = fix_wrong_lats(vy)
    cs = m.quiver(down_sample(xi), down_sample(yi), down_sample(vy), down_sample(vx), down_sample(magitude), width=0.001,
             headlength=7, headwidth=5, pivot='tail', clim=[0,1.5])
    #m.scatter(down_sample(xi),down_sample(yi),marker='.',color='r')
    
    m.contour( xi.reshape((361,181)), yi.reshape((361,181)), topo.reshape((361,181)),
                            colors ='k', linewidths= 0.3 )

    m.drawparallels(np.arange(-90., 90., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 180., 45.), labels=[0,0,0,1], fontsize=10)
    #m.drawcoastlines()   

    cbar = m.colorbar(cs, location='bottom', pad="10%")
    plt.title("{}Ma Velocity".format(time))
    plt.savefig(output_dir + '/{}Ma_velocity.png'.format(time), bbox_inches='tight')
    print(output_dir + '/{}Ma_velocity.png has been saved!'.format(time))
    #plt.show()
    plt.close()

if  __name__ == "__main__":
    start_time = 0
    end_time  = 5
    time_step = 5
    data_dir = "../benchmark/output"
    topo_suffix = 'Golonka'
    output_dir = "./velocity"
    try:
        start_time = int(sys.argv[1])
        end_time  = int(sys.argv[2])
        time_step = int(sys.argv[3])
        data_dir = sys.argv[4]
        topo_suffix = sys.argv[5]
    except:
        print("Usage 'python "+ sys.argv[0] + " 0 100 5 ../benchmark/output Golonka'")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for time in range(start_time, end_time+1, time_step):
        draw_velocity(time, output_dir, data_dir, topo_suffix)
