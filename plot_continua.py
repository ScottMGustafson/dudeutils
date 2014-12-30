import dudeutils
from dump_parser import FitData
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import lineid_plot
import os
import numpy as np



#some defaults
dbase = 'database.xml'
cut=[4845.,4853.]
colors=[]
multiplier = 10.E12
output_fname='continua.png'
color_cycle=['r','b','g','y','c','m','k']



def plot_cut(lst,xr=None):

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=14)
    fig,ax = plt.subplots(figsize=(9.,7.))

    x=lst[0].waves

    fits = [multiplier*item.abs for item in lst]
    conts= [multiplier*item.cont for item in lst]
    ax.plot(x,multiplier*lst[0].flux, '-k',label='data',linestyle='steps-')
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    for i in range(0,len(conts)):
        ax.plot(x, conts[i],linestyle='default',color='black')
        ax.plot(x, fits[i],linestyle='dashed',color='black')

    if not xr is None:
        ax.set_xlim([xr[0],xr[-1]])
    ax.set_xlabel(r"Wavelength (\r{A})")

    ax.set_ylim(0.,0.65)
    ax.set_yticks(np.arange(0.1,0.8,0.1))
    plt.ylabel(r"$F_{\lambda}~(10^{12}~$erg~s$^{-1}~$cm$^{-2}~$\r{A}$^{-1})$")

    line_wave = [4847.28, 4847.52, 4848.60]
    line_label1 = ['D I','H I','H I']
    label1_size = np.array([12, 12, 12])

    #fig, ax = lineid_plot.plot_line_ids(x, lst[0].flux, line_wave, line_label1, label1_size,fig=fig,ax=ax)

    plt.grid(False)
    plt.subplots_adjust(top=0.87,bottom=0.13,left=0.1,right=0.95)
    plt.savefig(output_fname)

    return

if __name__ == '__main__':
    #dumpfiles = [os.path.splitext(item[-1])[0]+'.dump' for item in dudeutils.all_conts(dbase)]
    dumpfiles = [os.path.splitext(item[-1])[0]+'.dump' for item in dudeutils.cont_check_pipeline()]
    print(len(dumpfiles))
    #waves = FitData(dumpfiles[0]).waves
    #ind = FitData.get_ind(waves=waves,start=cut[0],end=cut[1])
    #continua = [ FitData(item).get_cut(indices=ind) for item in dumpfiles]
    continua = [ FitData(item) for item in dumpfiles]
    plot_cut(continua,xr=[cut[0],cut[1]])

    

