Sep_interval = 27
Njob=20
import numpy as np
import math
import numpy.polynomial.legendre as lgd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
filename = 'BinHist'
dir_name = 'weighed_uniform'
filedir = '/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/BinHist/'+dir_name+'/'
plot_title = 'Obiwan VS eboss NGC'
def JKnife_CorrFunc(Njob, Jacknife=-1,k0=1, order=0,name = filename):
    #    pp = PdfPages('Correlation_Function_of_order_No'+str(order)+'.pdf')
    #    b2=np.loadtxt('xi0isoChallenge_matterpower6.0.dat').transpose()
    #    plt.plot(b2[0],b2[1]*b2[0]*b2[0]*k0)
    filenameDD=[]
    filenameDR=[]
    filenameRR=[]
    for i in range(0,Njob):
        for j in range(i,Njob):
            filenameDD.append(filedir+name+'D'+str(i)+'D'+str(j)+'.dat')
    for i in range(0,Njob):
        for j in range(0,Njob):
                filenameDR.append(filedir+name+'D'+str(i)+'R'+str(j)+'.dat')
    for i in range(0,Njob):
        for j in range(i,Njob):
            filenameRR.append(filedir+name+'R'+str(i)+'R'+str(j)+'.dat')

    FilelistDD = []
    FilelistRR = []
    FilelistDR = []
    COUNT = 0
    void_pair = np.zeros((Sep_interval, 2)) #interval of angular scales
    for i in range(0,Njob):
        for j in range(i,Njob):
            if i!=Jacknife and j!=Jacknife :
                aDD = np.loadtxt(filenameDD[COUNT])
                FilelistDD.append(aDD)
            else:
                FilelistDD.append(void_pair)
            COUNT+=1
    COUNT = 0
    for i in range(0,Njob):
        for j in range(0,Njob):
            if i!=Jacknife and j!=Jacknife :
                aDR = np.loadtxt(filenameDR[COUNT])
                FilelistDR.append(aDR)
            else:
                FilelistDR.append(void_pair)
            COUNT+=1
    COUNT = 0
    for i in range(0,Njob):
        for j in range(i,Njob):
            if i!=Jacknife and j!=Jacknife :
                aRR = np.loadtxt(filenameRR[COUNT])
                FilelistRR.append(aRR)
            else:
                FilelistRR.append(void_pair)
            COUNT+=1
    DD_total=np.zeros(Sep_interval)
    DR_total=np.zeros(Sep_interval)
    RR_total=np.zeros(Sep_interval)
    for i in range(0,(Njob+1)*Njob/2):
        for j in range(0,len(aDD)):
            DD_total[j]+=FilelistDD[i][j][1]
    for i in range(0,Njob*Njob):
        for j in range(0,len(aDD)):
                DR_total[j]+=FilelistDR[i][j][1]
    for i in range(0,Njob*(Njob+1)/2):
        for j in range(0,len(aDD)):
                RR_total[j]+=FilelistRR[i][j][1]
   
    TotalPoints=np.loadtxt(filedir+'TotalPoints.txt')
    DD_total_num=0
    DR_total_num=0
    RR_total_num=0
    for i in range(0,len(TotalPoints)):
        if TotalPoints[i][3]!=Jacknife and TotalPoints[i][4]!=Jacknife :
            DD_total_num+=TotalPoints[i][0]
            DR_total_num+=TotalPoints[i][1]
            RR_total_num+=TotalPoints[i][2]
    #print 'DD_total_num='+str(DD_total_num)+' DR_total_num'+str(DR_total_num)+' RR_total_num'+str(RR_total_num)
    Final_total=np.zeros(Sep_interval)
    for i in range(0,Sep_interval):
             Final_total[i] = ((DD_total[i]/DD_total_num)-(DR_total[i]/DR_total_num)*2+(RR_total[i]/RR_total_num))/(RR_total[i]/RR_total_num)
    return Final_total

def CorrFunc(Njob=Njob,k0=1, order=0,name_corr = filename):
    Final_total = []
    Final_total = JKnife_CorrFunc(Njob,-1,k0,order,name = name_corr)

    Jacknife_list = []
    for i in range(0,Njob):
        Jacknife_list.append(JKnife_CorrFunc(Njob,i,k0,order,name = name_corr))
        
    CorrFunc_Err = np.zeros(Sep_interval)
    for i in range(0,Sep_interval):
        for j in range(0,Njob):
            CorrFunc_Err[i] += (Jacknife_list[j][i]-Final_total[i])**2
        CorrFunc_Err[i] = CorrFunc_Err[i]*(Njob-1)/Njob
    CorrFunc_Err = np.sqrt(CorrFunc_Err)
    from math import log10
    from math import pi
    x = np.zeros(Sep_interval)
    MaxAng = log10(1.*pi/180.)
    MinAng = log10(0.01*pi/180.)
    for i in range(0,Sep_interval):
        x[i] = (MaxAng-MinAng)*(i+0.5)/Sep_interval+MinAng
    x = np.zeros(Sep_interval)
    from math import pi
    from math import pow
    MaxAngle = log10(5.*pi/180.)
    MinAngle = log10(0.01*pi/180.)
    AngleInteval = (MaxAngle-MinAngle)/(Sep_interval-1.)
    for i in range(0,Sep_interval):
        x[i]=pow(10,i*AngleInteval+MinAngle)*180./pi
    
    np.savetxt('/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/'+dir_name+'.out',zip(x,Final_total,CorrFunc_Err))
    #plt.errorbar(x, Final_total, yerr=CorrFunc_Err)
    plt.gca().set_xscale("log", nonposx='clip')
    Obiwan = plt.errorbar(x, Final_total,CorrFunc_Err,marker = 'o',label = 'Obiwan')
    plt.title(plot_title)
    plt.xlabel(r'$\theta$'+'(Degree)')
    plt.ylabel('w('+r'$\theta$'+')')
    
    #compare with power law:
    pwlx = np.linspace(0.02,1)
    pwly = 1.*0.01*np.power(pwlx,-0.8)
    powerlaw, = plt.plot(pwlx,pwly,label = 'power law')
    
    
    #compare with Johan's result:
    johan_random = np.loadtxt('/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/rawdata/corrfunc_data/Johan_random.txt').transpose()
    Johan = plt.errorbar(johan_random[0],johan_random[1],johan_random[2],label = 'Johan')
    
    #plt.gca().set_xlim((0.009, 1))
    #plt.gca().set_ylim((-0.3, 0.25))
    plt.legend(handles = [Obiwan,powerlaw,Johan])

    
    plt.show()
    return True

CorrFunc()
'''    
def CorrFunc_Add(Njob=20,k0=1, order=0):
    b2=np.loadtxt('xi0isoChallenge_matterpower6.0.dat').transpose()
    plt.plot(b2[0],b2[1]*b2[0]*b2[0]*k0)
    plt.xlabel('Mpc',size=16)
    CorrFunc(Njob,name_corr = 'NewData_subfiles_chunk21')
    CorrFunc(Njob,name_corr = 'NewData_subfiles_chunk22')
    CorrFunc(Njob,name_corr = 'NewData_subfiles_chunk23')
    plt.show()     

def CorrFunc_ALL_sub(index,Njob=20,k0=1,order=0):
    Final_total = []
    Final_total = JKnife_CorrFunc(Njob,index,k0,order)
    d=[0]*len(Final_total[0])
    d[order]+=1
    b=[0]*len(Final_total)
    for i in range(0,len(Final_total)): 
        for j in range(0,len(Final_total[0])):
             b[i]=b[i]+0.01*Final_total[i][j]*lgd.legval(Final_total[i][j],d)
        b[i]=(2*order+1)*b[i]
    c=[(i+0.5)*5 for i in range(0,len(Final_total))]
    for i in range (0,len(Final_total)):
        b[i]=b[i]*c[i]*c[i]
    plt.plot(c,b)
    return True

def CorrFunc_ALL(Njob=20):
    for i in range(0,20):
        CorrFunc_ALL_sub(i,Njob)
    plt.show()
    return True

def JKnife_show(index):
    import sys
    CorrFunc_ALL_sub(index)
    plt.show()
    print 'continue?'
    sys.stdin.readline()
    plt.clf()
    return True

def weighed_tot(k0=1):
    chunk21 = np.loadtxt('./data/NewData_subfiles_chunk21_wCut.3_upweight.txt').transpose()
    chunk22 = np.loadtxt('./data/NewData_subfiles_chunk22_wCut.3_upweight.txt').transpose()
    chunk23 = np.loadtxt('./data/NewData_subfiles_chunk23_wCut.3_upweight.txt').transpose()
    chunk_tot = [0]*len(chunk21[0])
    chunk_err_tot = [0]*len(chunk21[0])
    for i in range(0,len(chunk_tot)):
        chunk_tot[i]+=chunk21[1][i]/chunk21[2][i]**(2)
        chunk_tot[i]+=chunk22[1][i]/chunk22[2][i]**(2)
        chunk_tot[i]+=chunk23[1][i]/chunk23[2][i]**(2)
        chunk_tot[i]=chunk_tot[i]/(chunk21[2][i]**(-2)+chunk22[2][i]**(-2)+chunk23[2][i]**(-2))
        chunk_err_tot[i] = (chunk21[2][i]**(-2)+chunk22[2][i]**(-2)+chunk23[2][i]**(-2))**(-0.5)
    p21 = plt.errorbar(chunk21[0],chunk21[1],yerr=chunk21[2])
    p22 = plt.errorbar(chunk22[0],chunk22[1],yerr=chunk22[2])
    p23 = plt.errorbar(chunk23[0],chunk23[1],yerr=chunk23[2])
    plt.axis([0,200,-100,100])
    ptot = plt.errorbar(chunk21[0],chunk_tot,yerr=chunk_err_tot)
    b2 = np.loadtxt('xi0isoChallenge_matterpower6.0.dat').transpose()
    ptheory = plt.plot(b2[0],b2[1]*b2[0]*b2[0]*k0)
    plt.legend((p21, p22,p23,ptot,ptheory), ('chunk21', 'chunk22','chunk23','total','theory'))
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages('./data/output/combined_wCut.3') as pdf:
         plt.title('combined correlation function for three chunks')
         plt.xlabel('Mpc')
         pdf.savefig()
         plt.show()
    f = open('./data/combined_wCut.3_upweight.txt','w')
    for i in range(0,len(chunk21[0])):
        f.write(str(chunk21[0][i])+' '+str(chunk_tot[i])+' '+str(chunk_err_tot[i])+'\n')
    return True

if __name__=='__main__':
     chunk21 = np.loadtxt('./data/NewData_subfiles_chunk21_wCut.3.txt').transpose()
     chunk22 = np.loadtxt('./data/NewData_subfiles_chunk22_wCut.3.txt').transpose()
     chunk23 = np.loadtxt('./data/NewData_subfiles_chunk23_wCut.3.txt').transpose()
     chunkcb = np.loadtxt('./data/combined_wCut.3.txt').transpose()
     plt.plot(chunk21[0],chunk21[1],color = 'red',alpha=0.3)
     plt.fill_between(chunk21[0],chunk21[1]-chunk21[2],chunk21[1]+chunk21[2],color = 'salmon',alpha=0.3)
     plt.plot(chunk22[0],chunk22[1],color = 'green',alpha=0.3)
     plt.fill_between(chunk22[0],chunk22[1]-chunk22[2],chunk22[1]+chunk22[2],color = 'lime',alpha=0.3)
     plt.plot(chunk23[0],chunk23[1],color = 'mediumblue',alpha=0.3)
     plt.fill_between(chunk23[0],chunk23[1]-chunk23[2],chunk23[1]+chunk23[2],color = 'blue',alpha=0.3)
     plt.plot(chunkcb[0],chunkcb[1],color = 'black')
     plt.fill_between(chunkcb[0],chunkcb[1]-chunkcb[2],chunkcb[1]+chunkcb[2],color = 'k',alpha=0.6)
     import matplotlib.patches as mpatches
     red_patch = mpatches.Patch(color='red', label='chunk21')
     green_patch = mpatches.Patch(color='green', label='chunk22')
     blue_patch  = mpatches.Patch(color='blue', label= 'chunk23')
     black_patch = mpatches.Patch(color='black',label='combined')
     plt.legend(handles=[red_patch,green_patch,blue_patch,black_patch])
     from matplotlib.backends.backend_pdf import PdfPages
     with PdfPages('./data/CorrelationFunction.pdf') as pdf:
          #plt.rc('xtick', labelsize=26)
          #plt.rc('ytick', labelsize=16)
          #plt.rc('axes', titlesize=12)
          #plt.rc('axes', labelsize=16)  
          plt.tick_params(axis='x', labelsize=14)
          plt.tick_params(axis='y', labelsize=14)
          plt.title('Correlation Function of three chunks and their combination')
          plt.xlabel('Mpc')
          plt.axis([0,200,-50,100])
          pdf.savefig()
          #plt.show()
          #pdf.close()
'''
