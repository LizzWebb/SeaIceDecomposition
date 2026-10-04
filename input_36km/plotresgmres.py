#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
######################## -*- coding: utf-8 -*-
"""Usage: plotres.py INPUTFILE(S)
"""
import sys
from getopt import gnu_getopt as getopt
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
#matplotlib.use("TkAgg")

# Using LaTeX in figures
#plt.rc('text', usetex=True)
#plt.rc('font', family='sans')
#plt.rc('font', size='14')

# parse command-line arguments
try:
    # looking to see the short and long versions of the command-line options/arguments
    optlist,args = getopt(sys.argv[1:], 'hf:n:i:m:', ['help','files=','names=','fig_name=','maxNLit=','verbose'])
except getopt.GetoptError:
    # If error, print the instructions line 
    print('Use: plotresgmres -f <files> -n <names for plot> -i <figure name> -m <max NonLin its>')

# What was entered from the command-line options/arguments    
print('optlist= ',optlist)
# Whatever was entered but did not match one of the above defined arguments 
print('args= ',args)

# Initializing defaults
files = []
names = []
fig_name = ''
maxNLit = 1

# If no arguments/options are specified, then assume the command line is a list of files to plot
# Ex if you wrote python thiscript.py STDOUT1 STOUT2
if len(optlist)==0:
    files=args
# If you used only the options in the command line 
elif len(args)==0:
    # Loops through all the options and prints the settings used/applied 
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            print('plotresgmres -f <files> -n <names> -i <figure name> -m <max NonLin its>')
            sys.exit()
        elif opt in ('-f','--files'):
            files.append(str(arg))
        elif opt in ('-n','--names'):
            names.append(str(arg))
        elif opt in ('-i','--fig_name'):
            fig_name = str(arg)
        elif opt in ('-m','--maxNLit'):
            maxNLit = int(arg)
        else:
            print('unknown option')
            
# If not file is given/defined, then it assumed STDOUT.0000        
if len(files)==0:
    files=['STDOUT.0000']

print(files)
print(names)
print(fig_name)

def get_output (fname, mystring):
    """parse fname and get some numbers out"""
    import re
    r = re.compile(r'SEAICEnonLinIterMax')
    iters = []
    itersLSR = []
    res   = []
    lsr   = 0
    nMax  = 1
    tempLinLSR = 0
    try:
        f=open(fname)
    except:
        print(fname + " does not exist, continuing")
    else:
        nextline = False
        for line in f:

            # if nextline:
            #     ll = line.split()
            #     nMax = int(ll[-1].replace(',',''))
            #     nextline = False

            # if r.search(line):
            #     nextline = True
            nMax = 20

            if mystring in line:
                ll = line.split()
                if 'gamma_lin' in line and 'gamma_lin' in mystring:
                    res.append(float(ll[-1].replace('D','e').replace(',','')))
                    iters.append(int(ll[-3].replace(',','')))
                elif 'FGMRES' in line:
                    res.append(float(ll[-1].replace('D','e').replace(',','')))
                    iters.append(int(ll[-7].replace(',','')))
            elif 'SEAICE_LSR: Residual Initial ipass,Uice,Vice' in line:
                if 'FGMRES' in mystring:
                    pass
                else:
                    lsr = 1
                    ll = line.split()
                    ures = float(ll[-1].replace('D','e').replace(',',''))
                    vres = float(ll[-2].replace('D','e').replace(',',''))
                    res.append(np.sqrt(ures*ures+vres*vres))
                    iters.append(int(ll[-3].replace(',','')))
                    LSRlinIter = True
            elif 'SEAICE_LSR (ipass=' in line and LSRlinIter == True:
                    ll = line.split()
                    # print(ll)
                    if 'dU' in line:
                        tempLinLSR = int(ll[-3])
                    elif 'dV' in line:
                        tempLinLSR = tempLinLSR + int(ll[-3])
                        itersLSR.append(tempLinLSR)
                        tempLinLSR = 0
                        LSRlinIter = False
            elif 'SEAICE_JFNK: Newton iterate / total, JFNKgamma_lin' in line \
                and 'gamma_lin' in mystring:
                    ll = line.split()
                    res.append(float(ll[-1].replace('D','e').replace(',','')))
                    iters.append(int(ll[-3].replace(',','')))
            elif 'SEAICE_JFNK: Newton iterate ' in line \
                 and 'Nb. of FGMRES iterations' in line \
                 and 'gamma_lin' not in mystring:
                    ll = line.split()
                    res.append(float(ll[-1].replace('D','e').replace(',','')))
                    iters.append(int(ll[-7].replace(',','')))

        if lsr==1:
            iters = np.asarray(iters)
            for it in range(1,len(iters)):
                if iters[it-1]>iters[it]:
                    maxits = iters[it-1]-iters[it]
                    iters[it:] = iters[it:]+maxits

        f.close()

    return iters, res, nMax, itersLSR
# done


# everything OK, lets start
fig, ax =plt.subplots(2,1,sharex=True) #,sharey=True)
newtonstr='SEAICE_KRYLOV: Picard iterate / total, KRYLOVgamma_lin, initial norm'
fgmresstr='SEAICE_KRYLOV: Picard iterate / total'
#fgmresstr = 'Nb. of FGMRES iterations'
lstr = ['.-','x-','+-']
i, j = 0, -1

for infile in files:
    if np.mod(i,7)==0: j = j + 1
    i = i + 1
    # get the data
    iters, fres, nMax, LSRlinIter = get_output(infile, newtonstr)
    fgmiters, fgmres, nfMax, itersLSR = get_output(infile, fgmresstr)
    nMax, nfMax = maxNLit, maxNLit

    if len(fgmiters)!=len(itersLSR): # in case we wait for the LSR
        fgmiters=fgmiters[:-1]

    if len(names)==len(files):
        lab = str(names[i-1]+' - tot LSR its='+str(np.sum(itersLSR)))
        #lab = str(names[i-1])
    else:
        print(infile,' Total number of linear iterations: ',np.sum(itersLSR))
        lab = infile.split('/')[-1]
        if len(infile.split('/'))>1:
            lab = infile.split('/')[-2]
    fres = np.asarray(fres) #/fres[0]
    if len(iters)>0:
        # now plot everything
        ax[0].semilogy(np.asarray(iters)/float(nMax), fres, lstr[j],linewidth=1.0, label = lab)
        if len(itersLSR)>0:
            ax[1].semilogy(np.asarray(fgmiters)/float(nfMax), itersLSR, lstr[j], linewidth=1.0, label = lab)
        else:
            ax[1].plot(np.asarray(fgmiters)/float(nfMax), fgmres, lstr[j], linewidth=1.0, label = lab)

ax[0].set_ylabel('scaled residual')
#ax.set_title('JFNK')
plt.legend(loc = 7, markerscale=2, prop={'size': 12})
if len(itersLSR)>0:
    ax[1].set_ylabel('LSR linear iterations')
else:
    ax[1].set_ylabel('FMGRES (linear) iterations')

if nMax > 1:
    ax[1].set_xlabel('time step')
else:
    ax[1].set_xlabel('Nonlinear iterations')

for axitem in ax:
#    axitem.set_xlim([35400,36000])
#    axitem.set_xlim([69000,70200])
#    axitem.set_xlim([36000,38000])
#    axitem.set_xlim([0,100])
#    axitem.set_xlim([100,420])
#    axitem.set_xlim([np.max([0,iters[-1]-10*100]),iters[-1]])
    axitem.grid(True)

plt.show()

if len(fig_name)==0:
    figname = 'res'
else:
    figname = fig_name

plt.tight_layout()
fig.savefig(figname,dpi=600,bbox_inches='tight')
print('figure saved with name:',fig_name)

