import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import stelardatafile as sdf

path=os.path.join(os.path.curdir,'data')
polymer=sdf.StelarDataFile('297K.sdf',path)
polymer.sdfimport()

### preprocessing of the data experiment number ie=1
ie=1
parameters=polymer.getparameter(ie)
bs=int(parameters['BS'])
try:
    nblk=int(parameters['NBLK'])
except:
    nblk=1;
ns=int(parameters['NS'])
try:
    dw=parameters['DW']*1e-6 #dwell time is in [mu s]
except:
    dw=1
temperature=parameters['TEMP']
fid=pd.DataFrame(polymer.getdata(ie),index=np.linspace(dw,dw*bs*nblk,bs*nblk),
                 columns=['real', 'im'])/ns
fid['abs']=( fid['real']**2 + fid['im']**2 )**0.5 # last two lines may represent a seperate method
fid.plot(marker='x',linestyle='none') #better: bokeh server
#plt.show() #only for illustration use, better: bokeh server


#now consider a relaxation experiment with relaxation time tau_decay~T1MX:


#first we need to create the tau (time) axis
###not working yet
##if ('NP/S' or 'PP/S') in parameters['EXP']:
##    try:
##        T1MX=parameters['T1MX'] #capital letters needed for exec below
##    except KeyError:
##        print(['no t1mx found in experiment number ', ie, '!! no relaxation experiment?'])
##    if parameters['BGRD']=='LOG':
##        bini=exec(parameters['BINI']) #usually this will read sth. like '(0.3*T1MX)'
##        bend=exec(parameters['BEND']) #here also
##        tau=np.logspace(bini,bend,nblk)
##    if parameters['BGRD']=='LIN':
##        bini=exec(parameters['BINI']) #usually this will read sth. like '(0.3*T1MX)'
##        bend=exec(parameters['BEND']) #here also
##        tau=np.linspace(bini,bend,nblk)
##    #this listformat is so inhumane it should be changed by STELAR
##    if parameters['BGRD']=='LIST':
##        str=parameters['BLST']
##        str=str.replace(';',':')
##        sep=str.split(':')
##        print('list not supported right now')

tau=np.logspace(-3,np.log10(5*parameters['T1MX']),nblk) #as a dummy

#second we need to create the magnetisation (phi) axis (essentially the amplitude of the fid)

startpoint=int(0.05*bs)-1
endpoint=int(0.1*bs)
phi=np.zeros(nblk)
for blk in range(nblk):
    start=startpoint + blk * bs
    end=endpoint + blk * bs
    phi[blk]=fid['abs'].iloc[start:end].sum() / (endpoint-startpoint)

plt.figure()
plt.plot(tau,phi)
plt.show()





