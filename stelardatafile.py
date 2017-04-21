import os
import sys
import re
import numpy as np
import pandas as pd
from utils import freq_shift

class StelarDataFile:
    def __init__(self, FileName, PathName):
        self.FileName=FileName
        self.PathName=PathName
        self.options='standard'
        self.datas={}

    def file(self):
        return os.path.join(self.PathName,self.FileName)
                    
    def options(self, text):
        self.options=text
        
    def adddata(self, ie, parameter, data):
        self.datas[ie]=(parameter,data)
        
    def getexp(self,ie):
        return self.datas[ie]

    def getparameter(self,ie):
        parameter, data =self.datas[ie]
        return parameter

    def getdata(self,ie):
        parameter, data = self.datas[ie]
        return data

    def getparvalue(self,ie,par):
        return self.getparameter(ie)[par]

    def getfid(self,ie):
        parameters=self.getparameter(ie)
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

        #calculate series of fid
        fid=pd.DataFrame(self.getdata(ie),index=np.linspace(dw,dw*bs*nblk,bs*nblk),
                         columns=['real', 'im'])/ns
        fid['magnitude']=( fid['real']**2 + fid['im']**2 )**0.5
        return fid

    def get_number_of_experiments(self):
        return len(self.datas)

    def addparameter(self, ie, par, val):
        parameter, data = self.datas[ie]
        parameter[par] = val
        self.adddata(ie,parameter,data)

    def rephase_fids(self):
        for ie in range(1,self.get_number_of_experiments()):
            fid=self.getfid(ie)
            phase=0 #rephase in deg
            freq=self.getparvalue(ie,'F1')
            dt=self.getparvalue(ie,'DW')*1e-6
            
            fid['real']=np.real(
                (fid['real'] + fid['im'] * 1.j) * np.exp(1j*phase/360*2*np.pi ) * np.exp(-1j* freq *np.array(fid.index))
                )
            fid['im']=np.imag(
                (fid['real'] + fid['im'] * 1.j) * np.exp(1j*phase/360*2*np.pi ) * np.exp(-1j* freq *np.array(fid.index))
                )
            #fid['im']=np.imag((fid['real'] + fid['im'] * 1.j) * np.exp(1j*(phase/360+freq*fid.index)*2*np.pi))
            self.addparameter(ie,'rephased_fid',fid)

    def scan_parameters(self,quantiles=20):
        par_df=pd.DataFrame()
        for ppiepigpii in range(1):
            par=[]
            for ie in range(self.get_number_of_experiments()):
                par.append(self.getparameter(ie+1))
            par_df=pd.DataFrame(par)

        # categorize the data
        columns=sorted(par_df.columns)
        discrete = [x for x in columns if (par_df[x].dtype == object or par_df[x].dtype == str)]
        continuous = [x for x in columns if x not in discrete]
        time = [x for x in continuous if x=='TIME']
        quantileable = [x for x in continuous if len(par_df[x].unique()) > quantiles]
        return par_df, columns, discrete, continuous, time, quantileable

    def sdfimport(self):
        ie=1
        olddir=os.getcwd()
        os.chdir(self.PathName)

        with open(self.FileName,"r") as fid:
            line=fid.readline() #skip first line"
            words=['bla']
            parameters=dict()
            while 'DATA' not in words[0]:
                words = fid.readline().split('\t')
                if not words[0]: #probably eof reached
                    print('probably end of file reached')
                    break #escape the while loop
                #read the parameters of file
                if len(words)==2:
                    words[0]=re.sub('[=]','',words[0]) #get rid of '='
                    words[1]=re.sub('[\n]','',words[1])#get rid of trailing \n
                    words[1]=re.sub('[\t]','',words[1])#get rid of \t
                    try:
                        words[1]=int(words[1])                                        #ints are converted
                        exec('parameters[\''+words[0].rstrip()+'\'] = '+str(words[1]))  #and stored as int
                    except ValueError:
                        try:
                            words[1]=float(words[1])                                        #floats are converted
                            exec('parameters[\''+words[0].rstrip()+'\'] = '+str(words[1]))  #and stored as float
                        except ValueError:
                            exec('parameters[\''+words[0].rstrip()+'\'] = \''+words[1]+'\'')#else its stored as string
                    if words[0]=='TIME':
                        try:
                            parameters['TIME']=pd.to_datetime(words[1])
                        except:
                            print('TIME is not in datetime Format')
                        

                try:
                     if 'DATA' in words[0]:
                        data=[]
                        if parameters['NBLK']==0:
                            parameters.update({'NBLK': 1})
                        for i in range(0,int(parameters['NBLK']*parameters['BS'])):
                            columns=fid.readline().replace('\n','').split('\t')
                            data.append([])
                            for c, ii in zip(columns, range(0,10)):
                                data[i].append(int(c))
                        self.adddata(ie,parameters,data)
                        #print(self.getexp(ie))
                        #print(ie)
                        ie = ie + 1
                        words[0]='bla'
                        parameters=dict()
                except KeyError:
                    print('fastforward') #no experiment found
        print(str(ie-1)+' experiments read')
        os.chdir(olddir)
