import os
import sys
import re
import numpy as np
import numpy.matlib
import pandas as pd
from num_string_eval import NumericStringParser
from utils import freq_shift
import logging

class StelarDataFile:
    def __init__(self, FileName, PathName, logger=None):
        self.logger = logger or logging.getLogger(__name__)
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
        fid['magnitude'] = ( fid['real']**2 + fid['im']**2 )**0.5
##        complex_fid = fid['real']+ 1j * fid['im']
##        complex_fid_corr=np.mean(complex_fid)
##        blocks=[]
##        for g,df in fid.groupby(np.arange(len(fid))//bs):
##            blocks.append(df)
##        for block in blocks:
##            rephase
        return fid

    def get_number_of_experiments(self):
        return len(self.datas)
    
    def get_tau_axis(self, ie):
        parameters = self.getparameter(ie)
        nsp = NumericStringParser()
        nblk = parameters['NBLK']
        T1MX = parameters['T1MX'] # T1MX is used in the 'eval' expressions below
        # TODO: not tested yet for all cases
        if parameters['BGRD'] == 'LIST':
            temp = parameters['BLST']
            temp.replace(';', ':')
            sep_indices = [pos for pos, char in enumerate(temp) if char == ':'] # find indices of ':'
            Tini = nsp.eval(temp[:sep_indices[0]].replace('T1MX', str(T1MX)))
            Tend = nsp.eval(temp[sep_indices[0]+1:sep_indices[1]].replace('T1MX', str(T1MX)))
            npts = nsp.eval(temp[sep_indices[2]+1:].replace('T1MX', str(T1MX))) # number of points selected: can be ~= NBLK    
            if temp[sep_indices[1]+1:sep_indices[2]] == 'LIN':
                listx = np.linspace(Tini,Tend,npts);
            elif temp[sep_indices[1]+1:sep_indices[2]] == 'LOG':
                listx = np.logspace(np.log10(Tini),np.log10(Tend),npts);
            nrep = np.ceil(nblk/npts) # find if the time vector needs to be longer
            x = numpy.matlib.repmat(listx,1,nrep) # re-create the time vector
            x = x[:nblk] # select the portion corresponding to the number of blocs (needed if npts~=nblk)
        elif parameters['BGRD'] == 'LIN':
            Tini = nsp.eval(parameters['BINI'].replace('T1MX', str(T1MX)))
            Tend = nsp.eval(parameters['BEND'].replace('T1MX', str(T1MX)))
            x = np.linspace(Tini, Tend, nblk) # re-create the time vector
        elif parameters['BGRD'] == 'LOG':
            Tini_Tend = [0, 0]
            # This has to be so complicated b/c of mixed datatype that is possible for parameters['BGRD'].
            # Documentation would be beneficial, what can occur
            for nb, b in enumerate(['BINI', 'BEND']):
                if type(parameters[b]) == type(0.0):
                    Tini_Tend[nb] = parameters[b]
                elif type(parameters[b]) == type('asdf'):
                    if 'T1MX' in parameters[b]:
                        Tini_Tend[nb] = nsp.eval(parameters[b].replace('T1MX', str(T1MX)))
                    else:
                        Tini_Tend[nb] = nsp.eval(parameters[b])
            Tini, Tend = Tini_Tend
            x = np.logspace(np.log10(Tini),np.log10(Tend),nblk) # re-create the time vector
        return x

    def addparameter(self, ie, par, val):
        parameter, data = self.datas[ie]
        parameter[par] = val
        self.adddata(ie,parameter,data)

    def rephase_fids(self):
        self.logger.debug('called rephase_fids')
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
                    self.logger.info('probably end of file reached')
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
                            self.logger.warning('TIME is not in datetime Format')
                try:
                    if 'DATA' in words[0]:
                        data=[]
                        if parameters['NBLK']==0:
                            parameters.update({'NBLK': 1})
                        for i in range(0,int(parameters['NBLK']*parameters['BS'])):
                            columns=fid.readline().replace('\n','').split('\t')
                            data.append([])
                            for c in columns:
                                data[i].append(int(c))
                        self.adddata(ie,parameters,data)
                        #print(self.getexp(ie))
                        #print(ie)
                        ie = ie + 1
                        words[0]='bla'
                        parameters=dict()
                except KeyError:
                    self.logger.warning('fastforward') #no experiment found
        self.logger.info('{} experiments read'.format(ie-1))
        os.chdir(olddir)
