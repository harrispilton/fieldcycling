import numpy as np
import numpy.matlib
import scipy.signal
from num_string_eval import NumericStringParser

def fit_exp_linear(t, y, C=0):
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    popt = [A, K, C]
    return popt

def model_exp_dec(t, A, K, C):
    return A * np.exp(- K * t) + C
def fun_exp_dec(par,t,y):
    A, K, C = par
    return model_exp_dec(t, A,K,C) - y

def magnetization_fit(df,p0,fit_option):
    # fit exponential decay. 
    # Options:
    #    1) Linearize the system, and fit a line to the log of the data.
    #        - would be prefered, but needs the y-axes offset.
    #    2) Use a non-linear solver (e.g. scipy.optimize.curve_fit
    if fit_option ==1:
        C0 = 0 # offset
        popt = fit_exp_linear(df.tau, df.phi_normalized, C0)
    elif fit_option == 2:
        from scipy.optimize import leastsq
        popt, _ = leastsq(fun_exp_dec, p0  , args=(np.array(df.tau),np.array(df.phi_normalized)) )
    df['fit_phi'] = model_exp_dec(df.tau, *popt)
    df['phi_normalized'] = (df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
    return df, popt
    
def get_mag_amplitude(fid,startpoint, endpoint, nblk, bs):
    phi=np.zeros(nblk)
    for blk in range(nblk):
        start=startpoint + blk * bs-1
        end=endpoint + blk * bs
        phi[blk]=fid['magnitude'].iloc[start:end].sum() / (endpoint-startpoint)
    return phi

def get_x_axis(parameters):
    nsp = NumericStringParser()
    nblk = parameters['NBLK']
    T1MX = parameters['T1MX'] # T1MX is used in the 'eval' expressions below
    # TODO: not tested yet for all cases
    if parameters['BGRD'] == 'LIST':
        print('BGRD = LIST')
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
        print('BGRD = LIN')
        Tini = nsp.eval(parameters['BINI'].replace('T1MX', str(T1MX)))
        Tend = nsp.eval(parameters['BEND'].replace('T1MX', str(T1MX)))
        x = np.linspace(Tini, Tend, nblk) # re-create the time vector
    elif parameters['BGRD'] == 'LOG':
        print('BGRD = LOG')
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


##nextpow2 and freq_shift were copied from the internet. not carefully checked
def nextpow2(x):
    """Return the first integer N such that 2**N >= abs(x)"""

    return int(np.ceil(np.log2(np.abs(x))))


##the algorithm can produce artifacts!! use with care
def freq_shift(x, f_shift, dt):
    """
    Shift the specified signal by the specified frequency.
    """

    # Pad the signal with zeros to prevent the FFT invoked by the transform from
    # slowing down the computation:
    N_orig = len(x)
    N_padded = 2**nextpow2(N_orig)
    t = np.arange(0, N_padded)
    return (
        scipy.signal.hilbert(
            np.hstack(
                (x, np.zeros(N_padded-N_orig, x.dtype))
                )
            )*np.exp(2j*np.pi*f_shift*dt*t))[:N_orig].real
