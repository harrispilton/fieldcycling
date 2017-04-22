import numpy as np
import scipy.signal

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
    """ fits exponential decay. 

    :param df: data
    :type df: pandas dataframe
    :param p0: initial set of parameters
    :type p0: list of 3 elements
    :param fit_option: 1) Linearize the system, and fit a line to the log of the data.
                            - would be prefered, but needs the y-axes offset.
                       2) Use a non-linear solver (e.g. scipy.optimize.curve_fit
    :type fit_option: int, 1 or 2
    """
    ''' 
    '''
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
