import numpy as np
from numba import njit

def create_parameters(parameters, seed_signal=1, seed_network=1, error_points_t=3000, dt=1e-4, time_c=5):
    d_sim, d_pert, d_seed = {}, {}, {}
    if parameters == 'Biological':
        # Signal parameters
        d_sim['s_mag'] = 3 # signal magnitude
        d_sim['signal_noise'] = 0.5 # parameter that controls the maximum absolute amplitude of the *slow* noise for the signal

        # Network parameters
        d_sim['reset'] = 0.02 # Multiplicative term for actual resets Omega[i, i] *= 2 ** reset
        d_sim['threshold'] = 0.55
        d_sim['lambd'] = 100 # Time constant for the membrane Voltage
        d_sim['ref_t'] = 0.002 # refractory period in seconds

        #Simulation Parameters
        d_sim['noise'] = 0.5

    # elif parameters == 'Plain':
    #     # Signal parameters
    #     d_sim['s_mag'] = 3 # signal magnitude
    #     d_sim['signal_noise'] = 0 # parameter that controls the maximum absolute amplitude of the *slow* noise for the signal

    #     # Network parameters
    #     d_sim['reset'] = 0 # Multiplicative term for actual resets Omega[i, i] *= 2 ** reset
    #     d_sim['threshold'] = 0.5
    #     d_sim['lambd'] = 10 # Time constant for the membrane Voltage
    #     d_sim['ref_t'] = 0 # refractory period in seconds

    #     #Simulation Parameters
    #     d_sim['noise'] = 0

    # elif parameters == 'OnlySaturation':
    #     # Signal parameters
    #     d_sim['s_mag'] = 3 # signal magnitude
    #     d_sim['signal_noise'] = 0 # parameter that controls the maximum absolute amplitude of the *slow* noise for the signal

    #     # Network parameters
    #     d_sim['reset'] = 0 # Multiplicative term for actual resets Omega[i, i] *= 2 ** reset
    #     d_sim['threshold'] = 0.5
    #     d_sim['lambd'] = 10 # Time constant for the membrane Voltage
    #     d_sim['ref_t'] = 0.002 # refractory period in seconds

    #     #Simulation Parameters
    #     d_sim['noise'] = 0

    # elif parameters[:6] == 'Movies':
    #     # Signal parameters
    #     d_sim['s_mag'] = 3 # signal magnitude
    #     d_sim['signal_noise'] = 0 # parameter that controls the maximum absolute amplitude of the *slow* noise for the signal

    #     # Network parameters
    #     d_sim['reset'] = 0 # Multiplicative term for actual resets Omega[i, i] *= 2 ** reset
    #     d_sim['threshold'] = 0.5
    #     d_sim['lambd'] = 1 # Time constant for the membrane Voltage
    #     d_sim['ref_t'] = 0.002 # refractory period in seconds

    #     #Simulation Parameters
    #     d_sim['noise'] = 0
    #     if parameters == 'MoviesPaperBigThreshold':
    #         d_sim['ref_t'] = 0.0002 # refractory period in seconds
    #         d_sim['threshold'] = 0.7
    #         d_sim['lambd'] = 100 # Time constant for the membrane Voltage
    #     if parameters == 'MoviesBigThreshold':
    #         d_sim['threshold'] = 0.55
    #     if parameters == 'MoviesBigBIGThreshold':
    #         d_sim['threshold'] = 1
    # elif parameters == 'Animations':
    #     # Signal parameters
    #     d_sim['s_mag'] = 3 # signal magnitude
    #     d_sim['signal_noise'] = 0 # parameter that controls the maximum absolute amplitude of the *slow* noise for the signal

    #     # Network parameters
    #     d_sim['reset'] = 0 # Multiplicative term for actual resets Omega[i, i] *= 2 ** reset
    #     d_sim['threshold'] = 0.7
    #     d_sim['lambd'] = 100 # Time constant for the membrane Voltage
    #     d_sim['ref_t'] = 0.0002 # refractory period in seconds

    #     #Simulation Parameters
    #     d_sim['noise'] = 0

    d_seed['signal'] = seed_signal
    d_seed['network'] = seed_network
    # Simulation Parameters
    d_sim['time_r'] = 0.4
    d_sim['time_c'] = time_c
    d_sim['dt'] = dt
    d_sim['delay'] = 0
        

    d_sim['T_r'] = int(d_sim['time_r'] / d_sim['dt'])
    d_sim['T_c'] = int(d_sim['time_c'] / d_sim['dt'])
    d_sim['step_t'] = (d_sim['T_c'] + 1) // error_points_t
    d_sim['error_points_t'] = ((d_sim['T_c'] - 1) // d_sim['step_t']) + 1
    return d_sim, d_pert, d_seed


# def error_cube(dimension, N_aux = 100000, max_m = 100):
#     try:
#         res = np.load('aux/hypercube_distances_{}.npy'.format(N_aux))
#     except:
#         res = np.zeros(max_m)
#         for m in range(max_m):
#             print(m)
#             X = np.random.rand(N_aux, m + 1)
#             X -= 0.5 * np.ones(m + 1)
#             Y = np.linalg.norm(X, axis=1)
#             res[m] = np.mean(Y)
#         np.save('aux/hypercube_distances_{}'.format(N_aux), res)
#     out = res[dimension - 1]
#     return out

# def rotate_rnd(u, theta):
#     N = u.size
#     norm_u = np.linalg.norm(u)
#     if norm_u > 0:
#         ux = u / norm_u
#         v = np.random.randn(N)
#         v -= np.dot(v, ux) * ux
#         v /= np.linalg.norm(v)
#         w = np.cos(theta) * ux + np.sin(theta) * v
#         w /= np.linalg.norm(w)
#         return w
#     else:
#         return u

def ramping_signal(M, T, dt, s_mag=3, seed_signal=1):
    '''
    Inputs: 
        M: positive integer. Dimension of the signal being represented by the network.
        T: positive integer. Number of time steps sampled for the stimulus.
        dt: positive scalar. Magnitude of time step.
        s_mag: positive scalar. Standard Deviation of each signal component
        seed_signal: positive integer. Seed for defining what the actual point will be
    Outputs:
        X: (T, M) array. Signal to be represented.
        Xdot: (T, M) array. Time derivative of signal to be represented.
    '''
    np.random.seed(seed_signal) # Initialize the random seed
    p1 = s_mag * np.random.randn(M) # Draw the M-dimensional signal from a multivariate gaussian with mean 0 and standard deviation s_mag
    X = np.outer(np.linspace(0, 1, T + 1), p1) # Create the ramp array 
    Xdot = np.diff(X.T).T / dt # compute the exact derivative as the difference between time steps
    return X[:-1, :], Xdot

def constant_signal(x, T, dt, signal_noise=0, seed=1, frac=0.1):
    '''
    Inputs: 
        x: (M,) array. Signal value
        T: positive integer. Number of time steps sampled for the stimulus.
        dt: positive scalar. Magnitude of time step.
        seed: positive integer. Seed for defining what the actual noise will be
        signal_noise: positive scalar. Magnitude of slow noise
        frac: proportion. how much to attenuate the noise such that beggining and end match with x 
    Outputs:
        X: (T, M) array. Signal to be represented.
    '''
    M = x.shape[0]
    X = np.zeros((T, M))
    for i in range(M):
        X[:, i] = x[i] # populate the signal array
    if signal_noise != 0:
        np.random.seed(seed) # Initialize the random seed
        N = int(1 / dt) # perform running averages of size 1 second
        dX = np.random.randn(T + 2 * N, M) # create the initial white noise
        dXt = np.empty((T, M)) # initialize the signal noise vector that will be added to the signal
        plateau = np.ones(T) # this array will ensure that at both ends the signal noise vanishes
        no_ramp = int(T * frac)
        ramp = np.linspace(0, 1, no_ramp)
        plateau[:no_ramp] = ramp
        plateau[-no_ramp:] = ramp[::-1]
        for i in range(M):
            csdX = np.cumsum(dX[:, i])
            csdX = (csdX[N:] - csdX[:-N]) / float(N) # this step and the previous is an efficient way to perform a running average
            csdX = np.cumsum(csdX)
            dXt[:, i] = (csdX[N:] - csdX[:-N]) / float(N)
            dXt[:, i] *= plateau # here we make sure the signal error term vanishes at the beggining and end
            dXt[:, i] *= signal_noise / np.max(np.abs(dXt[:, i])) # here we make sure that the maximum deviation of the signal error term is signal_noise
        X += dXt # Now we add the noise term to the signal
    return X

# def signal(M, T_r, T_c, dt, s_mag=10, seed_signal=1, signal_noise=0):
#     X_r, Xdot_r = ramping_signal(M, T_r, dt, s_mag=s_mag, seed_signal=seed_signal)
#     X_c = constant_signal(X_r[-1, :], T_c, dt, signal_noise=signal_noise, seed=seed_signal)
#     X = np.zeros((T_r + T_c, M))
#     X[:T_r, :] = X_r
#     X[T_r:, :] = X_c
#     Xdot = np.zeros((T_r + T_c, M))
#     Xdot[:T_r, :] = Xdot_r
#     return X, Xdot

def sim(F, Omega, Th, X, dt, dV, Lambd, lambd, delay, R0, V0, ref_t=0.003, jump=1, voltages=False, thresholds=False, tmpspikefailure=0, synapticnoise=0, seed=1):
    '''
    Inputs: 
        F: (N, M) array. Encoding Weights. N is the number of neurons and M is the signal dimension
        Omega: (N, N) array. Recurrent Connections. N is the number of neurons
        Th: (N,) array. Neuron Thresholds. N is the number of neurons
        X: (T, M) array. Signal to be represented. T is the number of time steps and M is the signal dimension
        dt: scalar. Time step of the simulation
        dV: (T, N) array. Additive term on the voltages traces (should already be multiplied by the time step)
        Lambd: (N,) array of positive scalars. Neuronal time constants
        lambd: positive scalar. Decoder time constant
        delay: non-negative integer. Number of time steps of the network delay
        R0: (N,) array. Initial Rate Values
        V0: (N,) array. Initial Voltage values
        ref_t: non-negative scalar. refractory period duration in seconds
        jump: positive integer. time step jumps to store the Rates
        voltages: boolean. True, saves voltage history
        thresholds: boolean. True, Thresholds are variable
        tmpspikefailure: probability os a synapse fail at each spike.
    Outputs:
        R: ((T + jump - 1) // jump, N) array. Firing Rates approximation of each neuron. T is the number of time steps and N is the number of neurons
        Spikes: list of lists. Each element in Spikes is a list containing the simulation time where that neuron fired
        Rn: (N,) array. Final Rates of the simulation
        Vn: (N,) array. Final Voltages of the simulation
    '''
    np.random.seed(seed)
    ref = max(int(np.ceil(ref_t / dt)), 1) # minimum number of time steps between same-neuron spikes
    dtlamb = 1 - dt * lambd
    dtLamb = 1 - dt * Lambd
    N = F.shape[0]  # N is the number of Neurons
    T = X.shape[0]  # T is the total time iterations
    Vc = V0.copy() # Array with the current Voltages
    Rc = R0.copy() # Array with the current Rates
    if voltages:
        V = np.zeros(((T + jump - 1) // jump, N)) # Array with the Voltages
    R = np.zeros(((T + jump - 1) // jump, N))  # Array with the approximation of the firing rates
    Spikes = [[] for i in range(N)]  # list of lenght Neurons of list of lenght no_spikes
    FX = np.dot(X, F.T)  # Network inputs
    Ref = -1 * np.ones(N)
    Th_c = Th
    saturation = np.zeros(N)
    if delay > 0:
        Omega_resets = np.diag(np.diag(Omega))
        Omega_delay = Omega - Omega_resets
        # Array with spikes for the last delay time steps
        S = np.zeros((delay, N))
        for i in range(T):
            if thresholds:
                Th_c = Th[i, :]
            if synapticnoise == 0:
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            else:
                Omega_delay_mod = Omega_delay * (1 - synapticnoise) ** (2 * np.random.rand(N, N) - 1)
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay_mod, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay_mod * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            S[(i % delay), :] = 0
            SN = np.where((Vc > Th_c) & (Ref <= 0))[0] # SN is a list containing neurons that should fire
            for j in range(SN.size):
                Rc[SN[j]] += 1
                Spikes[SN[j]].append(i)
                S[(i % delay), SN[j]] = 1
                Ref[SN[j]] = ref
                # Adding the resets of each spike through Omega_resets
                Vc -= Omega_resets[:, SN[j]]
            Rc *= dtlamb            
            Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
            Ref -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
                if voltages:
                    V[i // jump, :] = Vc
    else:
        for i in range(T):
            if thresholds:
                Th_c = Th[i, :]
            # SN is a list containing neurons that should fire
            VR = ((Vc > Th_c) & (Ref <= 0))
            while np.any(VR):
                SN = np.where(VR)[0]
                SNN = np.where(Vc[SN] == np.max(Vc[SN]))[0]
                sn = SN[np.random.choice(SNN)]
                Rc[sn] += 1
                Spikes[sn].append(i)
                Ref[sn] = ref
                # Adding the resets of each spike through Omega
                if synapticnoise == 0:
                    if tmpspikefailure == 0:
                        Vc -= Omega[:, sn]
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega[:, sn] * (aux >= tmpspikefailure)
                else:
                    scaling_Omega = (1 - synapticnoise) ** (2 * np.random.rand(N) - 1)
                    scaling_Omega[sn] = 1
                    if tmpspikefailure == 0:
                        Vc -= Omega[:, sn] * scaling_Omega
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega[:, sn] * (aux >= tmpspikefailure) * scaling_Omega
                VR = ((Vc > Th_c) & (Ref <= 0))
            if np.any(Vc > Th_c):
                for j in range(N):
                    if Vc[j] > Th_c[j]:
                        saturation[j] += 1
            Rc *= dtlamb
            Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
            Ref -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
                if voltages:
                    V[i // jump, :] = Vc
    print('dt', dt, 'ref', ref, 'delay', delay, 'N', N, 'most saturated neuron: ', np.max(saturation) / T, 'average saturated neuron:', np.mean(saturation) / T)
    if voltages:
        return R, Spikes, V, Rc, Vc
    else:
        return R, Spikes, Rc, Vc


def sim_slow(A, B, F, Omega_fast, Omega_slow, Th, dt, dV, lambd, delay, R0, V0, ref_t=0.003, jump=1, voltages=False, tmpspikefailure=0, synapticnoise=0, thresholds=False, readout_correction=True, D=None, seed=1):
    '''
    Inputs: 
        A: (M, M) array. Dynamical system matrix A. M is the signal dimension
        B: (T, M) array. Dynamical system vector b. T is the number of time steps and M is the signal dimension
        F: (N, M) array. Encoding Weights. N is the number of neurons and M is the signal dimension
        Omega_fast: (N, N) array. Fast Recurrent Connections. N is the number of neurons
        Omega_slow: (N, N) array. Slow Recurrent Connections. N is the number of neurons
        Th: (N,) array. Neuron Thresholds. N is the number of neurons
        dt: scalar. Time step of the simulation
        dV: (T, N) array. Additive term on the voltages traces (should already be multiplied by the time step)
        lambd: positive scalar. Decoder time constant
        delay: non-negative integer. Number of time steps of the network delay
        R0: (N,) array. Initial Rate Values
        V0: (N,) array. Initial Voltage values
        ref_t: non-negative scalar. refractory period duration in seconds
        jump: positive integer. time step jumps to store the Rates
        voltages: boolean. True, saves voltage history
        readout_correction: boolean. True, re-scales readout. If True, the decoder should be provided
    Outputs:
        R: ((T + jump - 1) // jump, N) array. Firing Rates approximation of each neuron. T is the number of time steps and N is the number of neurons
        Spikes: list of lists. Each element in Spikes is a list containing the simulation time where that neuron fired
        Rn: (N,) array. Final Rates of the simulation
        Vn: (N,) array. Final Voltages of the simulation
    '''
    np.random.seed(seed)
    ref = max(int(np.ceil(ref_t / dt)), 1) # minimum number of time steps between same-neuron spikes
    dtlamb = 1 - dt * lambd
    N = F.shape[0]  # N is the number of Neurons
    T = B.shape[0]  # T is the total time iterations
    Vc = V0.copy() # Array with the current Voltages
    Rc = R0.copy() # Array with the current Rates
    if voltages:
        V = np.zeros(((T + jump - 1) // jump, N)) # Array with the Voltages
    R = np.zeros(((T + jump - 1) // jump, N))  # Array with the approximation of the firing rates
    Spikes = [[] for i in range(N)]  # list of lenght Neurons of list of lenght no_spikes
    FB = np.dot(B, F.T)  # Network inputs
    Ref = -1 * np.ones(N)
    Th_c = Th
    saturation = np.zeros(N)
    if delay > 0:
        Omega_resets = np.diag(np.diag(Omega_fast))
        Omega_delay = Omega_fast - Omega_resets
        # Array with spikes for the last delay time steps
        S = np.zeros((delay, N))
        for i in range(T):
            if thresholds:
                Th_c = Th[i, :]
            if synapticnoise == 0:
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            else:
                Omega_delay_mod = Omega_delay * (1 - synapticnoise) ** (2 * np.random.rand(N, N) - 1)
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay_mod, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay_mod * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            S[(i % delay), :] = 0
            SN = np.where((Vc > Th_c) & (Ref <= 0))[0] # SN is a list containing neurons that should fire
            for j in range(SN.size):
                Rc[SN[j]] += 1
                Spikes[SN[j]].append(i)
                S[(i % delay), SN[j]] = 1
                Ref[SN[j]] = ref
                # Adding the resets of each spike through Omega_resets
                Vc -= Omega_resets[:, SN[j]]
            Rc *= dtlamb
            if readout_correction:
                norm_xh = max(2, np.linalg.norm(np.dot(D, Rc)))
                Vc = dtlamb * Vc + dt * FB[i, :] + dt * np.dot(Omega_slow, Rc) * (1 + (Th[0] - 0.5) / norm_xh) + dV[i, :]
            else:
                Vc = dtlamb * Vc + dt * FB[i, :] + dt * np.dot(Omega_slow, Rc) + dV[i, :]
            Ref -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
                if voltages:
                    V[i // jump, :] = Vc
    else:
        for i in range(T):
            if thresholds:
                Th_c = Th[i, :]
            VR = ((Vc > Th_c) & (Ref <= 0))
            while np.any(VR):
                SN = np.where(VR)[0]
                SNN = np.where(Vc[SN] == np.max(Vc[SN]))[0]
                sn = SN[np.random.choice(SNN)]
                Rc[sn] += 1
                Spikes[sn].append(i)
                Ref[sn] = ref
                # Adding the resets of each spike through Omega
                if synapticnoise == 0:
                    if tmpspikefailure == 0:
                        Vc -= Omega_fast[:, sn]
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega_fast[:, sn] * (aux >= tmpspikefailure)
                else:
                    scaling_Omega = (1 - synapticnoise) ** (2 * np.random.rand(N) - 1)
                    scaling_Omega[sn] = 1
                    if tmpspikefailure == 0:
                        Vc -= Omega_fast[:, sn] * scaling_Omega
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega_fast[:, sn] * (aux >= tmpspikefailure) * scaling_Omega
                VR = ((Vc > Th_c) & (Ref <= 0))
            if np.any(Vc > Th_c):
                for j in range(N):
                    if Vc[j] > Th_c[j]:
                        saturation[j] += 1
            Rc *= dtlamb
            if readout_correction:
                norm_xh = max(2, np.linalg.norm(np.dot(D, Rc)))
                Vc = dtlamb * Vc + dt * FB[i, :] + dt * np.dot(Omega_slow, Rc) * (1 + (Th[0] - 0.5) / norm_xh) + dV[i, :]
            else:
                Vc = dtlamb * Vc + dt * FB[i, :] + dt * np.dot(Omega_slow, Rc) + dV[i, :]
            Ref -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
                if voltages:
                    V[i // jump, :] = Vc
    print('dt', dt, 'ref', ref, 'delay', delay, 'N', N, 'most saturated neuron: ', np.max(saturation) / T, 'average saturated neuron:', np.mean(saturation) / T)
    if voltages:
        return R, Spikes, V, Rc, Vc
    else:
        return R, Spikes, Rc, Vc



# def sim_movies(F, Omega, Th, X, dt, dV, Lambd, lambd, delay, R0, V0, ref_t=0.003, jump=1, voltages=False, thresholds=False, tmpspikefailure=0, synapticnoise=0, seed=1):
#     '''
#     Inputs: 
#         F: (N, M) array. Encoding Weights. N is the number of neurons and M is the signal dimension
#         Omega: (N, N) array. Recurrent Connections. N is the number of neurons
#         Th: (N,) array. Neuron Thresholds. N is the number of neurons
#         X: (T, M) array. Signal to be represented. T is the number of time steps and M is the signal dimension
#         dt: scalar. Time step of the simulation
#         dV: (T, N) array. Additive term on the voltages traces (should already be multiplied by the time step)
#         Lambd: (N,) array of positive scalars. Neuronal time constants
#         lambd: positive scalar. Decoder time constant
#         delay: non-negative integer. Number of time steps of the network delay
#         R0: (N,) array. Initial Rate Values
#         V0: (N,) array. Initial Voltage values
#         ref_t: non-negative scalar. refractory period duration in seconds
#         jump: positive integer. time step jumps to store the Rates
#         voltages: boolean. True, saves voltage history
#         thresholds: boolean. True, Thresholds are variable
#         tmpspikefailure: probability os a synapse fail at each spike.
#     Outputs:
#         R: ((T + jump - 1) // jump, N) array. Firing Rates approximation of each neuron. T is the number of time steps and N is the number of neurons
#         Spikes: list of lists. Each element in Spikes is a list containing the simulation time where that neuron fired
#         Rn: (N,) array. Final Rates of the simulation
#         Vn: (N,) array. Final Voltages of the simulation
#     '''
#     np.random.seed(seed)
#     ref = max(int(np.ceil(ref_t / dt)), 1) # minimum number of time steps between same-neuron spikes
#     dtlamb = 1 - dt * lambd
#     dtLamb = 1 - dt * Lambd
#     N = F.shape[0]  # N is the number of Neurons
#     T = X.shape[0]  # T is the total time iterations
#     Vc = V0.copy() # Array with the current Voltages
#     Rc = R0.copy() # Array with the current Rates
#     if voltages:
#         V = np.zeros(((T + jump - 1) // jump, N)) # Array with the Voltages
#     R = np.zeros(((T + jump - 1) // jump, N))  # Array with the approximation of the firing rates
#     Spikes = [[] for i in range(N)]  # list of lenght Neurons of list of lenght no_spikes
#     FX = np.dot(X, F.T)  # Network inputs
#     Ref = -1 * np.ones(N)
#     Th_c = Th
#     if delay > 0:
#         Omega_resets = np.diag(np.diag(Omega))
#         Omega_delay = Omega - Omega_resets
#         # Array with spikes for the last delay time steps
#         S = np.zeros((delay, N))
#         for i in range(T):
#             if thresholds:
#                 Th_c = Th[i, :]
#             if synapticnoise == 0:
#                 if tmpspikefailure == 0:
#                     Vc -= np.dot(Omega_delay, S[(i % delay), :])
#                 else:
#                     Vc -= np.dot(Omega_delay * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
#             else:
#                 Omega_delay_mod = Omega_delay * (1 - synapticnoise) ** (2 * np.random.rand(N, N) - 1)
#                 if tmpspikefailure == 0:
#                     Vc -= np.dot(Omega_delay_mod, S[(i % delay), :])
#                 else:
#                     Vc -= np.dot(Omega_delay_mod * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
#             S[(i % delay), :] = 0
#             SN = np.where((Vc > Th_c) & (Ref <= 0))[0] # SN is a list containing neurons that should fire
#             for j in range(SN.size):
#                 Rc[SN[j]] += 1
#                 Spikes[SN[j]].append(i)
#                 S[(i % delay), SN[j]] = 1
#                 Ref[SN[j]] = ref
#                 # Adding the resets of each spike through Omega_resets
#                 Vc -= Omega_resets[:, SN[j]]
#             Rc *= dtlamb
#             Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
#             Ref -= 1
#             if i % jump == 0:
#                 R[i // jump, :] = Rc
#                 if voltages:
#                     V[i // jump, :] = Vc
#     else:
#         for i in range(T):
#             if thresholds:
#                 Th_c = Th[i, :]
#             # SN is a list containing neurons that should fire
#             VR = ((Vc > Th_c) & (Ref <= 0))
#             if np.any(VR):
#                 SN = np.where(VR)[0]
#                 SNN = np.where(Vc[SN] == np.max(Vc[SN]))[0]
#                 sn = SN[np.random.choice(SNN)]
#                 Rc[sn] += 1
#                 Spikes[sn].append(i)
#                 Ref[sn] = ref
#                 # Adding the resets of each spike through Omega
#                 if synapticnoise == 0:
#                     if tmpspikefailure == 0:
#                         Vc -= Omega[:, sn]
#                     else:
#                         aux = np.random.rand(N)
#                         aux[sn] = 1
#                         Vc -= Omega[:, sn] * (aux >= tmpspikefailure)
#                 else:
#                     scaling_Omega = (1 - synapticnoise) ** (2 * np.random.rand(N) - 1)
#                     scaling_Omega[sn] = 1
#                     if tmpspikefailure == 0:
#                         Vc -= Omega[:, sn] * scaling_Omega
#                     else:
#                         aux = np.random.rand(N)
#                         aux[sn] = 1
#                         Vc -= Omega[:, sn] * (aux >= tmpspikefailure) * scaling_Omega
#             Rc *= dtlamb
#             Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
#             Ref -= 1
#             if i % jump == 0:
#                 R[i // jump, :] = Rc
#                 if voltages:
#                     V[i // jump, :] = Vc
#     print('dt', dt, 'ref', ref, 'delay', delay, 'N', N)
#     if voltages:
#         return R, Spikes, V, Rc, Vc
#     else:
#         return R, Spikes, Rc, Vc




@njit
def simnumba(F, Omega, Th, X, dt, dV, Lambd, lambd, delay, R0, V0, ref_t=0.003, jump=1, voltages=False, tmpspikefailure=0, synapticnoise=0, seed=1):
    '''
    Inputs: 
        F: (N, M) array. Encoding Weights. N is the number of neurons and M is the signal dimension
        Omega: (N, N) array. Recurrent Connections. N is the number of neurons
        Th: (N,) array. Neuron Thresholds. N is the number of neurons
        X: (T, M) array. Signal to be represented. T is the number of time steps and M is the signal dimension
        dt: scalar. Time step of the simulation
        dV: (T, N) array. Additive term on the voltages traces (should already be multiplied by the time step)
        Lambd: (N,) array of positive scalars. Neuronal time constants
        lambd: positive scalar. Decoder time constant
        delay: non-negative integer. Number of time steps of the network delay
        R0: (N,) array. Initial Rate Values
        V0: (N,) array. Initial Voltage values
        ref_t: non-negative scalar. refractory period duration in seconds
        jump: positive integer. time step jumps to store the Rates
        voltages: boolean. True, saves voltage history
        tmpspikefailure: probability os a synapse fail at each spike.
    Outputs:
        R: ((T + jump - 1) // jump, N) array. Firing Rates approximation of each neuron. T is the number of time steps and N is the number of neurons
        Spikes: list of lists. Each element in Spikes is a list containing the simulation time where that neuron fired
        Rn: (N,) array. Final Rates of the simulation
        Vn: (N,) array. Final Voltages of the simulation
    '''
    np.random.seed(seed)
    ref = max(int(np.ceil(ref_t / dt)), 1) # minimum number of time steps between same-neuron spikes
    dtlamb = 1 - dt * lambd
    dtLamb = 1 - dt * Lambd
    N = F.shape[0]  # N is the number of Neurons
    T = X.shape[0]  # T is the total time iterations
    Vc = V0.copy() # Array with the current Voltages
    Rc = R0.copy() # Array with the current Rates
    R = np.zeros(((T + jump - 1) // jump, N))  # Array with the approximation of the firing rates
    Spikes = []
    for i in range(N):
        Spikes.append([np.int64(x) for x in range(0)])
    FX = np.dot(X, F.T)  # Network inputs
    Ref = np.array([-1 for i in range(N)])
    Th_c = Th
    saturation = np.zeros(N)
    if delay > 0:
        Omega_resets = np.diag(np.diag(Omega))
        Omega_delay = Omega - Omega_resets
        # Array with spikes for the last delay time steps
        S = np.zeros((delay, N))
        for i in range(T):
            if synapticnoise == 0:
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            else:
                Omega_delay_mod = Omega_delay * (1 - synapticnoise) ** (2 * np.random.rand(N, N) - 1)
                if tmpspikefailure == 0:
                    Vc -= np.dot(Omega_delay_mod, S[(i % delay), :])
                else:
                    Vc -= np.dot(Omega_delay_mod * (np.random.rand(N, N) >= tmpspikefailure), S[(i % delay), :])
            S[(i % delay), :] = 0
            SN = np.where((Vc > Th_c) & (Ref <= 0))[0] # SN is a list containing neurons that should fire
            for j in range(SN.size):
                Rc[SN[j]] += 1
                Spikes[SN[j]].append(i)
                S[(i % delay), SN[j]] = 1
                Ref[SN[j]] = ref
                # Adding the resets of each spike through Omega_resets
                Vc -= Omega_resets[:, SN[j]]
            for i_neuron in range(N):
                Rc[i_neuron] *= dtlamb
            Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
            for i_neuron in range(N):
                Ref[i_neuron] -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
    else:
        for i in range(T):
            VR = ((Vc > Th_c) & (Ref <= 0))
            while np.any(VR):
                SN = np.where(VR)[0]
                SNN = np.where(Vc[SN] == np.max(Vc[SN]))[0]
                sn = SN[np.random.choice(SNN)]
                Rc[sn] += 1
                Spikes[sn].append(i)
                Ref[sn] = ref
                # Adding the resets of each spike through Omega
                if synapticnoise == 0:
                    if tmpspikefailure == 0:
                        Vc -= Omega[:, sn]
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega[:, sn] * (aux >= tmpspikefailure)
                else:
                    scaling_Omega = (1 - synapticnoise) ** (2 * np.random.rand(N) - 1)
                    scaling_Omega[sn] = 1
                    if tmpspikefailure == 0:
                        Vc -= Omega[:, sn] * scaling_Omega
                    else:
                        aux = np.random.rand(N)
                        aux[sn] = 1
                        Vc -= Omega[:, sn] * (aux >= tmpspikefailure) * scaling_Omega
                VR = ((Vc > Th_c) & (Ref <= 0))
            if np.any(Vc > Th_c):
                for j in range(N):
                    if Vc[j] > Th_c[j]:
                        saturation[j] += 1
            for i_neuron in range(N):
                Rc[i_neuron] *= dtlamb
            Vc = dtLamb * Vc + dt * FX[i, :] + dV[i, :]
            for i_neuron in range(N):
                Ref[i_neuron] -= 1
            if i % jump == 0:
                R[i // jump, :] = Rc
    print('dt', dt, 'ref', ref, 'delay', delay, 'N', N, 'most saturated neuron: ', np.max(saturation) / T, 'average saturated neuron:', np.mean(saturation) / T)
    return R, Spikes, Rc, Vc


def create_param(M, N, threshold=0.5, reset=0, seed_network=1, hypercube=False, regular=False):
    np.random.seed(seed_network)
    D = np.random.randn(M, N)
    if regular and M == 2:
        angles = np.arange(N) * 2 * np.pi / N
        if N == 6:
            angles = np.array([0, np.pi / 2, -12 * np.pi / 10, -9 * np.pi / 10, -6 * np.pi / 10, -3 * np.pi / 10])
        D[0, :] = np.sin(angles)
        D[1, :] = np.cos(angles)
    if hypercube:
        Q = np.random.randn(M, M)
        Q, R = np.linalg.qr(Q)
        if N >= M:
            D[:, :M] = Q
        else:
            D[:, :N] = Q[:, :N]
        if N >= 2 * M:
            D[:, M:(2 * M)] = -Q
        for i in range(2 * M, N):
            D[:, i] /= np.linalg.norm(D[:, i])
    else:
        for i in range(N):
            D[:, i] /= np.linalg.norm(D[:, i])
    F = D.copy().T
    Omega = np.dot(D.T, D)
    for i in range(N):
        Omega[i, i] *= 2 ** reset
    Th = threshold * np.ones(N)
    return F, Omega, Th, D

def modify_params(T, F, Omega, Th, D, d, dt, seed=1, random_death = False):
    '''
    Inputs:
        d: dictionary
            d['deaths'] = proportion of neurons that are eliminated
            d['births'] = proportion of neurons that are added
            d['noise'] = standard deviation of noise
            d['th-one-random'] = amount of current injected on a single random neuron
            d['th-one-star'] = amount of current injected on the most active neuron that is still alive 
            d['enc-dec-mismatch'] = angle of mismatch. Note that this overwrites any previous recurrent matrix definitions
            d['weak-sparse'] = eliminating this value% of the weakest synapses of the alive neurons
            d['rec-noise'] = every synapse (except resets) is multiplied by a factor up to (1 - value) or its inverse
            and many more
    '''
    np.random.seed(seed)
    M, N = D.shape
    previously_alive_neurons = np.where(Th < np.inf)[0]
    N_previously_alive_neurons = previously_alive_neurons.shape[0]
    previously_dead_neurons = np.where(Th == np.inf)[0]
    N_previously_dead_neurons = previously_dead_neurons.shape[0]

    # Kill neurons, add neurons and add voltage noise
    if random_death:
        N_kill = int(np.round(N_previously_alive_neurons * np.around(d['deaths'], 2)))
        neuron_order = np.random.permutation(previously_alive_neurons)
        killed_neurons = neuron_order[:N_kill]
        dead_neurons = np.concatenate((killed_neurons, previously_dead_neurons))
        alive_neurons = neuron_order[N_kill:]
    else:
        dead_neurons = np.array(d['deaths']).astype(int)
        alive_neurons = np.array([i for i in range(N) if i not in d['deaths']])
    N_dead_neurons = dead_neurons.shape[0]
    N_alive = alive_neurons.shape[0]
    N_birth = int(np.round(N_alive * np.around(d['births'], 2)))
    N_t = N + N_birth
    alive_neurons = np.concatenate((alive_neurons, np.arange(N, N_t)))
    N_alive = alive_neurons.shape[0]
    D_m = np.empty((M, N_t))
    D_m[:, :N] = D
    D_m[:, N:] = np.random.randn(M, N_birth)
    for i in range(N, N_t):
        D_m[:, i] /= np.linalg.norm(D_m[:, i])
    F_m = np.empty((N_t, M))
    F_m[:N, :] = F
    F_m[N:, :] = D_m[:, N:].T
    F_m[dead_neurons, :] = 0
    Omega_m = np.empty((N_t, N_t))
    Omega_m[:N, :N] = Omega
    for i_reset in d['reset']:
        amount_reset = 0.4
        Omega_m[i_reset, i_reset] *= 2 ** amount_reset
    Omega_m[N:, :] = np.dot(F_m[N:, :], D_m)
    Omega_m[:N, N:] = np.dot(F_m[:N, :], D_m[:, N:])
    Omega_m[dead_neurons, :] = 0
    Omega_m[:, dead_neurons] = 0
    Th_m = np.empty(N_t)
    Th_m[:N] = Th * d['threshold']    
    Th_m[N:N_t] = np.mean(Th[previously_alive_neurons]) * np.ones(N_birth)
    Th_m[dead_neurons] = np.inf
    # if d['inhi-random'] != 0:
    #     N_opto = int(np.round(N_alive * np.around(d['inhi-random'], 2)))
    #     new_neuron_order = np.random.permutation(alive_neurons)
    #     neurons_opto = new_neuron_order[:N_opto]
    #     Th_m[neurons_opto] += 0.3
    # if d['exci-random'] != 0:
    #     N_opto = int(np.round(N_alive * np.around(d['exci-random'], 2)))
    #     new_neuron_order = np.random.permutation(alive_neurons)
    #     neurons_opto = new_neuron_order[:N_opto]
    #     Th_m[neurons_opto] -= 0.3
    # if d['inhi-active'] != 0:
    #     N_active = d['active_neurons_sorted'].shape[0]
    #     N_opto = int(np.round(N_active * np.around(d['inhi-active'], 2)))
    #     neurons_opto = d['active_neurons_sorted'][-N_opto:]
    #     Th_m[neurons_opto] += 0.3
    # if d['exci-active'] != 0:
    #     N_active = d['active_neurons_sorted'].shape[0]
    #     N_opto = int(np.round(N_active * np.around(d['exci-active'], 2)))
    #     neurons_opto = d['active_neurons_sorted'][-N_opto:]
    #     Th_m[neurons_opto] -= 0.3
    # if d['inhi-alignd'] != 0:
    #     N_alignd = d['aligned_neurons_sorted'].shape[0]
    #     N_opto = int(np.round(N_alignd * np.around(d['inhi-alignd'], 2)))
    #     neurons_opto = d['aligned_neurons_sorted'][-N_opto:]
    #     Th_m[neurons_opto] += 0.3
    # if d['exci-alignd'] != 0:
    #     N_alignd = d['aligned_neurons_sorted'].shape[0]
    #     N_opto = int(np.round(N_alignd * np.around(d['exci-alignd'], 2)))
    #     neurons_opto = d['aligned_neurons_sorted'][-N_opto:]
    #     Th_m[neurons_opto] -= 0.3

    dV = np.sqrt(dt) * d['noise'] * np.random.randn(T, N_t)
    p_opto = 170
    dV[:, d['excitation']] += p_opto * dt
    dV[:, d['inhibition']] -= p_opto * dt
    dV[:, dead_neurons] = 0

    # Change Connections
    Omega_m_alive = Omega_m[np.ix_(alive_neurons, alive_neurons)]
    if d['enc-dec-mismatch'] != 0:
        for i in alive_neurons:
            F_m[i, :] = rotate_rnd(D_m[:, i], d['enc-dec-mismatch'])
        Omega_m = np.dot(F_m, D_m)
    if d['weak-sparse'] != 0:
        Omega_m_alive[np.abs(Omega_m_alive) <= np.quantile(np.abs(Omega_m_alive), d['weak-sparse'])] = 0
        Omega_m[np.ix_(alive_neurons, alive_neurons)] = Omega_m_alive
    if d['rec-noise'] != 0:
        aux = np.diag(Omega_m_alive)
        Omega_m_alive *= (1 - d['rec-noise']) ** (2 * np.random.rand(N_alive, N_alive) - 1)
        Omega_m_alive[np.arange(N_alive), np.arange(N_alive)] = aux
        Omega_m[np.ix_(alive_neurons, alive_neurons)] = Omega_m_alive
    if d['remove-excitation']:
        Omega_m = np.maximum(Omega_m, 0)
    if d['reduced-exci'] != 0:
        exci_t = np.percentile(Omega_m[Omega_m < 0], d['reduced-exci'] * 100)
        Omega_m[Omega_m <= exci_t] = 0

    return F_m, Omega_m, Th_m, D_m, dV, N_t, N_dead_neurons

def build_dict(noise=0, random_death=False):
    d = {}
    d['baseline'] = 0
    if random_death:
        d['deaths'] = 0
    else:
        d['deaths'] = []
    d['births'] = 0
    d['noise'] = noise
    d['reset'] = []
    d['enc-dec-mismatch'] = 0
    d['weak-sparse'] = 0
    d['rec-noise'] = 0
    d['tmp-spike-failure'] = 0
    d['synaptic-noise'] = 0
    d['remove-excitation'] = False
    d['delay'] = 0
    # d['inhi-random'] = 0
    # d['inhi-active'] = 0
    # d['inhi-alignd'] = 0
    # d['exci-random'] = 0
    # d['exci-active'] = 0
    # d['exci-alignd'] = 0
    d['active_neurons_sorted'] = []
    d['aligned_neurons_sorted'] = []
    d['excitation'] = []
    d['inhibition'] = []
    d['reduced-exci'] = 0
    d['threshold'] = 1
    d['tau'] = 0
    return d


# def balance_trials(Dimensions, Expansions, starting_seedS, no_seedsS, min_seedsS):
#     no_dim = len(Dimensions)
#     no_exp = len(Expansions)
#     decay_exp = 0.3
#     decay_dim = 1
#     starting_seed = (starting_seedS * np.ones((no_dim, no_exp))).astype(int)
#     no_seeds = (no_seedsS * np.ones((no_dim, no_exp))).astype(int)
#     for i, s in enumerate(Dimensions):
#         for j, e in enumerate(Expansions):
#             no_seeds[i, j] = int(max(min_seedsS, np.exp(-decay_dim * s / 100) * np.exp(-decay_exp * e / 50) * no_seedsS))
#     return no_seeds, starting_seed


def ramp(dim, cov, d_sim, d_seed, recycle_runs=True, more_info=False, numba=True):
    s_mag, threshold, reset = d_sim['s_mag'], d_sim['threshold'], d_sim['reset']
    lambd, ref_t, delay = d_sim['lambd'], d_sim['ref_t'], d_sim['delay']
    noise, dt, time_r, T_r = d_sim['noise'], d_sim['dt'], d_sim['time_r'], d_sim['T_r']
    seed_signal, seed_network = d_seed['signal'], d_seed['network']

    str_ramp = str((dim, cov, s_mag, threshold, reset, lambd, ref_t, delay, noise, dt, time_r, numba, seed_signal, seed_network))
    str_ramp = clean_string(str_ramp)

    N = int(np.around(dim * cov))
    run_r = True
    if recycle_runs:
        try:
            if more_info:
                res_r = np.load('aux/start/' + str_ramp + '_more_info.npz')
                X_r = res_k['X_r']
                Xh_nc = res_k['Xh_nc']
                Xh_c = res_k['Xh_c']
                Spikes_r = res_k['Spikes_r']
            else:
                res_r = np.load('aux/start/' + str_ramp + '.npz')
            R0 = res_r['R0']
            V0 = res_r['V0']
            F = res_r['F']
            Omega = res_r['Omega']
            Th = res_r['Th']
            D = res_r['D']
            x = res_r['x']
            run_r = False
            Lambd = lambd * np.ones(N)
        except:
            pass
    if run_r:
        # Either we didn't have or we are forcing it to recompute (the initial ramping activity)
        print(str_ramp)
        X_r, C_r = ramping_signal(dim, T_r, dt, s_mag=s_mag, seed_signal=seed_signal)
        x = X_r[-1, :]
        C_r += lambd * X_r
        F, Omega, Th, D = create_param(dim, N, threshold=threshold, reset=reset, seed_network=seed_network, hypercube=False)
        dV = np.sqrt(dt) * noise * np.random.randn(T_r, N)
        delay_dt = int(np.ceil((delay / 1000) / dt))
        Lambd = lambd * np.ones(N)
        if numba:
            R_r, Spikes_r, R0, V0 = simnumba(F, Omega, Th, C_r, dt, dV, Lambd, lambd, delay_dt, np.zeros(N), np.zeros(N), ref_t=ref_t, jump=1, seed = seed_signal * seed_network)
        else:
            R_r, Spikes_r, R0, V0 = sim(F, Omega, Th, C_r, dt, dV, Lambd, lambd, delay_dt, np.zeros(N), np.zeros(N), ref_t=ref_t, jump=1, seed = seed_signal * seed_network)
        if more_info:
            Xh_nc = np.dot(R_r, D.T)
            Xh_c = readout_correction(R_r, D, threshold)
            np.savez('aux/start/' + str_ramp + '_more_info', R0=R0, V0=V0, F=F, Omega=Omega, Th=Th, D=D, x=x, X_r=X_r, Xh_nc=Xh_nc, Xh_c=Xh_c, Spikes_r=Spikes_r) # Save this variables for future usage
        else:
            np.savez('aux/start/' + str_ramp, R0=R0, V0=V0, F=F, Omega=Omega, Th=Th, D=D, x=x) # Save this variables for future usage
    if more_info:
        return R0, V0, F, Omega, Th, Lambd, D, x, str_ramp, X_r, Xh_nc, Xh_c, Spikes_r
    else:
        return R0, V0, F, Omega, Th, Lambd, D, x, str_ramp


def keep(d, R0, V0, F, Omega, Th, Lambd, D, x, str_ramp, d_sim, d_seed, manipulation='baseline', value=0, recycle_runs=True, more_info=False, no_sig=2, numba=True):
    lambd, ref_t, threshold = d_sim['lambd'], d_sim['ref_t'], d_sim['threshold']
    dt, time_c, T_c = d_sim['dt'], d_sim['time_c'], d_sim['T_c']
    step_t, signal_noise, seed = d_sim['step_t'], d_sim['signal_noise'], d_seed['trial']
    if more_info:
        step_t = 1
    str_trial = str(('{value:.{no_sig}f}'.format(value=value, no_sig=no_sig), signal_noise, time_c, step_t, numba, str_ramp, seed))
    str_trial = clean_string(str_trial)
    run_k = True
    if recycle_runs:
        try:
            if more_info:
                res_k = np.load('aux/{}/{}_more_info.npz'.format(manipulation, str_trial))
                SpikeDots_nc = res_k['SpikeDots_nc']
                SpikeDots_c = res_k['SpikeDots_c']
                X = res_k['X']
                Xh_nc = res_k['Xh_nc']
                Xh_c = res_k['Xh_c']
                Spikes = res_k['Spikes']
            else:
                res_k = np.load('aux/{}/{}.npz'.format(manipulation, str_trial))
            Rs = res_k['Rs']
            CVs = res_k['CVs']
            Es_nc = res_k['Es_nc']
            Errors_nc = res_k['Errors_nc']
            Es_c = res_k['Es_c']
            Errors_c = res_k['Errors_c']
            total_spikes = res_k['total_spikes']
            volley_count = res_k['volley_count']
            run_k = False
        except:
            pass
    if run_k:
        print(str_trial)
        X = constant_signal(x, T_c, dt, signal_noise=signal_noise, seed=seed)
        C = lambd * X
        X = X[::step_t, :]
        F_m, Omega_m, Th_m, D_m, Lambd_m, dV, N_t, N_dead_neurons = modify_params(T_c, F, Omega, Th, Lambd, D, d, dt, seed=seed)
        N = F.shape[0]
        delay_dt = int(np.ceil((d['delay'] / 1000) / dt))
        R0_big = np.zeros(N_t)
        R0_big[:N] = R0
        V0_big = np.zeros(N_t)
        V0_big[:N] = V0
        V0_big[N:] = np.dot(F_m[N:, :], X[0, :] - np.dot(D_m, R0_big))
        if numba:
            R, Spikes, Rn, Vn = simnumba(F_m, Omega_m, Th_m, C, dt, dV, Lambd, lambd, delay_dt, R0_big, V0_big, ref_t=ref_t, jump=step_t, tmpspikefailure=d['tmp-spike-failure'], synapticnoise=d['synaptic-noise'], seed=seed)
        else:
            R, Spikes, Rn, Vn = sim(F_m, Omega_m, Th_m, C, dt, dV, Lambd, lambd, delay_dt, R0_big, V0_big, ref_t=ref_t, jump=step_t, tmpspikefailure=d['tmp-spike-failure'], synapticnoise=d['synaptic-noise'], seed=seed)

        Xh_nc = np.dot(R, D_m.T)
        Epsilon = X - Xh_nc
        Epsilon_abs = np.abs(Epsilon)
        amp = np.max(np.maximum(Epsilon, 0), axis = 0) + np.max(np.maximum(-Epsilon, 0), axis = 0)
        l2 = np.linalg.norm(Epsilon, axis=1)
        l1 = np.linalg.norm(Epsilon, ord=1, axis=1)

        Es_nc = Epsilon_abs.reshape(-1)
        Es_max_amp = np.max(amp)
        Es_mean_amp = np.mean(amp)
        Es_max = np.max(Epsilon_abs)
        Es_mean_max = np.mean(np.max(Epsilon_abs, axis=0))
        Es_mean = np.mean(Epsilon_abs)
        Es_mean_l2 = np.mean(l2)
        Es_max_l2 = np.max(l2)
        Es_mean_l1 = np.mean(l1)
        Es_max_l1 = np.max(l1)

        Epsilon_abs_dead = np.abs(X)
        amp_dead = np.max(np.maximum(X, 0), axis = 0) + np.max(np.maximum(-X, 0), axis = 0)
        l2_dead = np.linalg.norm(X, axis=1)
        l1_dead = np.linalg.norm(X, ord=1, axis=1)

        Es_max_amp_dead = np.max(amp_dead)
        Es_mean_amp_dead = np.mean(amp_dead)
        Es_max_dead = np.max(Epsilon_abs_dead)
        Es_mean_max_dead = np.mean(np.max(Epsilon_abs_dead, axis=0))
        Es_mean_dead = np.mean(Epsilon_abs_dead)
        Es_mean_l2_dead = np.mean(l2_dead)
        Es_max_l2_dead = np.max(l2_dead)
        Es_mean_l1_dead = np.mean(l1_dead)
        Es_max_l1_dead = np.max(l1_dead)

        Errors_nc = np.zeros((9, 2))
        Errors_nc[0, 0] = Es_max_amp
        Errors_nc[1, 0] = Es_mean_amp
        Errors_nc[2, 0] = Es_max
        Errors_nc[3, 0] = Es_mean_max
        Errors_nc[4, 0] = Es_mean
        Errors_nc[5, 0] = Es_mean_l2
        Errors_nc[6, 0] = Es_max_l2
        Errors_nc[7, 0] = Es_mean_l1
        Errors_nc[8, 0] = Es_max_l1
        Errors_nc[0, 1] = Es_max_amp_dead
        Errors_nc[1, 1] = Es_mean_amp_dead
        Errors_nc[2, 1] = Es_max_dead
        Errors_nc[3, 1] = Es_mean_max_dead
        Errors_nc[4, 1] = Es_mean_dead
        Errors_nc[5, 1] = Es_mean_l2_dead
        Errors_nc[6, 1] = Es_max_l2_dead
        Errors_nc[7, 1] = Es_mean_l1_dead
        Errors_nc[8, 1] = Es_max_l1_dead

        if more_info:
            SpikeDots_nc = []
            for neuron, spikes in enumerate(Spikes):
                for spike in spikes:
                    SpikeDots_nc.append(np.dot(D[:, neuron], Epsilon[spike - 1, :] / np.linalg.norm(Epsilon[spike - 1, :])))

        Xh_c = readout_correction(R, D_m, threshold)

        Epsilon = X - Xh_c
        Epsilon_abs = np.abs(Epsilon)
        amp = np.max(np.maximum(Epsilon, 0), axis = 0) + np.max(np.maximum(-Epsilon, 0), axis = 0)
        l2 = np.linalg.norm(Epsilon, axis=1)
        l1 = np.linalg.norm(Epsilon, ord=1, axis=1)

        Es_c = Epsilon_abs.reshape(-1)
        Es_max_amp = np.max(amp)
        Es_mean_amp = np.mean(amp)
        Es_max = np.max(Epsilon_abs)
        Es_mean_max = np.mean(np.max(Epsilon_abs, axis=0))
        Es_mean = np.mean(Epsilon_abs)
        Es_mean_l2 = np.mean(l2)
        Es_max_l2 = np.max(l2)
        Es_mean_l1 = np.mean(l1)
        Es_max_l1 = np.max(l1)

        Epsilon_abs_dead = np.abs(X)
        amp_dead = np.max(np.maximum(X, 0), axis = 0) + np.max(np.maximum(-X, 0), axis = 0)
        l2_dead = np.linalg.norm(X, axis=1)
        l1_dead = np.linalg.norm(X, ord=1, axis=1)

        Es_max_amp_dead = np.max(amp_dead)
        Es_mean_amp_dead = np.mean(amp_dead)
        Es_max_dead = np.max(Epsilon_abs_dead)
        Es_mean_max_dead = np.mean(np.max(Epsilon_abs_dead, axis=0))
        Es_mean_dead = np.mean(Epsilon_abs_dead)
        Es_mean_l2_dead = np.mean(l2_dead)
        Es_max_l2_dead = np.max(l2_dead)
        Es_mean_l1_dead = np.mean(l1_dead)
        Es_max_l1_dead = np.max(l1_dead)

        Errors_c = np.zeros((9, 2))
        Errors_c[0, 0] = Es_max_amp
        Errors_c[1, 0] = Es_mean_amp
        Errors_c[2, 0] = Es_max
        Errors_c[3, 0] = Es_mean_max
        Errors_c[4, 0] = Es_mean
        Errors_c[5, 0] = Es_mean_l2
        Errors_c[6, 0] = Es_max_l2
        Errors_c[7, 0] = Es_mean_l1
        Errors_c[8, 0] = Es_max_l1
        Errors_c[0, 1] = Es_max_amp_dead
        Errors_c[1, 1] = Es_mean_amp_dead
        Errors_c[2, 1] = Es_max_dead
        Errors_c[3, 1] = Es_mean_max_dead
        Errors_c[4, 1] = Es_mean_dead
        Errors_c[5, 1] = Es_mean_l2_dead
        Errors_c[6, 1] = Es_max_l2_dead
        Errors_c[7, 1] = Es_mean_l1_dead
        Errors_c[8, 1] = Es_max_l1_dead

        if more_info:
            SpikeDots_c = []
            for neuron, spikes in enumerate(Spikes):
                for spike in spikes:
                    SpikeDots_c.append(np.dot(D[:, neuron], Epsilon[spike - 1, :] / np.linalg.norm(Epsilon[spike - 1, :])))

        alive_neurons = np.where(Th_m < np.inf)[0]
        N_alive = alive_neurons.shape[0]
        if N_t - N_dead_neurons != N_alive:
            raise SystemExit(0)
        CVs = np.empty(N_alive)
        Rs = np.empty(N_alive)
        for index, i_neuron in enumerate(alive_neurons):
            length = len(Spikes[i_neuron]) # Spikes[i_neuron] is a list with the index of time steps where there was a spike for neuron # i_neuron
            Rs[index] = length / time_c # Firing rate of that neuron is the total number of spikes for the entire trial divided by the time of trial
            if length != 0:
                if length != 1:
                    if length > 3:
                        isi = np.diff(Spikes[i_neuron]) #instead of taking the ISI in units of time, we take it in units of simulation time
                        CVs[index] = np.std(isi) / np.mean(isi)
                    else:
                        CVs[index] = -1 # These neurons only spiked less than 3 spikes (so CV is prone to error)
                else:
                    CVs[index] = -2 # These neurons only spiked once (so CV is undefined)
            else:
                CVs[index] = -3 # These neurons never spiked (so CV is undefined)
        # Total number of spikes
        total_spikes = np.sum([len(i) for i in Spikes])
        # distribution of spike volleys
        spike_times, volley_count = np.unique([s_t for neuron in Spikes for s_t in neuron], return_counts=True)
        if more_info:
            np.savez('aux/' + manipulation + '/' + str_trial + '_more_info', Rs=Rs, CVs=CVs, Es_nc=Es_nc, Errors_nc=Errors_nc, Es_c=Es_c, Errors_c=Errors_c, total_spikes=total_spikes, volley_count=volley_count, SpikeDots_nc=SpikeDots_nc, SpikeDots_c=SpikeDots_c, X=X, Xh_nc=Xh_nc, Xh_c=Xh_c, Spikes=Spikes)
        else:
            np.savez('aux/' + manipulation + '/' + str_trial, Rs=Rs, CVs=CVs, Es_nc=Es_nc, Errors_nc=Errors_nc, Es_c=Es_c, Errors_c=Errors_c, total_spikes=total_spikes, volley_count=volley_count)
    if more_info:
        return Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count, SpikeDots_nc, SpikeDots_c, X, Xh_nc, Xh_c, Spikes
    else:
        return Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count



def run_trial(dim, cov, d_sim, d_seed, manipulation='baseline', value=0, recycle_runs=False, more_info=False, no_sig=2, numba=True):
    '''
        General function that runs one trial
    '''
    # Section 1: Ramping Activity
    R0, V0, F, Omega, Th, Lambd, D, x, str_ramp = ramp(dim, cov, d_sim, d_seed, recycle_runs=recycle_runs, numba=numba)
    # Section 2: Baseline Activity
    d = build_dict()
    d['noise'] = d_sim['noise']
    d['delay'] = d_sim['delay']
    Rs_b, CVs_b, Es_nc_b, Errors_nc_b, Es_c_b, Errors_c_b, total_spikes_b, volley_count_b = keep(d, R0, V0, F, Omega, Th, Lambd, D, x, str_ramp, d_sim, d_seed, recycle_runs=recycle_runs, numba=numba)
    active_neurons = np.where(Rs_b > 0)[0]
    d['active_neurons_sorted'] = active_neurons[np.argsort(Rs_b[active_neurons])]
    d['neurons_sorted'] = np.argsort(Rs_b)
    d['aligned_neurons_sorted'] = np.argsort(np.dot(F, x))
    # Section 3: Run Perturbation trial
    d[manipulation] = value
    if more_info:
        Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count, SpikeDots_nc, SpikeDots_c, X, Xh_nc, Xh_c = keep(d, R0, V0, F, Omega, Th, Lambd, D, x, str_ramp, d_sim, d_seed, manipulation=manipulation, value=value, recycle_runs=recycle_runs, more_info=True, no_sig=no_sig, numba=False)
        return Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count, SpikeDots_nc, SpikeDots_c, X, Xh_nc, Xh_c
    else:
        Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count = keep(d, R0, V0, F, Omega, Th, Lambd, D, x, str_ramp, d_sim, d_seed, manipulation=manipulation, value=value, recycle_runs=recycle_runs, no_sig=no_sig, numba=numba)
        return Rs, CVs, Es_nc, Errors_nc, Es_c, Errors_c, total_spikes, volley_count


def clean_string(string):
    return string.replace(", ", '_').replace('(', '').replace(')', '').replace('\'', '').replace('.', 'd')


# # def readout_correction(R, D, threshold):
# #     Xh_nc = np.dot(R, D.T)
# #     av_norm_xh = np.mean(np.linalg.norm(Xh_nc, axis=1))
# #     if av_norm_xh > 0:
# #         return Xh_nc * (av_norm_xh + threshold - 0.5) / av_norm_xh
# #     else:
# #         return Xh_nc

def readout_correction(R, D, threshold):
    Xh_c = np.zeros((R.shape[0], D.shape[0]))
    for i in range(R.shape[0]):
        xh = np.dot(D, R[i, :])
        if np.linalg.norm(xh) > 0:
            Xh_c[i, :] = xh * (1 + (threshold - 0.5) / np.linalg.norm(xh))
    return Xh_c