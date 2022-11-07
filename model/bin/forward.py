# Start with plotting settings, need to override my defaults
import matplotlib
import matplotlib.pyplot as plt

# load the relevant portions of the scientific python stack
from math import exp
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat

# Use the random package rather than numpy.random because it is much faster for single pulls and we cannot preallocate
# Profiling indicates that generating the random numbers is a pretty significant chunk of execution time
from random import randint, gauss, random as rand
from scipy import stats
# for simple profiling
import time


def model(layerthickness, acc_t, rho_o, total_time, rho_snow, eta0, n_tsteps=200):
    """ Return the surface height of the reflector
    
    Parameters
    ---------
    layerthickness: np.ndarray(len(depth) x 1)
        The thickness of the intial layers (e, form Navarre)
    acc_t: np.ndarray(n_tsteps x 1)
        daily accumulation rate
    initial density: np.ndarray(len(depth) x 1)
        estimate for the initial density of the snow pack (assumed to be nearly linear)
    total_time: float
        total time of the model run
    rho_snow: float
        new snow density
    eta0: float
        snow viscosity
    n_tsteps: int
        number of time steps in model run

    Returns
    ---------    
    reflector height: nparray
    
    """


    timestep_size = total_time / n_tsteps
    rho_t=rho_o
    RH = np.zeros(len(acc_t))
    for step in range(n_tsteps):
        RH[step],rho_t,layerthickness = model_step(rho_t, acc_t[step],layerthickness,eta0,rho_snow=rho_snow, timestep_size=timestep_size,)

    return RH
        
def model_step(rho, acc_current, layerthickness,eta,timestep_size=86400,rho_snow=270.0, Hpole = 6.120):
    """simple model for snow densification (assuming snow is a viscous fluid).
    
    
    Parameters
    ---------
    rho: np.ndarray(variable)
        updated density profile
    layerthickness: np.ndarray(varibale)
        depth of snow parcels
    acc_current: float
        history of daily accumulation
    eta: float
        snow viscosity
    timestep_size: float
        size of timestep
    rho_snow: float
        new snow density
    Hpole: float
        total length of pole used to anchor antenna

        
        
    Returns
    ---------    
    reflector height: float
    
    """
    g = 9.8
    D_new=layerthickness.copy()
    rho_new=rho.copy()
    if acc_current != 0.0:
        D_new=np.insert(D_new,0,acc_current)
        rho_new=np.insert(rho_new,0,rho_snow)
        

    
    M = np.multiply(D_new,rho_new) 

    sig_z = np.cumsum(-g * M)
    dDD = sig_z/eta * timestep_size
    Dnew=np.multiply(1.0+dDD,D_new)
    rhonew=np.multiply(M,1.0/Dnew)

    RH = Hpole - np.sum(Dnew)

    return RH,rhonew,Dnew



def cost(modeled_surface, measured_surface, acc,rho0,eta,rhonew, regularization=0.0):
    """Return the cost function given the measured layer and our model output
    
    Parameters
    ----------
    modeled_surface: np.ndarray(n x 1)
        The model output
    measured_surface: np.ndarray(n x 1)
        The data
    regularization: float
        A scaling of the total variation to add to the misfit when determining the cost function (really just for the initial density estimate.
        This will need to be tuned if used. Turn off by setting to zero.
        
    Returns
    cost_reg: float
        The misfit/objective function/cost function
    cost_unreg: float
        The misfit/objective function/cost function without regularization. Useful for comparing regularizations
    """
    #create regularization vector for accumulation, initial density variations.
    
    # need to split the regulariztion into accumulation and velocity variations
    x_cost = np.linspace(1,200,200)
    if regularization > 0:
        #slope, intercept, r_value, p_value, std_err = stats.linregress(x_cost,modeled_surface)
        #slope_measured, intercept, r_value, p_value, std_err=  stats.linregress(x_cost,measured_surface)
        reg = abs(modeled_surface[-1]-measured_surface[-1])
    
    else:
        reg = 0.0
        
    unreg_cost = np.sqrt(np.nansum((x_cost*(modeled_surface - measured_surface))**2))
    return unreg_cost + regularization * reg, unreg_cost



def monte_carlo_run(rho0,acc0,rho_new,eta,D,RH_data,
                    total_time=200.0*86400,
                    n_tsteps=200,
                    n_iterates=2 * 10 ** 4,
                    stepsize=0.001,
                    regularization=0.0):
    """Perform a suite of runs to get our distribution. We start with an initial guess, and move around.
    
    Parameters
    ----------
    rho0: np.ndarray(lengh(depth)x1)
        Initial density
    acc0: np.ndarray(n_tstepsx1)
        initial accumulation history
    rhow_new: float
        new snow density
    total_time: float, optional
        Time over which to run the model
    n_tsteps: int, optional
        Number of degrees of freedom each for accumulation variable
    n_iterates: int, optional
        Number of accepted models to run to.
    stepsize: float, optional
        Perturbation size for each new model.
    regularization: float, optional
        Regularize the output by scaling the square total variation by this amount. Should be >=0.0, and =0 if unregularized.
    
    Returns
    -------
    outputs: np.ndarray (n_tsteps * (len(acc0)+ len(rho0)+2) x n_iterates)
        The model inputs that have been accepted. The first will be the initial guess.
    misfits: np.ndarray (n_iterates x 2)
        Cost function with and without regularization. First column is regularized, second isn't.
    """
    # this is just to check how long it took
    now = time.time()
    
    l=int(len(rho0)+len(acc0)+2)
    # Preallocate locations for output--we want all our accepted models saved
    outputs = np.zeros((l,n_iterates))
    # misfits will have both a regularized and an unregularized row
    misfits = np.zeros((n_iterates, 2))
    
    # variable to increment--we use a while loop rather than a for loop since sometimes we dont step, and thus
    # we do no know how many times we will actually execute the loop
    i = 0
    vector_inputs=np.append(rho0,acc0)
    vector_inputs=np.append(vector_inputs,rho_new)
    vector_inputs=np.append(vector_inputs,eta)
    # fencepost; this is so we don't have errors comparing to previous state
    RH = model(D,acc0,rho0,total_time,rho_new,eta, n_tsteps=n_tsteps)
    outputs[:, i] = vector_inputs
    misfits[i] = cost(RH,RH_data,acc0,rho0,eta,rho_new,regularization=regularization)
    # This is the real bulk of the work--keep trying to make steps downward
    while i < n_iterates - 1:
        # start by perturbing--first copy our vector to avoid modifying original in case we dont step
 
        vector_inputs_pert = vector_inputs.copy()

        # randint gives us a random component of the vector to perturb
        rand_ind = randint(0, len(vector_inputs)-1)
        # We use a normal distribution with variance "stepsize" for the perturbation
        if rand_ind>len(rho0)-1 and rand_ind<len(vector_inputs_pert)-2:
            vector_inputs_pert[rand_ind] = vector_inputs_pert[rand_ind] + gauss(0, stepsize*10.)
        elif rand_ind<len(rho0) or rand_ind ==len(vector_inputs_pert)-2:
            vector_inputs_pert[rand_ind] = vector_inputs_pert[rand_ind] + gauss(0, stepsize*100.)
        else:
            vector_inputs_pert[rand_ind] = vector_inputs_pert[rand_ind] + gauss(0, stepsize*1.*10e10)
            
        # enforce that accumulation density and viscosity are strictly positive
        if vector_inputs_pert[rand_ind] < 0:
            vector_inputs_pert[rand_ind] = 0.0
            # if not positive, we don't save the perturbed state and we just restart this iteration of the while loop
        
        #evaluate model
        
        rho=vector_inputs_pert[:len(rho0)]
        acc=vector_inputs_pert[len(rho0):][:len(acc0)]
        rho_new=vector_inputs_pert[len(rho0):][len(acc0):][0]
        eta=vector_inputs_pert[len(rho0):][len(acc0):][1]
        
        #RH = model(x, vx_and_acc_of_t_pert, acc_rate, vels, total_time, n_tsteps=n_tsteps)
        RH = model(D,acc,rho,total_time,rho_new,eta, n_tsteps=n_tsteps)
        # see if model is any good---store cost function in vector, will be overwritten if not chosen
        misfits[i + 1, :] = cost(RH,RH_data,acc,rho,eta,rho_new, regularization=regularization)
        
        # Decide whether to accept our new guess
        if i < 1 or rand() < min(1, exp(-(misfits[i + 1, 0] - misfits[i, 0]))):
            # We have accepted the model! We need to store it, and save this state as the new model to step from
            outputs[:, i + 1] = vector_inputs_pert
            vector_inputs = vector_inputs_pert.copy()
            
            # increment only upon acceptance of the guess
            i += 1

    print('Run took {:5.1f}s'.format(time.time() - now))
    return outputs, misfits