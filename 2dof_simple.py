import control.matlab as ctr
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def hit_ground(t,y):
    y[0] = 0

def f(t,y,g,m,lt,J,y_engine_fire):
    
    x, xdot, y, ydot, theta, thetadot = y
    
    thruster_switch         = 10*np.pi/180
    thruster_force          = 10
    thruster_vel_switch     = 20
    hysteresis_ang_switch   = 5 * thruster_switch
    
    
    #control law - Rotational Thruster - state machine
    if(thetadot < 0):
        if(theta > 0):                                      #hysteresis state
            if(abs(theta) > hysteresis_ang_switch):
                Ft = thruster_force * np.sign(theta)
            else: Ft = 0
        elif(theta < 0):                                    #non-hysteresis state
            if(abs(theta) > thruster_switch):
                Ft = thruster_force * np.sign(theta)
            else: Ft = 0
    elif(thetadot >= 0):
        if(theta < 0):                                      #hysteresis state
            if(abs(theta) > hysteresis_ang_switch):
                Ft = thruster_force * np.sign(theta)
            else: Ft = 0
        elif(theta > 0):                                    #non-hysteresis state
            if(abs(theta) > thruster_switch):
                Ft = thruster_force * np.sign(theta)
            else: Ft = 0
    

    
    #control law - Engine fire
    if(y < y_engine_fire and y > 0 and ydot < 0):    Fe = 200
    elif(y < 0): Fe = g
    else:                       Fe = 0
    
    #ground contraint - all of the sudden ydot becomes zero
    
    
    #though in an ideal world there is no torque from the main engine,
    #in the real world there is... So add a disturbance to the theta based 
    #main engine firing
    
    
    dydt = [
        xdot,
        np.cos(theta)*Ft/m - np.sin(theta)*Fe/m,
        ydot,
        np.sin(theta)*Ft/m + +np.cos(theta)*Fe/m - g,
        thetadot,
        -lt*Ft/J        
        ]
    
    return(dydt)

def main():
    
    x0          = 0
    xdot0       = 10
    y0          = 1000
    ydot0       = 0
    theta0      = np.pi
    thetadot0   = 0
    initial_vals = [x0,xdot0,y0,ydot0,theta0,thetadot0] 
    
    y_engine_fire = 500
    J = 0.1
    g = 9.81
    m = 10
    lt = 3
    tspan = (0,20)
    sol = sp.integrate.solve_ivp(f, tspan, initial_vals, args=(g, m, lt,J,y_engine_fire),max_step = 1e-3)
    x, xdot, y, ydot, theta, thetadot = (sol.y)
    t = sol.t

    # plt.ion
    # plt.plot(t,x)
    # plt.plot(t,y)
    fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharey=False)
    ax[0].plot(t,theta)
    ax[1].plot(t,thetadot)
    ax[2].plot(t,xdot)
    ax[3].plot(t,ydot)
    # plt.legend(['x','y','t'])
           

if __name__ == "__main__":
    main()
