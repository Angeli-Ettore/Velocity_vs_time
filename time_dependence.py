'''
taking as input a Comma Separated Values file with two columns (time, voltage) and one header, the program calculate 
the corresponding time dependence of the velocity as function of time, and plot it.
'''

import numpy as np
import pandas as pd
import csv

def read_csv(path, filename, time, voltage):
    with open(path+filename, "r", errors="ignore") as file_to_read:  #context manager with open() needed for safe closure of the file
        #rows = file_to_read.readlines()
        #headers = 0  # number of headers (starts with 0!!!)

        reader = csv.reader(file_to_read)
        headers = next(reader)

        for row in reader:  # skip the headers
            try:
                time.append(float(row[0]))
                voltage.append(float(row[1]))
            except ValueError:
                print(f"Skipping row due to invalid value: {row}") #if the value is invalid print the error and continue!
    return

def write_origin(path, filename, time, voltage, field, velocity, cum_vel, avg_vel, err_perc):
    headers = [['Time','(us)',''], ['Voltage','(V)',''], ['B','(mT)',''], ['Velocity','(m/s)',''], ['Cumulative velocity','(m/s)',''], ['Average velocity','(m/s)',''], ['Error percentage','','']]
    headers_t = [list(i) for i in zip(*headers)]
    
    data = [time, voltage, field, velocity, cum_vel, avg_vel, err_perc] 
    data_t = [list(i) for i in zip(*data)]

    data_frame = pd.DataFrame(headers_t + data_t)  
    name = path + 'output_' + filename
    print("Output file directory:", name)
    data_frame.to_csv(name, index=False, header=False)
    return

def from_voltage_to_field(conversion_factor,voltage,field):
    for i in voltage:
        field.append(conversion_factor * i)

def creep_regime(parameters, driving_force):
    return parameters[1] * np.exp(- parameters[2] * ((parameters[0]/driving_force)**0.25-1))

def depinning_transition(parameters, driving_force):
    return (parameters[1]/0.65) * parameters[2]**0.15 * ((driving_force-parameters[0])/parameters[0])**0.25

def polynomial_model():
    return

'''
def from_field_to_velocity(parameters,field,velocity):
    for f in field:
        if f > 0: #positive velocity
            if f < parameters[0]: #creep regime
                v = parameters[1] * np.exp(- parameters[2] * ((parameters[0]/f)**0.25-1))
            elif parameters[0] <= f < parameters[3]: #depinning transition
                v = (parameters[1]/0.65) * parameters[2]**0.15 * ((f-parameters[0])/parameters[0])**0.25
            elif f >= parameters[3]: #flow regime
                v = (parameters[1]/0.65) * f
        elif f < 0: #negative velocity
            if np.abs(f) < parameters[0]: #creep regime
                v = - parameters[1] * np.exp(- parameters[2] * ((parameters[0]/np.abs(f))**0.25-1))
            elif parameters[0] <= np.abs(f) < parameters[3]: #depinning transition
                v = - (parameters[1]/0.65) * parameters[2]**0.15 * ((np.abs(f)-parameters[0])/parameters[0])**0.25
            elif np.abs(f) >= parameters[3]: #flow regime
                v = - (parameters[1]/0.65) * np.abs(f)
        elif f == 0: #expected zero velocity for zero field!!
            v = 0
        velocity.append(v)'''

def from_field_to_velocity(parameters,field,velocity):
    for f in field:
        if f > 0: #positive velocity
            if f < parameters[0]: #creep regime
                v = creep_regime(parameters, f)
            elif  f >= parameters[0]: #depinning transition
                v = depinning_transition(parameters, f)
        elif f < 0: #negative velocity
            if np.abs(f) < parameters[0]: #creep regime
                v = - creep_regime(parameters, np.abs(f))
            elif np.abs(f) >= parameters[0]: #depinning transition
                v = - depinning_transition(parameters, np.abs(f))
        elif f == 0: #expected zero velocity for zero field!!
            v = 0
        velocity.append(v)

    
def cumulative_velocity_threshold(field, velocity, threshold=None): # normalized cumulative velocity
    if threshold is None or threshold < 0:
        threshold = np.mean(field)  # for ungiven or invalid valid, threshold corresponds to field mean value 

    cum_velocity = []
    vel_sum = 0
    valid_values = 0
    
    for i, f in enumerate(field):
        if f > threshold:
            valid_values += 1
            vel_sum += velocity[i]
            cum_velocity.append(vel_sum / valid_values)
        else:
            cum_velocity.append(0)
    return cum_velocity


def average_velocity_threshold(field, velocity, threshold=None):  
    if threshold is None or threshold < 0:
        threshold = np.mean(field)  # for ungiven or invalid valid, threshold corresponds to field mean value 
    
    threshold_indices = np.where(np.array(field) > threshold)[0] # indices for which field[i]>threshold
    velocity_indices = np.array(velocity)[threshold_indices] # velocities at the calculated indices
    avg_velocity = []
    
    for i, f in enumerate(field):
        if f > threshold:
            avg_velocity.append(np.mean(velocity_indices))
        else:
            avg_velocity.append(0)
        
    return avg_velocity


def error_percentage(cum_velocity, avg_velocity):
    err_velocity = []
    for cum, avg in zip(cum_velocity, avg_velocity):
        if avg > 0:
            err_velocity.append(100*(cum-avg)/avg)
        else:
            err_velocity.append(0)
    return err_velocity


def moving_average_velocity(velocity, window):
    if window <= 0 or window > len(velocity):
        raise ValueError("Invalid value of the moving average window")

    mvg_avg_velocity = []
    return mvg_avg_velocity

time=[]; voltage=[]; field=[]; velocity=[];

H_d = 0.005 #depinning field (mT)
v_d = 60 #depinning velocity (m/s)
T_d = 20000 #depinning temperature (K)
T_norm = T_d/293 #ambient temperature (K)
#H_t = 3 * H_d # Search the right value!!!
conversion_factor = 0.175 #conversion factor between applied voltage and resulting magnetic field

#parameters = [H_d, v_d, T_norm, H_t]
parameters = [H_d, v_d, T_norm]


path = 'C:/Users/Administrateur/Documents/PhD/samples/N18/osc/'
filename = 'N18_50.CSV'

read_csv(path, filename, time, voltage)

from_voltage_to_field(conversion_factor, voltage, field)

from_field_to_velocity(parameters, field, velocity)

cumulative_velocity = cumulative_velocity_threshold(field, velocity, 0.006)

average_plateau_velocity = average_velocity_threshold(field, velocity, 0.006)

error_percentage = error_percentage(cumulative_velocity, average_plateau_velocity)

print('velocity', velocity[:50])
print('cumulative velocity', cumulative_velocity[:50])
print('average plateau velocity', average_plateau_velocity[:50])

write_origin(path, filename, time, voltage, field, velocity, cumulative_velocity, average_plateau_velocity, error_percentage)

