import numpy as np

import imu.imu as imu
import vo.vo as vo
import vo_imu_fusion.vo_imu_fusion as fusion


class VioCore:
    def __init__(self, mode):
        self.mode = mode        
        

    def run(data):
        nof_samples = len(data)

        output = {'pos_x': np.zeros(nof_samples), 'pos_y': np.zeros(nof_samples), 'pos_z': np.zeros(nof_samples),
        'v_x': np.zeros(nof_samples), 'v_y': np.zeros(nof_samples), 'v_z': np.zeros(nof_samples),
        'roll': np.zeros(nof_samples), 'pitch': np.zeros(nof_samples), 'yaw': np.zeros(nof_samples)
        }

        previous_time = data[0].time
        initial_state = np.zeros(9).reshape(9,1)

        for idx, d_idx in enumerate(data):

            dt = d_idx.time - previous_time
            
            # Run IMU navigation equations
            imu_nav = imu.IMU(initial_state, dt)
            pos_state, vel_state, attitude_state = imu_nav.update(d_idx)

            # Run vo
            

            print("POS STATES: ", pos_state)
            print("VEL STATES: ", vel_state)
            print("ATT STATES: ",attitude_state)

            #exit()

            # Run vo+imu


            # Record
            # output['pos_x'][:, [idx]] = tire_force_estimator.f_x.reshape(-1,1)
            # output['pos_y'][:, [idx]] = tire_force_estimator.f_y.reshape(-1,1)
            # output['pos_z'][:, [idx]] = tire_force_estimator.f_z.reshape(-1,1)
            # output['v_x'][idx] =  vx_for_ecu.v_x
            # output['v_x_raw'][idx] = vx_raw_estimator.v_x
            # output['m'][idx] = vehicle_mass_estimator.m

            initial_state = (pos_state[0],pos_state[1],pos_state[2], vel_state[0], vel_state[1],vel_state[2],\
                attitude_state[0], attitude_state[1], attitude_state[2])
            previous_time = d_idx.time

        return output


