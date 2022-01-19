from cmath import pi
from matplotlib.pyplot import findobj
import numpy as np
import utils.transformations as TF

G = 9.80665 #m/s^2

class IMU:
    initial_pos_vel_alignment_ = False
    initial_attitude_alignment_ = False


    def __init__(self, initial_state, dt):
        self.pos_state = np.array([initial_state[0],
        initial_state[1],
        initial_state[2]]).reshape(3,1)
        self.vel_state =  np.array([initial_state[3],
        initial_state[4],
        initial_state[5]]).reshape(3,1)
        self.attitude_state = np.array([initial_state[6],
        initial_state[7],
        initial_state[8]]).reshape(3,1)        
        self.dt = dt
        '''Measurement interval dt = (t)-(t-1)
        '''
        

    def update(self, data):
        '''Inertial Navigation in Navigation-frame (n-frame): {North,Eeast,Down}
           position: {lat, long, height}_b
           velocity: {vx, vy, vz}_eb^n
           attitude: Cb^n
        '''
        # Get IMU measurements, here in b-frame
        fibb = np.array([data.ax, data.ay, data.az]).reshape(3,1)
        wibb = np.array([data.wx, data.wy, data.wz]).reshape(3,1)

        # Initial alignment, attitude is in _nb
        if self.initial_pos_vel_alignment_ == False:
            self.position_velocity_init(self.pos_state, self.vel_state)

        if self.initial_attitude_alignment_ == False:
            self.attitude_init(self.attitude_state, wibb, fibb)

        # Check aligment
        # print("Initialization P: ", self.pos_state)
        # print("Initialization V: ", self.vel_state)
        # print("Initialization A: ", self.attitude_state)
        
        # Conversions
        # Get initial C from previous or initialized attitude state
        Cbn = TF.Cbn_from_euler(self.attitude_state[0],self.attitude_state[1],self.attitude_state[2])
        # Check Cbn convertion
        # One way to check: Cbn * Cnb = I
        #print("Initial Cbn: ", Cbn)

        # Correction of raw data (calibration)

        # Attitude update

        # attitude change update
        #self.attitude_state = self.attitude_state + self.euler_rate(wibb)

        # Integration
        new_Cbn = self.attitude_update(self.pos_state, self.vel_state, wibb, Cbn)

        print("CBN: ", new_Cbn)

        # Transformation of specific force to navigation frame
        fibn = self.specific_force_transformation(fibb, Cbn, new_Cbn)


        # Velocity update
        # Acceleration change
        #self.vel_state = self.vel_state + self.acceleration_rate(fibb,wibb)
        new_vel = self.velocity_update(fibn, self.vel_state, self.pos_state)


        print("VEL: ", new_vel)
        exit()

        # Position and velocity update
        self.pos_state = self.pos_state + self.vel_state*self.dt

        # Finally output navigation
        return self.pos_state, self.vel_state, self.attitude_state


    def euler_rate(self,w):
        roll_rate = w[0] + np.sin(self.attitude_state[0])*np.tan(self.attitude_state[1])*w[1] + \
            np.cos(self.attitude_state[0])*np.tan(self.attitude_state[1])*w[2]
        pitch_rate = np.cos(self.attitude_state[0])*w[1] - np.sin(self.attitude_state[0])*w[2]
        yaw_rate = np.sin(self.attitude_state[0])/np.cos(self.attitude_state[1])*w[1] + \
            np.cos(self.attitude_state[0])/np.cos(self.attitude_state[1])*w[2]

        return np.array([roll_rate, pitch_rate, yaw_rate])

    def acceleration_rate(self, a, w):
        ax_rate = a[0] + w[2]*self.vel_state[1] - w[1]*self.vel_state[2] + G*np.sin(self.attitude_state[1])
        ay_rate = a[1] - w[2]*self.vel_state[0] + w[0]*self.vel_state[2] - \
            G*np.sin(self.attitude_state[0])*np.cos(self.attitude_state[1])
        az_rate = a[2] + w[1]*self.vel_state[0] - w[0]*self.vel_state[1] - \
            G*np.cos(self.attitude_state[0])*np.cos(self.attitude_state[1])
        
        return np.array([ax_rate, ay_rate, az_rate])

    #def integrate_rates(self, )

    def position_velocity_init(self, pos0, vel0):
        '''The INS position and velocity must be initialized using external information.
            Results are in local-frame: pos {lat,long,h} vel {vx,vy,vz}_eb^n.
        '''
        if np.all(self.pos_state) == 0:
            pos0[0] = self.pos_state[0]
            pos0[1] = self.pos_state[1]
            pos0[2] = self.pos_state[2]
        else:
            # Add a flag to detec other systems, such as GNSS_FLAG = True for instance to initialize pos and vel
            pos0[0] = 0.0
            pos0[1] = 0.0
            pos0[2] = 0.0

        
        if np.all(self.vel_state) == 0:
            vel0[0] = self.vel_state[0]
            vel0[1] = self.vel_state[1]
            vel0[2] = self.vel_state[2]
        else:
            # Add a flag to detec other systems, such as GNSS_FLAG = True for instance to initialize pos and vel
            vel0[0] = 0.0
            vel0[1] = 0.0
            vel0[2] = 0.0

            # If vel from GPS, convert it from _eb^e to _eb^n using equation (5.49)
            vel0 = np.matmul(TF.Cen_rotation(pos0[0,0], pos0[1,0]), vel0)


        self.initial_pos_vel_alignment_ = True
        

    def attitude_init(self, attitude_state0, wibb, fibb):
        '''Attitude is initialized using IMU measurements. The attitude after is in _nb.
           Reminder: The INS position and velocity must be initialized using external information.
        '''
        if np.all(self.attitude_state) == 0:
            # Intialize roll and pitch with levelling
            self.pitch_roll_init_levelling(fibb, attitude_state0[0],attitude_state0[1])
            
            #Initialize yaw
            attitude_state0[2] = self.yaw_init(wibb, attitude_state0[0], attitude_state0[1])

        else:
            # Or, whatever was in the previous state
            # Add a flag to detec other systems, such as GNSS_FLAG = True for instance to initialize pos and vel
            attitude_state0[0] = self.attitude_state[0]
            attitude_state0[1] = self.attitude_state[1]
            attitude_state0[2] = self.attitude_state[2]
            # attitude_state0[0] = 0.0
            # attitude_state0[1] = 0.0
            # attitude_state0[2] = 0.0

        self.initial_attitude_alignment_ = True


    def pitch_roll_init_levelling(self, fibb, roll, pitch):
        '''Pitch and roll initialization using accelerometers when INS is stationary (or constant velocity)
        '''
        pitch = np.arctan(fibb[0], np.sqrt((fibb[1])**2+(fibb[2])**2) )
        roll = np.arctan2(-fibb[1], -fibb[2])
        return pitch, roll

    def yaw_init(self, wibb, roll, pitch):
        cos_yaw = wibb[0]*np.cos(pitch)+ wibb[1]*np.sin(roll)*np.sin(pitch)+wibb[2]*np.cos(roll)*np.sin(pitch)
        sin_yaw = -wibb[1]*np.cos(roll) + wibb[2]*np.sin(roll)

        return np.arctan2(sin_yaw, cos_yaw)


    def attitude_update(self, pos, vel, wibb, prev_Cbn):
        '''Update Cbn using angular-rates, position_ebn, and velocity_ebn{NED} solution
        '''
        Skew_ien = TF.we * np.array([ 0, np.sin(pos[0,0]), 0,
                                      -np.sin(pos[0,0]), 0, -np.cos(pos[0,0]),
                                      0,      np.cos(pos[0,0]), 0 ]).reshape(3,3)
        w_enn = np.array([ vel[1,0]/(TF.N_geo(pos[0,0])+pos[2,0]), 
                          -vel[0,0]/(TF.M_geo(pos[0,0])+pos[2,0]),
                          -vel[1,0]*np.tan(pos[0,0])/(TF.N_geo(pos[0,0])+pos[2,0]) ]).reshape(3,1)
        
        # Cbn integration, page 178 Groves(2013) *****************DOUBLE CHECK THE MATRIX MULTIPLICATION ************************************************
        return (np.matmul(prev_Cbn,(np.identity(3) + (TF.vec_to_skew(wibb))*self.dt)) - np.matmul((Skew_ien + (TF.vec_to_skew(w_enn))),prev_Cbn*self.dt))

    def specific_force_transformation(self, fibb, prev_Cbn, Cbn):
        '''Transform specific force (fibb) into the same frame of integration (fibn)
        '''
        
        return (0.5*(np.matmul((prev_Cbn+Cbn), fibb)))

    def velocity_update(self, fibn, prev_vebn, pos):
        '''Update vebn using previous velocity, specific force transformed
        '''

        Skew_ien = TF.we * np.array([ 0, np.sin(pos[0,0]), 0,
                                      -np.sin(pos[0,0]), 0, -np.cos(pos[0,0]),
                                      0,      np.cos(pos[0,0]), 0 ]).reshape(3,3)
        w_enn = np.array([ prev_vebn[1,0]/(TF.N_geo(pos[0,0])+pos[2,0]), 
                          -prev_vebn[0,0]/(TF.M_geo(pos[0,0])+pos[2,0]),
                          -prev_vebn[1,0]*np.tan(pos[0,0])/(TF.N_geo(pos[0,0])+pos[2,0]) ]).reshape(3,1)

 
        # Check if vel is really in vebn, if it is from GPS initialization
        return (prev_vebn + fibn*self.dt + (TF.gravity_to_local(pos[0,0], pos[2,0]) - np.matmul((TF.vec_to_skew(w_enn) + 2*Skew_ien ),prev_vebn) )*self.dt)
    
    def position_update(self, prev_pos, prev_vel, current_vel):
        '''Update position {lat, long, h} using current velocity, previous velocity and position
        '''
        pass
        return new_pos


