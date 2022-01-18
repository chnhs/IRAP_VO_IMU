from cmath import pi
import numpy as np

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
        

    def update(self, data):
        '''Inertial Navigation in Earth-fixed frame
        '''
        # Get IMU measurements, here in b-frame
        fibb = np.array([data.ax, data.ay, data.az]).reshape(3,1)
        wibb = np.array([data.wx, data.wy, data.wz]).reshape(3,1)

        # Initial alignment
        if self.initial_pos_vel_alignment_ == False:
            self.position_velocity_init(self.pos_state, self.vel_state)
        
        if self.initial_attitude_alignment_ == False:
            self.attitude_init(self.attitude_state, wibb, fibb)
        
        print("Initialization P: ", self.pos_state)
        print("Initialization V: ", self.vel_state)
        print("Initialization A: ", self.attitude_state)
        
        # Conversions
        # Get initial Cbe from previous attitude state


        # Correction of raw data (calibration)

        # Attitude update
        # attitude change update
        self.attitude_state = self.attitude_state + self.euler_rate(wibb)

        # Transformation of specific force to navigation frame
        
        # Velocity update
        # Acceleration change
        self.vel_state = self.vel_state + self.acceleration_rate(fibb,wibb)

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
        '''The INS position and velocity must be initialized using external information
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

        self.initial_pos_vel_alignment_ = True
        

    def attitude_init(self, attitude_state0, wibb, fibb):
        '''The INS position and velocity must be initialized using external information
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


