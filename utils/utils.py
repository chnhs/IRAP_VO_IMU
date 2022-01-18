# Unit conversion and setting up the data-structure
import numpy as np
import pandas as pd

class PMcar:
    length = 4.484 #meters
    width = 1.736
    height_without_tire_deflection = 2.000 #meters
    wheelbase = 2.920 #meters
    track_width = 1.330 #meters
    grounf_clereance_without_tire_deflection = 0.282 #meters
    approach_angle = 57.7 #deg
    beparture_angle = 43.6 #deg
    breakover_angle = 12.0 #deg
    tires = '315/70R17'
    tire_radius = 0.422 #m
    dist_cog_front = wheelbase*0.6 #m
    dist_cog_rear = wheelbase*0.4 #m
    h_cog = .8 #m
    nof_wheels = 4 # Unitless
    m_nom = 2000. #kg
    wheel_inertia = .4
    yaw_inertia = 2000.
    processing_mode = 1 # 1 for IMU, 2 for VO, 3 for VIO 

class Data:
    def __init__(self, df, r_nom):
        self.time = df['Time']
        self.ax = df['AxG_S1']*9.81
        self.ay = df['AyG_S1']*9.81
        self.az = df['AzG_S1']*9.81
        self.wx = df['AVx_S1']*np.pi/180
        self.wy = df['AVy_S1']*np.pi/180
        self.wz = df['AVz_S1']*np.pi/180

        self.whl_spd = np.array([df['AVy_L1']*0.10472*r_nom,
        df['AVy_R1']*0.10472*r_nom,
        df['AVy_L2']*0.10472*r_nom,
        df['AVy_R2']*0.10472*r_nom])

        self.drv_tq = np.array([df['My_Dr_L1'],
        df['My_Dr_R1'],
        df['My_Dr_L2'],
        df['My_Dr_R2']])

        self.brk_tq  = np.array([df['My_Bk_L1'],
        df['My_Bk_R1'],
        df['My_Bk_L2'],
        df['My_Bk_R2']])

        self.whl_ang_front = df['StrK_L1']*np.pi/180
        self.whl_ang_rear = df['StrK_R1']*np.pi/180

def get_ground_truth(df):
    return  {'v_x': df['sv_vx'], 'v_y': df['sv_vy'], 'v_z': df['sv_vz'],
     'roll': df['Roll_E']*np.pi/180, 'pitch': df['Pitch']*np.pi/180,
     'yaw': df['Yaw']*np.pi/180}


def parse_data(file, r_nom):
    """Create list of data suitable for sequential processing
    """
    df = pd.read_csv(file)
    data = [ Data(row[1], r_nom) for row in df.iterrows()]
    ground_truth = get_ground_truth(df)
    return data, ground_truth