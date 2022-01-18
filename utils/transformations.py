from cv2 import sqrt
import numpy as np

# All angles are in radians

# Constants *****************************************************************************
we = 7.292115E-5 # Earth rotation rate rad/s
# Earth-rotation vector resolved in an ECI or ECEF frame
wiee_vector = wiei_vector = np.array([0, 0, we]).reshape(3,1) 

# Geodetic constants and quantities *********************************************************
# WGS84 ellipsoid
a_geo = 6378137.0 # Earth semi major axis (m)
f_geo = 1/298.257223563 # Earth flatening - For GRS80: 1/298.257222101
e2_geo = 2*f_geo - f_geo*f_geo # Squared of the first eccentricity of the ellipsoid

def N_geo(lat):
    """Radius of curvature of the ellipsoid in the prime vertical plane
    """
    return a_geo / np.sqrt(1 - e2_geo*(np.sin(lat)*np.sin(lat)))

def M_geo(lat):
    """Radius of curvature of the ellipsoid in the prime vertical plane,
    """
    return a_geo*(1 - e2_geo) / (1 - e2_geo*(np.sin(lat))^2)^(3/2) 


# Rotations *****************************************************************************

def R1(phi):
    """Rotation about x-axis by the angle phi (roll), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule).
    """
    return np.mat([[1.0, 0.0, 0.0], 
                    [0.0, np.cos(phi), np.sin(phi)],
                    [0.0, -np.sin(phi), np.cos(phi)] ])
def R2(theta):
    """Rotation about y-axis by the angle theta (pitch), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule).
    """
    return np.mat([[np.cos(theta), 0.0, -np.sin(theta)], 
                    [0.0, 1.0, 0.0 ],
                    [np.sin(theta), 0.0, np.cos(theta)] ])
def R3(yaw):
    """Rotation about z-axis by the angle phi (yaw), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule).
    """
    return np.mat([[np.cos(yaw), np.sin(yaw), 0.0], 
                    [-np.sin(yaw), np.cos(yaw), 0.0 ],
                    [0.0,     0.0,       1.0] ])

def wien(lat):
    """Earth-rotation vector resolved into local navigation frame
    """
    return np.array([we*np.cos(lat), 0.0, -we*np.sin(lat)]).reshape(3,1)

# Inertial and Earth Frames
def Cie_rotation(t_since_start):
    """ Inertial to Earth Frames rotation matrix
    """
    return np.array([np.cos(we*t_since_start), np.sin(we*t_since_start), 0.0,
                    -np.sin(we*t_since_start), np.cos(we*t_since_start), 0.0,
                     0.0,                      0.0,                          1]).reshape(3,3)

def Cei_rotation(t_since_start):
    """ Earth to Inertial Frames rotation matrix
    """
    return Cie_rotation(t_since_start).transpose()

# Earth and Local Navigation Frames
def Cne_rotation(lat, long):
    """Local Navigation to Earth Frames rotation matrix
    """
    return np.array([-np.sin(lat)*np.cos(long), -np.sin(long), -np.cos(lat)*np.cos(long),
                    -np.sin(lat)*np.sin(long), np.cos(long),   -np.cos(lat)*np.sin(long),
                     np.cos(lat),              0.0,              -np.sin(lat)]).reshape(3,3)

def Cen_rotation(lat, long):
    """Earth to Local Navigation Frames rotation matrix
    """
    return Cne_rotation(lat,long).transpose()


def Cnb_from_euler(phi, theta, yaw):
    """Local to Body Navigation Frames rotation matrix from Euler {roll, pitch, yaw}
    """
    return np.array([np.cos(theta)*np.cos(yaw), np.cos(theta)*np.sin(yaw), - np.sin(theta),
             -np.cos(phi)*np.sin(yaw)+ np.sin(phi)*np.sin(theta)*np.cos(yaw),
             np.cos(phi)*np.cos(yaw)+ np.sin(phi)*np.sin(theta)*np.sin(yaw), np.sin(phi)*np.cos(theta),
             np.sin(phi)*np.sin(yaw)+np.cos(phi)*np.sin(theta)*np.cos(yaw),
             -np.sin(phi)*np.cos(yaw)+np.cos(phi)*np.sin(theta)*np.sin(yaw), np.cos(phi)*np.cos(theta) ]).reshape(3,3)

def Cbn_from_euler(phi, theta, yaw):
    """Body to Local Navigation Frames rotation matrix from Euler {roll, pitch, yaw}
    """
    return Cnb_from_euler(phi, theta, yaw).transpose()
    
def vec_to_skew(vector):
    """Transform a vector |omega ^| into a symmetric matrix
    """
    return np.array([[0.0, -vector[2], vector[1]], 
                    [vector[2], 0, -vector[0] ],
                    [-vector[1], vector[0], 0] ])

def Cbn_to_euler(Cbn):
    '''Euler {roll, pitch, yaw} from Local Navigation rotation matrix
    '''
    #print("TESTING: ", Cbn[2,1], Cbn[2,2])
    return np.array([np.arctan2(Cbn[2,1],Cbn[2,2]), -np.arcsin(Cbn[2,0]), np.arctan2(Cbn[1,0],Cbn[0,0]) ]).reshape(3,1)


# Geodetic transformation *************************************************************

# def cart2geod(x,y,z):
#     lat0 = np.arctan( (z /sqrt(x*x + y*y) ) * (1 / (1- e2_geo) ) )

#     lat = np.arctan( (z /sqrt(x*x + y*y) ) * (1 + (e2_geo*N_geo(lat0)) ) )
#     long = np.arctan(y/x)
#     h    = sqrt(x*x + y*y) / np.cos(lat) - N_geo(lat)
#     return np.array([lat, long, h]).reshape(3,1)

# Converts geodetic coordinates {lat, long, h} into Cartesian coordinates {X,Y,Z}
def geod2cart(lat, long, h):
    """Converts geodetic coordinates {lat, long, h} into Cartesian coordinates {X,Y,Z} 
    """
    return np.array( [(N_geo(lat)+h)*np.cos(lat)*np.cos(long), 
             (N_geo(lat)+h)*np.cos(lat)*np.sin(long), 
             (N_geo(lat)*(1 - e2_geo) + h)*np.sin(lat)] ).reshape(3,1)
