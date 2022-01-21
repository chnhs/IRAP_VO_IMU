from cv2 import sqrt
import numpy as np

# All angles are in radians

# Constants *****************************************************************************

we = 7.292115E-5 
'''Earth rotation rate in (rad/s)
'''

wiee_vector = wiei_vector = np.array([0, 0, we]).reshape(3,1) 
''' Earth-rotation vector resolved in an ECI or ECEF frame {3x1}
'''


# Geodetic constants and quantities *********************************************************
# WGS84 ellipsoid
a_geo = 6378137.0
'''WGS84 semi-major axis in (m)
'''
b_geo = 6356752.31425; 
'''WGS84 semi-minor axis in (m)
'''

f_geo = 1/298.257223563 
'''Earth flatening - For GRS80: 1/298.257222101 (uniteless)
'''
e2_geo = 2*f_geo - f_geo*f_geo 
'''Squared of the first eccentricity of the ellipsoid (uniteless)
'''

e_geo = 0.0818191908425 
'''WGS84 eccentricity
'''

mu = 3.986004418E14; 
'''WGS84 Earth gravitational constant (m^3 s^-2)
'''

def N_geo(lat):
    """Radius of curvature of the ellipsoid in the prime vertical plane in (m).
       Groves(2013) uses RE
    """
    return a_geo / np.sqrt(1 - e2_geo*(np.sin(lat)*np.sin(lat)))

def M_geo(lat):
    """Meridian radius of curvature (m)
       Groves(2013) uses RN
    """
    return a_geo*(1 - e2_geo) / (1 - e2_geo*(np.sin(lat))**2)**(3/2) 


# Rotations *****************************************************************************

def R1(phi):
    """Rotation about x-axis by the angle phi (roll), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule). 
    Rout{3x3}
    """
    return np.mat([[1.0, 0.0, 0.0], 
                    [0.0, np.cos(phi), np.sin(phi)],
                    [0.0, -np.sin(phi), np.cos(phi)] ])
def R2(theta):
    """Rotation about y-axis by the angle theta (pitch), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule). 
    Rout{3x3}
    """
    return np.mat([[np.cos(theta), 0.0, -np.sin(theta)], 
                    [0.0, 1.0, 0.0 ],
                    [np.sin(theta), 0.0, np.cos(theta)] ])
def R3(yaw):
    """Rotation about z-axis by the angle phi (yaw), positive in the
    counterclockwise sense as viewed along the axis toward the origin (right-hand rule). 
    Rout{3x3}
    """
    return np.mat([[np.cos(yaw), np.sin(yaw), 0.0], 
                    [-np.sin(yaw), np.cos(yaw), 0.0 ],
                    [0.0,     0.0,       1.0] ])

def wien(lat):
    """Earth-rotation vector resolved into local navigation frame {3x1}
    """
    return np.array([we*np.cos(lat), 0.0, -we*np.sin(lat)]).reshape(3,1)

# Inertial and Earth Frames
def Cie_rotation(t_since_start):
    """ Inertial to Earth Frames rotation matrix {3x3}
    """
    return np.array([np.cos(we*t_since_start), np.sin(we*t_since_start), 0.0,
                    -np.sin(we*t_since_start), np.cos(we*t_since_start), 0.0,
                     0.0,                      0.0,                          1]).reshape(3,3)

def Cei_rotation(t_since_start):
    """ Earth to Inertial Frames rotation matrix {3x3}
    """
    return Cie_rotation(t_since_start).transpose()

# Earth and Local Navigation Frames
def Cne_rotation(lat, long):
    """Local Navigation to Earth Frames rotation matrix {3x3}
    """
    return np.array([-np.sin(lat)*np.cos(long), -np.sin(long), -np.cos(lat)*np.cos(long),
                    -np.sin(lat)*np.sin(long), np.cos(long),   -np.cos(lat)*np.sin(long),
                     np.cos(lat),              0.0,              -np.sin(lat)]).reshape(3,3)

def Cen_rotation(lat, long):
    """Earth to Local Navigation Frames rotation matrix {3x3}
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
    """Transform a vector {3x1} into a symmetric matrix {3x3}
    """
    return np.array([[0.0, -vector[2,0], vector[1,0]], 
                    [vector[2,0], 0, -vector[0,0] ],
                    [-vector[1,0], vector[0,0], 0] ]).reshape(3,3)

def Cbn_to_euler(Cbn):
    '''Euler {roll, pitch, yaw} from Local Navigation rotation matrix {3x3}
    '''
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

def gravity_to_local(lat, h):
    '''Calculates acceleration due to gravity in Local coordinates {N,E,D}
    '''
    # Calculate surface gravity using the Somigliana model, (2.134)
    sinsqL = np.sin(lat)**2
    g_0 = 9.7803253359 * (1 + 0.001931853 * sinsqL) / np.sqrt(1 - e2_geo * sinsqL)

    # Calculate gravity vector using (2.140) and (2.139)
    return np.array([-8.08E-9*h*np.sin(2 * lat), 
                    0.0,
                    g_0 * (1 - (2 / a_geo) * (1 + f_geo * (1 - 2 * sinsqL) + \
        (we**2 * a_geo**2 * b_geo / mu)) * h + (3 * h**2 / a_geo**2))]).reshape(3,1)
