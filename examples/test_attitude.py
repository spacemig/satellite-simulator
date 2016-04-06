from numpy import array,pi
from mod_attitude import *

euler_0 = array([300,-10,150])*pi/180
q_0 = quaternion_from_euler(euler_0)
#print q_0

euler_1 = euler_from_quaternion(q_0)*180/pi
print euler_1