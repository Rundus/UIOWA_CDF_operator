from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as R
import numpy as np

vec = np.array([[0,0,1],[0,0,2],[1,0,1]])
r = R.from_rotvec([0, np.pi/2,0])
print(r.apply(vec))


# key_rots = R.random(5,random_state=1)
# key_times = [0, 1, 2, 3, 4]
# slerp = Slerp(key_times, key_rots) #create interpolator object
#
# times = [0, 0.5, 0.75, 1]
# interp_rots = slerp(times)
#
# print(key_rots.as_euler('xyz',degrees=True))
# print(interp_rots.as_euler('xyz',degrees=True))