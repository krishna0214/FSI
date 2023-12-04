import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# Your data
x = np.array([0, 0.010827067669172935, 0.02165413533834587, 0.035187969924812046, 0.048721804511278194, 0.060902255639097735, 0.0730827067669173, 0.0825563909774436, 0.0893233082706767, 0.09744360902255639, 0.10421052631578948, 0.11097744360902256, 0.1177443609022556, 0.12586466165413535, 0.13263157894736843, 0.13804511278195492, 0.14616541353383455, 0.15699248120300752, 0.16646616541353382, 0.1786466165413534, 0.1894736842105263, 0.20030075187969926, 0.20977443609022556, 0.22195488721804513, 0.23007518796992482, 0.23954887218045112, 0.24766917293233082, 0.2557894736842105, 0.2666165413533835, 0.2760902255639098, 0.2855639097744361, 0.2923308270676692, 0.30045112781954886, 0.30992481203007516, 0.31939849624060146, 0.3275187969924812, 0.3369924812030075, 0.3464661654135339, 0.35458646616541356, 0.36406015037593986, 0.3735338345864663, 0.38300751879699246, 0.39248120300751876, 0.40601503759398505, 0.4195488721804511, 0.43443609022556384, 0.44661654135338347, 0.4615037593984962, 0.47233082706766927, 0.487218045112782, 0.5007518796992482, 0.5142857142857142, 0.5291729323308271, 0.5413533834586466, 0.5535338345864661, 0.5684210526315789, 0.5792481203007519, 0.5914285714285715, 0.6063157894736843, 0.6171428571428571, 0.6320300751879699, 0.6442105263157895, 0.6590977443609023, 0.6712781954887218, 0.684812030075188, 0.695639097744361, 0.709172932330827, 0.7227067669172933, 0.7362406015037594, 0.7497744360902256, 0.7619548872180452, 0.7754887218045113, 0.7876691729323309, 0.8025563909774437])
y = np.array([0.004, 0.0786026200873362, 0.1179039301310052, 0.1965065502183414, 0.3930131004366828, 0.6681222707423586, 0.9432314410480345, 1.3362445414847173, 1.6899563318777293, 2.082969432314412, 2.5152838427947604, 2.8689956331877724, 3.2620087336244534, 3.6550218340611345, 3.969432314410481, 4.323144104803495, 4.598253275109171, 4.794759825327512, 4.951965065502183, 4.951965065502183, 4.794759825327512, 4.598253275109171, 4.323144104803495, 3.969432314410481, 3.7336244541484724, 3.458515283842795, 3.222707423580786, 2.9475982532751104, 2.5938864628820966, 2.2794759825327517, 1.9650655021834051, 1.650655021834062, 1.4541484716157207, 1.179039301310043, 0.9432314410480345, 0.6681222707423586, 0.3144104803493448, 0.0786026200873362, -0.1572052401746724, -0.3930131004366828, -0.6288209606986896, -0.7467248908296966, -0.9039301310043673, -0.9432314410480327, -0.9432314410480327, -0.8646288209606983, -0.8253275109170328, -0.6288209606986896, -0.5502183406113534, -0.43231441048034824, -0.3144104803493448, -0.23580786026200862, -0.1572052401746724, -0.0786026200873362, 0, 0, 0.0786026200873362, 0.0786026200873362, 0.1179039301310052, 0.1179039301310052, 0.1179039301310052, 0.1179039301310052, 0.0786026200873362, 0.1179039301310052, 0.0786026200873362, 0.0786026200873362, 0.0786026200873362, 0.0786026200873362, 0, -0.0786026200873362, -0.0786026200873362, -0.0786026200873362, -0.0786026200873362, -0.1179039301310052])

# Fit a spline to the data
spline = UnivariateSpline(x, y, s=0)

def velocity(Area,time):
    y_value = spline(time)
    u_in=(((y_value*10**-3)/60)/(Area))
    return u_in