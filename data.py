"""
Atomic data for common metal line multiplets are defined here. You can
create additional multiplets using the Multiplets class.
"""
from .state import State
from .transitions import Multiplet

__all__ = ["atomic_data"]

atomic_data = {}

#Na D
atomic_data['Na D'] = Multiplet(
    Upper_states = [
        State(16973.37, 3/2, 1, 1/2),
        State(16956.17, 1/2, 1, 1/2),
    ],
    Lower_states = [
        State(0.0, 1/2, 0, 1/2),
    ],
    log_gf = {
        (1/2, 3/2) :  0.108,
        (1/2, 1/2) : -0.194,
    }
)

#Na 8200AA
atomic_data['Na 82'] = Multiplet(
    Upper_states = [
        State(29172.88, 3/2, 2, 1/2),
        State(29172.84, 5/2, 2, 1/2),
    ],
    Lower_states = [
        State(16973.37, 3/2, 1, 1/2),
        State(16956.17, 1/2, 1, 1/2),
    ],
    log_gf = {
        (1/2, 3/2) :  0.237,
        (3/2, 3/2) : -0.462,
        (3/2, 5/2) :  0.492,
    }
)

#Mg-b
atomic_data['Mg b'] = Multiplet(
    Upper_states = [
        State(41197.403, 1, 0, 1),
    ],
    Lower_states = [
        State(21850.405, 0, 1, 1),
        State(21870.464, 1, 1, 1),
        State(21911.178, 2, 1, 1),
    ],
    log_gf = {
        (0, 1) : -0.931,
        (1, 1) : -0.450,
        (2, 1) : -0.239
    }
)

#AlI 3950AA
atomic_data['Al i'] = Multiplet(
    Upper_states = [
        State(25347.756, 1/2, 0, 1/2),
    ],
    Lower_states = [
        State(  0.000, 1/2, 1, 1/2),
        State(112.061, 3/2, 1, 1/2),
    ],
    log_gf = {
        (1/2, 1/2) : -0.623,
        (3/2, 1/2) : -0.323,
    }
)

#KI doublet
atomic_data['KI'] = Multiplet(
    Upper_states = [
        State(13042.90, 3/2, 1, 1/2),
        State(12985.19, 1/2, 1, 1/2),
    ],
    Lower_states = [
        State(0.0, 1/2, 0, 1/2),
    ],
    log_gf = {
        (1/2, 3/2) :  0.149,
        (1/2, 1/2) : -0.154,
    }
)

#Ca i resonance line
atomic_data['Ca 4227'] = Multiplet(
    Upper_states = [
        State(23652.304, 1, 1, 0),
    ],
    Lower_states = [
        State(0.0, 0, 0, 0),
    ],
    log_gf = {
        (0, 1) : 0.265,
    }
)

#Ca HK
atomic_data['Ca HK'] = Multiplet(
    Upper_states = [
        State(25414.40, 3/2, 1, 1/2),
        State(25191.51, 1/2, 1, 1/2),
    ],
    Lower_states = [
        State(0.0, 1/2, 0, 1/2),
    ],
    log_gf = {
        (1/2, 3/2) :  0.092,
        (1/2, 1/2) : -0.213,
    }
)

#Ca triplet
atomic_data['Ca triplet'] = Multiplet(
    Upper_states = [
        State(25414.40, 3/2, 1, 1/2),
        State(25191.51, 1/2, 1, 1/2),
    ],
    Lower_states = [
        State(13710.88, 5/2, 2, 1/2),
        State(13650.19, 3/2, 2, 1/2),
    ],
    log_gf = {
        (3/2, 3/2) : -1.429,
        (5/2, 3/2) : -0.476,
        (3/2, 1/2) : -0.736,
    }
)

#CrI triplet
atomic_data['Cr i'] = Multiplet(
    Upper_states = [
        State(26801.9009, 1, 1, 2),
        State(26796.2691, 2, 1, 2),
        State(26787.4640, 3, 1, 2),
    ],
    Lower_states = [
        State(7593.4184, 2, 0, 2),
    ],
    log_gf = {
        (2, 1) : -0.208,
        (2, 2) :  0.019,
        (2, 3) :  0.158,
    }
)

#Fe 4300AA 3F->3G
atomic_data['Fe 4300'] = Multiplet(
    Upper_states = [
        State(35379.208, 5, 4, 1),
        State(35767.564, 4, 4, 1),
        State(36079.372, 3, 4, 1),
    ],
    Lower_states = [
        State(11976.239, 4, 3, 1),
        State(12560.934, 3, 3, 1),
        State(12968.554, 2, 3, 1),
    ],
    log_gf = {
        (4, 3) : -2.104, #4148.84
        (4, 4) : -0.708, #4203.21
        (3, 3) : -0.714, #4251.98
        (4, 5) :  0.164, #4272.96
        (3, 4) : -0.072, #4309.11
        (2, 3) : -0.006, #4326.98
    }
)

#Fe 4400AA 3F->5G
atomic_data['Fe 4400'] = Multiplet(
    Upper_states = [
        State(34782.421, 5, 4, 2), #There is a J=6 state but can't transition to
        State(35257.324, 4, 4, 2), #lower states since Delta J = 0, ±1
        State(35611.625, 3, 4, 2),
        State(35856.402, 2, 4, 2),
    ],
    Lower_states = [
        State(11976.239, 4, 3, 1),
        State(12560.934, 3, 3, 1),
        State(12968.554, 2, 3, 1),
    ],
    log_gf = {
        (4, 3) : -3.427, #4230.94
        (4, 4) : -1.110, #4295.33
        (3, 3) : -1.695, #4338.27
        (2, 2) : -2.886, #4369.13
        (4, 5) :  0.200, #4384.78
        (3, 4) : -0.142, #4405.99
        (2, 3) : -0.615, #4416.36
        (3, 2) : -5.000, #dummy line (there is no 3->2 transition)
    }
)

#Fe 5300AA 5F->5D
atomic_data['Fe 5300'] = Multiplet(
    Upper_states = [
        State(25899.989, 4, 2, 2),
        State(26140.179, 3, 2, 2),
        State(26339.696, 2, 2, 2),
        State(26479.381, 1, 2, 2),
        State(26550.479, 0, 2, 2),
    ],
    Lower_states = [
        State(6928.268, 5, 3, 2),
        State(7376.764, 4, 3, 2),
        State(7728.060, 3, 3, 2),
        State(7985.785, 2, 3, 2),
        State(8154.714, 1, 3, 2),
    ],
    log_gf = {
        (5, 4) : -1.321, #5271.00
        (4, 3) : -1.466, #5328.52
        (3, 2) : -1.645, #5372.98
        (4, 4) : -1.993, #5398.63
        (2, 1) : -1.844, #5407.27
        (3, 3) : -1.879, #5431.20
        (1, 0) : -2.122, #5436.03
        (2, 2) : -1.914, #5448.43
        (1, 1) : -2.091, #5457.13
        (1, 2) : -2.849, #5499.04
        (3, 4) : -3.047, #5502.99
        (2, 3) : -2.797, #5508.31
    }
)
