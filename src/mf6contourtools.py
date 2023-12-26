import numpy as np

def get_contour_levels(lo, hi, n):
    """Return a float array of orderly levels for contouring.
    
    Parameters
    ----------
    hi: float
        high values (max value)
    lo: float
        low value (mean value)
    n:  int
        desired (approximate) number of levels
        
    Returns
    -------
    Given a maximum and minimum values of a data set, the function
    returns an array of n orderly values that include the hi and lo
    values and which is suitable for a set of contouring levels.
    
    >>>tests = [[0, 45], [0, 123], [-45, 201], [-1, 1], [0.01, 0.6], [-0.009, 0.004], [-0.8, 0.9], [-335, -12]]
    >>>for test in tests:
    >>>   print('test: ({}) --> ticks: {}'.format(test, get_orderly_levels(*test, 10)))
    test: ([0, 45]) --> ticks: [ 0.  5. 10. 15. 20. 25. 30. 35. 40. 45. 50.]
    test: ([0, 123]) --> ticks: [  0.  20.  40.  60.  80. 100. 120. 140. 160. 180. 200.]
    test: ([-45, 201]) --> ticks: [-100.  -60.  -20.   20.   60.  100.  140.  180.  220.  260.  300.]
    test: ([-1, 1]) --> ticks: [-1.  -0.8 -0.6 -0.4 -0.2 -0.   0.2  0.4  0.6  0.8  1. ]
    test: ([0.01, 0.6]) --> ticks: [0.   0.06 0.12 0.18 0.24 0.3  0.36 0.42 0.48 0.54 0.6 ]
    test: ([-0.009, 0.004]) --> ticks: [-0.01  -0.008 -0.006 -0.004 -0.002  0.     0.002  0.004  0.006  0.008
      0.01 ]
    test: ([-0.8, 0.9]) --> ticks: [-1.  -0.8 -0.6 -0.4 -0.2 -0.   0.2  0.4  0.6  0.8  1. ]
    test: ([-335, -12]) --> ticks: [-400. -360. -320. -280. -240. -200. -160. -120.  -80.  -40.    0.]
    """
    s = 10 ** (np.floor(np.log10(hi - lo)))
    start = s * np.floor(lo / s)
    end   = s * np.ceil( hi / s)
    s = (end - start) / n
    return np.round(np.arange(start, end + s, s), decimals=6)

    
if __name__ == '__main__':
    print('Hello')
    
    tests = [[0, 45], [0, 123], [-45, 201], [-1, 1], [0.01, 0.6], [-0.009, 0.004], [-0.8, 0.9], [-335, -12]]
    for test in tests:
       print('test: ({}) --> ticks: {}'.format(test, get_orderly_levels(*test, 10)))
