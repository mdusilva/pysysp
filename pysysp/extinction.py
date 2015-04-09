import numpy as np

"""Definitions of extinction laws"""

def cardelli(wavelength, A=0., Rv=3.1):
    """
    The extinction curve of Cardelli, Clayton & Mathis (1989). 
    The wavelength input should be given in Angstron
    """
    ans = []
    for w in wavelength:
        X = w * 1.e-4  #convert from angstron to micrometer
        X = 1. / X
        if X < 0.3:
            raise ValueError('Law not defined for wavelength larger than %f' % (1./0.3*1.e4))    
        elif X >= 0.3 and X < 1.1:  #infrared
            a = 0.574 * X**1.61
            b = -0.527 * X**1.61
        elif X >= 1.1 and X < 3.3:  #optical, NIR
            y = X - 1.82
            a = 1. + 0.17699 * y - 0.50447 * y**2. -0.02427 * y**3. + 0.72085 * y**4. \
                + 0.01979 * y**5. - 0.77530 * y**6. + 0.32999 * y**7.
            b = 1.41338 * y + 2.28305 * y**2. + 1.07233 * y**3. - 5.38434 * y**4. \
                - 0.62251 * y**5. + 5.30260 * y**6. - 2.09002 * y**7.
        elif X >= 3.3 and X < 8.:    #UV and Far-UV
            if X >= 5.9:
                Fa = -0.04473 * (X - 5.9)**2. - 0.009779 * (X - 5.9)**3.
                Fb = 0.2130 * (X - 5.9)**2. - 0.1207 * (X - 5.9)**3.
            else:
                Fa = 0.
                Fb = 0.
            a = 1.752 - 0.316 * X - 0.104 / ((X - 4.67)**2. + 0.341) + Fa
            b = -3.090 + 1.825  * X + 1.206 / ((X - 4.62)**2. + 0.263) + Fb
        elif X >= 8. and X < 10.:   #Far-UV
            a = -1.073 - 0.628 * (X - 8.) + 0.137 * (X - 8.)**2. - 0.070 * (X - 8.)**3.
            b = 13.670 + 4.257 * (X - 8.) - 0.420 * (X - 8.)**2. + 0.374 * (X - 8.)**3.
        else:
            raise ValueError('Law not defined for wavelength smaller than %f' % (1./10.*1.e4))    
        
        ans.append(A * (a + b / Rv))
        
    return np.array(ans)