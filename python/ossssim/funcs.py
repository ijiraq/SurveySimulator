"""
Functions us in analysis of OSSOS H distribution.

"""
import logging

import numpy
from scipy import stats
ln10 = numpy.log(10.)

# x = numpy.array([8.66,12,17.0])
# From Kavelaars et al. (2021)
# alpha = 3/5 alpha_SI => alpha_SI = 5*alpha/3
K_ALPHA_SI = 5 * 0.4 / 3
K_BETA_SI = 0.42
K_Ho = -2.6
K_Hb = 8.1

KAVELAARS_VARIABLY_TAPERED = {'Ho': -2.6, 'Hb': 8.1, 'alpha_SI': 0.67, 'beta_SI': 0.42}

def long_lat_r(x: numpy.array, y: numpy.array, z: numpy.array) -> {numpy.array, numpy.array, numpy.array}:
    """
    Return spherical long/lat/r given a cartesian position.
    """
    r = (x**2 + y**2 + z**2)**0.5
    long = numpy.arctan2(y, x)
    lat = numpy.arcsin(z / r)
    return long, lat, r


def broken_plaw(H, a1, a2, Ho, Hb):
    """
    LF from Fraser et al. 2014
    :param H: h magnitude
    :param a1: log slope bright of break)
    :param a2: log slope faint of break)
    :param Ho: normalization of bright component
    :param Hb: normalization of faint component
    :return: N

    As presented here this the nominal cummulative form of the function.
    """
    N = 10**(a1*(H-Ho))
    N[H > Hb] = 10**(a2*(H[H > Hb]-Ho)+(a1-a2)*(Hb-Ho))
    return N


def HBO(alpha_SI, beta_SI, A, B):
    """
    Convert the 'A' and 'B' parameters in variably tappered funciton to H1 and H2.
    :param alpha_SI:
    :param beta_SI:
    :param A:
    :param B:
    :return: H1, H2

    The functional form used in the analysis users constants 'A' and 'B' rather than 10**H1 and 10**H2 as the fits are more well
    behaved.  This function converts from A/B to the H normalization.
    """
    return -5/(3*alpha_SI) * numpy.log10(A), 5/(3*beta_SI) * numpy.log10(B)


def variably_tapered(h, A, B, alpha_SI, beta_SI):
    """
    Exponentially tappered powerlaw size distribution re-expressed in H magnitude space.
    :param A: The normalization of the asymtoptic exponential
    :param B: The point where the exponential taper begins.
    :param alpah_SI: exponential asymptotic value.
    :param beta_SI: exponent of the taper.

    Derived from Schmit et al.  This form of the LF mimics behaviour seen in the streaming instability planetesimal formation process.
    """
    return A*10.**(alpha_SI*3.*h/5.)*numpy.exp(-B*10.**(-beta_SI*3.*h/5.))


def variably_tapered2(h, Ho, Hb, alpha_SI, beta_SI):
    """
    Exponentially tapered exponential size distribution re-expressed in H magnitude space.
    """
    return 10.**(alpha_SI*3./5.*(h-Ho))*numpy.exp(-10.**(-beta_SI*3./5*(h-Hb)))


def variably_tapered2_diff(h, Ho=K_Ho, Hb=K_Hb, alpha_SI=K_ALPHA_SI, beta_SI=K_BETA_SI):
    """
    The derivative of the version 2 of the variably tapered function .

    """
    part1 = numpy.exp(-numpy.float_power(10, -0.6 * beta_SI * (h - Hb)))
    part1 *= beta_SI * 1.38155 * numpy.float_power(10, 0.6 * alpha_SI * (h - Ho) - 0.6 * beta_SI * (h - Hb))
    part2 = alpha_SI * 1.38155 * numpy.float_power(10, 0.6 * alpha_SI * (h - Ho))
    part2 *= numpy.exp(-numpy.float_power(10, -0.6 * beta_SI * (h - Hb)))
    return part1+part2


def variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI):
    """
    Differential form of the variably_tapere
    :param h:
    :param Dp:
    :param E:
    :param alpha_SI:
    :param beta_SI:
    :return:
    """

    return (3./5.*ln10)*Dp*10.**(alpha_SI*3.*h/5.)*(alpha_SI+beta_SI*E*10.**(-beta_SI*3.*h/5.))*numpy.exp(-E*10.**(-beta_SI*3.*h/5.))


def log_variably_tapered(h, Dp, E, alpha_SI, beta_SI):
    return numpy.log10(variably_tapered(h, Dp, E, alpha_SI, beta_SI))


def log_variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI):
    return numpy.log10(variably_tapered_diff(h, Dp, E, alpha_SI, beta_SI))


def likelihood(H, bias, A=0.88, B=450.0,  alpha_SI=0.67, beta_SI=0.65):
    """

    :param H: Vector of measured H values
    :param bias: Vector of bias on detection of a given object


    :return:
    """
    x = numpy.linspace(H.min()-3.0, H.max())
    dx = x[1]-x[0]
    cdf = variably_tapered(x, A, B, alpha_SI, beta_SI)
    pdf = numpy.diff(cdf)/dx
    ht = numpy.linspace(H.min(), H.max())
    dh = ht[1]-ht[0]
    dy = bias*numpy.interp(H, x[1:]-dx/2.0, pdf)
    dN = numpy.interp(ht, H, dy)
    if False:
        plt.clf()
        plt.plot(H, numpy.interp(H, x, cdf), '-')
        plt.plot(x[1:], (dx * pdf).cumsum(), ':k')
        plt.plot(ht, dh * dN.cumsum(), ':k')
        plt.plot(H, numpy.ones(len(H)).cumsum(), '-')
        plt.yscale('log')
        plt.ylim(1, 10000)
        plt.show()

    N = dh*dN.sum()
    l = -N + numpy.sum(numpy.log(dy/len(H)))
    return dy, N, l


def likelihood2(H, bias, Ho=0, Hb=7,  alpha_SI=0.67, beta_SI=0.65):
    """

    :param H: Vector of measured H values
    :param bias: Vector of bias on detection of a given object


    :return:
    """
    x = numpy.linspace(H.min()-3.0, H.max()+3.0)
    dx = x[1]-x[0]
    cdf = variably_tapered2(x, Ho, Hb, alpha_SI, beta_SI)
    pdf = numpy.diff(cdf)/dx
    ht = numpy.linspace(H.min(), H.max())
    dh = ht[1]-ht[0]
    dy = bias*numpy.interp(H, x[1:]-dx/2.0, pdf)
    dN = numpy.interp(ht, H, dy)

    N = dh*dN.sum()
    l = -N + numpy.sum(numpy.log(dy/len(H)))
    return dy, N, l


def fraser(H, a1, a2, Ho, Hb, dh):
    """
    LF from Fraser et al. 2014
    :param H: h magnitude
    :param a1: log slope bright of break)
    :param a2: log slope faint of break)
    :param Ho: normalization of bright component
    :param Hb: normalization of faint component
    :return: N
    """
    # print(f"alpha1: {a1}, alpha2: {a2}, Ho: {Ho}, Hb:{Hb}")
    # N[H > Hb] = a2*ln10*10**(a2*(H[H > Hb]-Ho)+(a1-a2)*(Hb-Ho))
    # 10**(a2*H-a2*Ho+a1*Hb-a1*Ho-a2*Hb+a2*Ho)
    # 10**(a2*(H-Hb)+a1*(Hb-Ho))
    # 10**(a1*(Hb-Ho)) * 10**(a2*(H-Hb))
    N = dh*a1*ln10*10**(a1*(H-Ho))
    C = 10**(a1*(Hb-Ho))
    N[H >= Hb] = dh*a1*ln10*C*10**(a2*(H[H > Hb]-Hb))
    return N


def double_plaw(H, simga_23=0.68, a1=1.36, a2=0.38, R_eq=22.8, dh=0.1):
    """
    Double powerlaw form from B14
    :param H:
    :param simga_23:
    :param a1:
    :param a2:
    :return: simga(H)
    """
    r = H + 10*numpy.log10(42.0)
    C = 10**((a2-a2)*(R_eq - 23))
    return dh*((1+C)*simga_23/(10**(-a1*(r-23)) + C*10**(-a2*(r-23))))


def rolling_plaw(H, sigma_23, a1, a2):
    """
    Rolling power law form form Bernstein et al. 2004
    :param H: numpy array of H mags.
    :param sigma_23: normalization at R=23
    :param a1: bright slope.
    :param a2: faint slope.
    :return: N
    """
    r = H + 10*numpy.log10(44)
    return sigma_23 * 10**(a1*(r-23)+a2*(r-23)**2)


def ll(params, H, bias, limits):
    """
    Compute the log-likelihood [setup to work with eemcc
    :param params: function parameters
    :param H: numpy array of H magnitudes
    :param bias: numpy array of detection bias
    :param limits: bounds on params, ll = -np.inf outside this range.
    :return: log of likelihood
    """
    for idx in ['Ho', 'Hb', 'alpha_SI', 'beta_SI']:
        if not limits[idx][0] < params[idx] < limits[idx][1]:
            return -numpy.inf
    return likelihood2(H, bias, alpha_SI=params[2], Hb=params[1], Ho=params[0], beta_SI=params[3])[2]


def poisson_range(k, prob):
    """
    Given an measured count rate, k, and an expectation of a poisson process return estimates of 'mu' that
    are consistent with 'k' being measured within probability prob.

    i.e. if we measure k objects compute a bunch of mu and return the range of mu that are consistent within prob.
    """
    # mu holds the  range of plausible estimates for actual mu
    mu = numpy.arange(max(0, k.min()-20*(k.min())**0.5), k.max()+20*(k.max())**0.5+10, 10+20*k.max()**0.5/100.)
    upper = stats.poisson.ppf(0.5-prob/2., mu)
    lower = stats.poisson.ppf(0.5+prob/2., mu)
    return numpy.interp(k, lower, mu), numpy.interp(k, upper, mu)


def implanted_pdf(H: (numpy.array, float), alpha_dwarf: float = 0.13,
                  H_dwarf: float = 3.2,
                  alpha_trans: float = 0.6,
                  H_trans: float = 6.0, 
                  H_divot: float = 9.0,
                  contrast: float = 0.7,
                  H_elbow:float = 16.5, alpha_elbow=0.15) -> numpy.array:
    """
    Describes the differential size frequency distribution taken from Kavelaars et al. (2021) and Petit et al. (2023)
    and Petit et al. (in prep)
    LF is a power law from -3 to H_dwarf, with slope alpha_dwarf then a powerlaw from H_dwarf
    H_trans and then  a tapered power law from H_trans to H_elbow (following parameters from Kavelaars et al. 2021) and
    then a power law with slope alpha_elbow from H_elbow to the end of the distribution.

    This is a cumulative distribution function, not a differential. The differential is the derivative of this function.

    :param H: numpy array of H magnitudes
    :param alpha_dwarf: slope of the power law for the dwarf planet region.
    :param H_dwarf: break point between the dwarf planet and transition region.
    :param alpha_trans: slope of the power law for the transition region.
    :param H_trans: break point between the transition region and the tapered region.
    :param alpha_elbow: slope of the power law for the small object region.
    :param H_elbow: break point between the tapered region and the small object power law region.
    """
    # Make sure H can be treated as an order numpy array.
    dN = H * 0.0
    H_norm = -3.0  # Normalization of the SFD
    # limit where we are considering dwarf planets. (Petit et al.)
    cond_dwarf = H < H_dwarf
    # transition region between dwarf and tapered SFD (Petit et al.)
    cond_trans = (~cond_dwarf) & (H < H_trans)
    # tapered region (Petit et al.)
    cond_taper = (~cond_dwarf) & (~cond_trans) & (H < H_elbow)
    # break at the 'elbow' (Singer et al.)
    cond_elbow = (~cond_dwarf) & (~cond_trans) & (~cond_taper)

    dN0 = 1*(H[1]-H[0])
    # d/dH(N×10^(a (H - X))) = a N log(10)×10^(a (H - X))
    dN[cond_dwarf] = dN0*alpha_dwarf*ln10*10**(alpha_dwarf*(H[cond_dwarf] - H_norm))
    dN0 = dN[cond_dwarf][-1]/(alpha_trans*ln10)
    # dN0 *= alpha_dwarf*ln10*10**(alpha_dwarf*(H_dwarf - H_norm))
    dN[cond_trans] = dN0*alpha_trans*ln10*10**(alpha_trans*(H[cond_trans] - H_dwarf))
    dN0 = dN[cond_trans][-1]/variably_tapered2_diff(H_trans, Hb = K_Hb) # +0.5)
    # dN0 *= alpha_trans*ln10*10**(alpha_trans*(H_trans - H_dwarf))/variably_tapered2_diff(H_trans)
    dN[cond_taper] = dN0*variably_tapered2_diff(H[cond_taper], Hb = K_Hb) # +0.5)  #, K_Ho, K_Hb, K_ALPHA_SI, K_BETA_SI)
    dN[H>H_divot] *= contrast
    dN0 = dN[cond_taper][-1]/(alpha_elbow*ln10)
    # dN0 *= variably_tapered2_diff(H_elbow) #, K_Ho, K_Hb, K_ALPHA_SI, K_BETA_SI)
    dN[cond_elbow] = dN0 * alpha_elbow*ln10*10**(alpha_elbow*(H[cond_elbow]-H_elbow))
    return dN


def implanted_cfd(H: (numpy.array, float), alpha_dwarf: float = 0.13,
                  H_dwarf: float = 3.2,
                  alpha_trans: float = 0.6,
                  H_trans: float = 6.2, H_elbow:float = 16.5, alpha_elbow=0.15) -> numpy.array:
    """
    Describes the cumulative size frequency distribution taken from Kavelaars et al. (2021) and Petit et al. (2023)
    and Petit et al. (in prep)
    LF is a power law from -3 to H_dwarf, with slope alpha_dwarf then a powerlaw from H_dwarf to H_trans with a slope of
    H_trans and then  a tapered power law from H_trans to H_elbow (following parameters from Kavelaars et al. 2021) and
    then a power law with slope alpha_elbow from H_elbow to the end of the distribution.

    This is a cumulative distribution function, not a differential. The differential is the derivative of this function.

    :param H: numpy array of H magnitudes
    :param alpha_dwarf: slope of the power law for the dwarf planet region.
    :param H_dwarf: break point between the dwarf planet and transition region.
    :param alpha_trans: slope of the power law for the transition region.
    :param H_trans: break point between the transition region and the tapered region.
    :param alpha_elbow: slope of the power law for the small object region.
    :param H_elbow: break point between the tapered region and the small object power law region.
    """
    # Make sure H can be treated as an order numpy array.
    N = H * 0.0
    H_norm = -3.0  # Normalization of the SFD
    # limit where we are considering dwarf planets. (Petit et al.)
    cond_dwarf = H < H_dwarf
    # transition region between dwarf and tapered SFD (Petit et al.)
    cond_trans = (~cond_dwarf) & (H < H_trans)
    # tapered region (Petit et al.)
    cond_taper = (~cond_dwarf) & (~cond_trans) & (H < H_elbow)
    # break at the 'elbow' (Singer et al.)
    cond_elbow = (~cond_dwarf) & (~cond_trans) & (~cond_taper)

    N0 = 1
    N[cond_dwarf] = N0*10**(alpha_dwarf*(H[cond_dwarf] - H_norm))
    N0 *= 10**(alpha_dwarf*(H_dwarf - H_norm))
    N[cond_trans] = N0*10**(alpha_trans*(H[cond_trans] - H_dwarf))
    N0 *= 10**(alpha_trans*(H_trans - H_dwarf))/cold_cfd(H_trans)
    N[cond_taper] = N0*cold_cfd(H[cond_taper])
    N0 *= cold_cfd(H_elbow)
    N[cond_elbow] = N0 * 10**(alpha_elbow*(H[cond_elbow]-H_elbow))

    return N


def cold_cfd(H:numpy.array, Hb=K_Hb) -> numpy.array:
    N = variably_tapered2(H, K_Ho, Hb, K_ALPHA_SI, K_BETA_SI)
    return N


def main():
    logging.info("This is the result from the OSSOS SFD papers.")
    from matplotlib import pyplot as plt
    Ho = -2.6
    Hb = 8.1
    beta_SI = 0.42
    alpha_SI = 0.677
    x = numpy.arange(-10, 20, 0.01)
    # cdf = variably_tapered2(x, Ho, Hb, alpha_SI, beta_SI)
    # pdf = variably_tapered2_diff(x, Ho, Hb, alpha_SI, beta_SI)
    cdf = implanted_cfd(x)
    pdf = implanted_pdf(x, H_dwarf=2.5, H_trans=5.8, alpha_trans=0.7, H_elbow=15.5)  # , H_dwarf=2)
    cold = 2.2*cold_cfd(x)
    napier = 6500*variably_tapered2(x, 2.22, 10.62, 0.26/0.6, 0.19/0.6)
    plt.plot(x, cold, label="K21")
    plt.plot(x, cdf, label="Hot; P23")
    # plt.plot(x, pdf)
    plt.plot(x, napier, label="N24")
    plt.plot(x, pdf.cumsum(), label="Hot; smooth N")
    plt.plot(x, pdf, label="Hot; smooth dN")
    plt.xlabel('H (r)')
    plt.ylabel('N(H<r)')
    plt.ylim(0.1,1E7)
    plt.xlim(-1, 11)
    plt.yscale('log')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.savefig('SFD.pdf')

if __name__ == '__main__':
    main()
