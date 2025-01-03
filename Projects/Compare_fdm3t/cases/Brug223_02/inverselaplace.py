# Copied moduele from mpmath
# The idea is to do away with arbitrary precision to speed up the inversions.
# contributed to mpmath by Kristopher L. Kuhlman, February 2017
# contributed to mpmath by Guillermo Navas-Palencia, February 2022
from mpmath.calculus import inverselaplace
import numpy as np
from scipy.special import factorial as fac

class InverseLaplaceTransform(object):
    r"""
    Inverse Laplace transform methods are implemented using this
    class, in order to simplify the code and provide a common
    infrastructure.

    Implement a custom inverse Laplace transform algorithm by
    subclassing :class:`InverseLaplaceTransform` and implementing the
    appropriate methods. The subclass can then be used by
    :func:`~mpmath.invertlaplace` by passing it as the *method*
    argument.
    """

    def __init__(self, method=None, **kwargs):
        self.method = method
        self.degree = kwargs.get("degree", 12)

    def calc_time_domain_solution(self, fp):
        r"""
        Compute the time domain solution, after computing the
        Laplace-space function evaluations at the abscissa required
        for the algorithm. Abscissa computed for one algorithm are
        typically not useful for another algorithm.
        
        First determine the vector of Laplace parameter values needed for an
        algorithm, this will depend on the choice of algorithm (de
        Hoog is default), the algorithm-specific parameters passed (or
        default ones), and desired time.
        """
        raise NotImplementedError


class FixedTalbot(InverseLaplaceTransform):

    def calc_time_domain_solution(self, fp, times, **kwargs):
        r"""The "fixed" Talbot method deforms the Bromwich contour towards
        `-\infty` in the shape of a parabola. Traditionally the Talbot
        algorithm has adjustable parameters, but the "fixed" version
        does not. The `r` parameter could be passed in as a parameter,
        if you want to override the default given by (Abate & Valko,
        2004).

        The Laplace parameter is sampled along a parabola opening
        along the negative imaginary axis, with the base of the
        parabola along the real axis at
        `p=\frac{r}{t_\mathrm{max}}`. As the number of terms used in
        the approximation (degree) grows, the abscissa required for
        function evaluation tend towards `-\infty`, requiring high
        precision to prevent overflow.  If any poles, branch cuts or
        other singularities exist such that the deformed Bromwich
        contour lies to the left of the singularity, the method will
        fail.

        **Optional arguments**

        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        recognizes the following keywords

        *tmax*
            maximum time associated with vector of times
            (typically just the time requested)
        *degree*
            integer order of approximation (M = number of terms)
        *r*
            abscissa for `p_0` (otherwise computed using rule
            of thumb `2M/5`)

        The working precision will be increased according to a rule of
        thumb. If 'degree' is not specified, the working precision and
        degree are chosen to hopefully achieve the dps of the calling
        context. If 'degree' is specified, the working precision is
        chosen to achieve maximum resulting precision for the
        specified degree.

        .. math ::

            p_0=\frac{r}{t}

        .. math ::

            p_i=\frac{i r \pi}{Mt_\mathrm{max}}\left[\cot\left(
            \frac{i\pi}{M}\right) + j \right] \qquad 1\le i <M

        where `j=\sqrt{-1}`, `r=2M/5`, and `t_\mathrm{max}` is the
        maximum specified time.

        The fixed Talbot time-domain solution is computed from the
        Laplace-space function evaluations using

        .. math ::

            f(t,M)=\frac{2}{5t}\sum_{k=0}^{M-1}\Re \left[
            \gamma_k \bar{f}(p_k)\right]

        where

        .. math ::

            \gamma_0 = \frac{1}{2}e^{r}\bar{f}(p_0)

        .. math ::

            \gamma_k = e^{tp_k}\left\lbrace 1 + \frac{jk\pi}{M}\left[1 +
            \cot \left( \frac{k \pi}{M} \right)^2 \right] - j\cot\left(
            \frac{k \pi}{M}\right)\right \rbrace \qquad 1\le k<M.

        Again, `j=\sqrt{-1}`.

        Before calling this function, call
        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        to set the parameters and compute the required coefficients.

        **References**

        1. Abate, J., P. Valko (2004). Multi-precision Laplace
           transform inversion. *International Journal for Numerical
           Methods in Engineering* 60:979-993,
           http://dx.doi.org/10.1002/nme.995
        2. Talbot, A. (1979). The accurate numerical inversion of
           Laplace transforms. *IMA Journal of Applied Mathematics*
           23(1):97, http://dx.doi.org/10.1093/imamat/23.1.97
        """
        if np.isscalar(times):
            times = np.array([times])
        tmax = times[-1]
        
        M = kwargs.get('degree', 12)

        # Abate & Valko rule of thumb for r parameter
        self.r = kwargs.get('r', 2 / 5 * M)

        theta = np.arange(M, dtype=np.complex128) * np.pi / M
        cot = np.zeros_like(theta)
        cot[1:] = np.tan(1 / theta[1:])
        
        # All but time-dependent part of p
        # This is p_i in the formula
        gamma  = np.zeros_like(theta)
        p      = np.zeros_like(theta)
        FP     = np.zeros_like(theta)
        results = np.zeros(len(times), dtype=np.complex128)
                     
        p[1:] = self.r / tmax * theta[1:] * (cot[1:] + 1j)        
        FP[1:]= np.array([fp(p_) for p_ in p[1:]], dtype=np.complex128)
        
        for it, t in enumerate(times):
            p[0] = self.r / t
            FP[0] = fp(p[0])
            
            gamma[0] = 0.5 * np.exp(self.r) * FP[0]
            gamma[1:] =np.exp(t * p[1:]) * (
                1 + 1j * theta[1:] * (1 + cot[1:] ** 2) - 1j * cot[1:]) 

            results[it] = 2 / (5 * t) * np.dot(gamma, FP)

        return results.real

# ****************************************

class Stehfest(InverseLaplaceTransform):
    r"""
    The Gaver-Stehfest method is a discrete approximation of the
    Widder-Post inversion algorithm, rather than a direct
    approximation of the Bromwich contour integral.
    The method abscissa along the real axis, and therefore has
    issues inverting oscillatory functions (which have poles in
    pairs away from the real axis).
    The working precision will be increased according to a rule of
    thumb. If 'degree' is not specified, the working precision and
    degree are chosen to hopefully achieve the dps of the calling
    context. If 'degree' is specified, the working precision is
    chosen to achieve maximum resulting precision for the
    specified degree.
    .. math ::
        p_k = \frac{k \log 2}{t} \qquad 1 \le k \le M
    """
    
    def _coeff(self, M):
        r"""Salzer summation weights (aka, "Stehfest coefficients")
        only depend on the approximation order (M) and the precision.
        
        Salzer summation weights
        get very large in magnitude and oscillate in sign,
        if the precision is not high enough, there will be
        catastrophic cancellation.
        """
        if M % 2 != 0:
            M = int(M) + 1
        M2 = M // 2  # checked earlier that M is even
        
        V = np.zeros(M)

        for k in range(1, M + 1):
            z = np.zeros(min(k, M2) + 1)
            for j in range((k + 1) // 2, min(k, M2) + 1):
                z[j] = (j ** M2 * fac(2 * j) / (
                    fac(M2 - j) * fac(j) * fac(j - 1) * fac(k - j) * fac(2*j - k)))
            V[k-1] = (-1) ** (k + M2) * np.sum(z)

        return V

    def calc_time_domain_solution(self, fp, times, **kwargs):
        r"""Compute time-domain Stehfest algorithm solution.

        .. math ::

            f(t,M) = \frac{\log 2}{t} \sum_{k=1}^{M} V_k \bar{f}\left(
            p_k \right)

        where

        .. math ::

            V_k = (-1)^{k + N/2} \sum^{\min(k,N/2)}_{i=\lfloor(k+1)/2 \rfloor}
            \frac{i^{\frac{N}{2}}(2i)!}{\left(\frac{N}{2}-i \right)! \, i! \,
            \left(i-1 \right)! \, \left(k-i\right)! \, \left(2i-k \right)!}

        As the degree increases, the abscissa (`p_k`) only increase
        linearly towards `\infty`, but the Stehfest coefficients
        (`V_k`) alternate in sign and increase rapidly in sign,
        requiring high precision to prevent overflow or loss of
        significance when evaluating the sum.

        **References**

        1. Widder, D. (1941). *The Laplace Transform*. Princeton.
        2. Stehfest, H. (1970). Algorithm 368: numerical inversion of
           Laplace transforms. *Communications of the ACM* 13(1):47-49,
           http://dx.doi.org/10.1145/361953.361969           
           
        Stef = 
        """
        # required
        if np.isscalar(times):
            times = [times]
            
        M = 2 * (kwargs.get('degree', 10) // 2)

        self.V = self._coeff(M)
        
        result = np.zeros(len(times), dtype=float)

        p_no_t = np.arange(1, M+1) * np.log(2.)
        for it, t in enumerate(times):
            p = p_no_t / t            
            result[it] = np.sum(self.V * np.array([fp(p_) for p_ in p])) * np.log(2.) / t

        return result.real

# ****************************************

class DeHoog(InverseLaplaceTransform):
    r"""the de Hoog, Knight & Stokes algorithm is an
    accelerated form of the Fourier series numerical
    inverse Laplace transform algorithms.
    .. math ::
        p_k = \gamma + \frac{jk}{T} \qquad 0 \le k < 2M+1
    where
    .. math ::
        \gamma = \alpha - \frac{\log \mathrm{tol}}{2T},
    `j=\sqrt{-1}`, `T = 2t_\mathrm{max}` is a scaled time,
    `\alpha=10^{-\mathrm{dps\_goal}}` is the real part of the
    rightmost pole or singularity, which is chosen based on the
    desired accuracy (assuming the rightmost singularity is 0),
    and `\mathrm{tol}=10\alpha` is the desired tolerance, which is
    chosen in relation to `\alpha`.`
    When increasing the degree, the abscissa increase towards
    `j\infty`, but more slowly than the fixed Talbot
    algorithm. The de Hoog et al. algorithm typically does better
    with oscillatory functions of time, and less well-behaved
    functions. The method tends to be slower than the Talbot and
    Stehfest algorithsm, especially so at very high precision
    (e.g., `>500` digits precision).
    print("This method does nothing)
    """

    def calc_time_domain_solution(self, fp, times, **kwargs):
        r"""Calculate time-domain solution for
        de Hoog, Knight & Stokes algorithm.

        The un-accelerated Fourier series approach is:

        .. math ::

            f(t,2M+1) = \frac{e^{\gamma t}}{T} \sum_{k=0}^{2M}{}^{'}
            \Re\left[\bar{f}\left( p_k \right)
            e^{i\pi t/T} \right],

        where the prime on the summation indicates the first term is halved.

        This simplistic approach requires so many function evaluations
        that it is not practical. Non-linear acceleration is
        accomplished via Pade-approximation and an analytic expression
        for the remainder of the continued fraction. See the original
        paper (reference 2 below) a detailed description of the
        numerical approach.

        **References**

        1. Davies, B. (2005). *Integral Transforms and their
           Applications*, Third Edition. Springer.
        2. de Hoog, F., J. Knight, A. Stokes (1982). An improved
           method for numerical inversion of Laplace transforms. *SIAM
           Journal of Scientific and Statistical Computing* 3:357-366,
           http://dx.doi.org/10.1137/0903022

        """
        if np.isscalar(times):
            times = np.array([times])
        
        M = kwargs.get('degree', 10)
        
        alpha = kwargs.get('alpha', 1e-15)
        
        tol = 10 * alpha
        
        NP = 2 * M + 1  # number of terms in approximation

        # scaling factor (likely tunable, but 2 is typical)
        scale = kwargs.get('scale', 2)
        T = kwargs.get('T', scale * times[-1])
        gamma = alpha - np.log(tol) / (scale * T)
        
        p = gamma + 1j / T * np.arange(NP, dtype=np.complex128)
        FP = np.array([fp(p_) for p_ in p], dtype=np.complex128)
                
        # would it be useful to try re-using
        # space between e&q and A&B?
        e = np.zeros((NP, M+1), dtype=np.complex128)
        q = np.zeros((2 * M, M), dtype=np.complex128)
        d = np.zeros((NP, 1),   dtype=np.complex128)
        
        A = np.zeros((NP+1, 1), dtype=np.complex128)
        B = np.ones((NP+1, 1),  dtype=np.complex128)

        # initialize Q-D table
        e[:, 0] = 0.0 + 0j
        q[0, 0] = FP[1] / (FP[0] / 2)
        for i in range(1, 2 * M):
            q[i, 0] = FP[i + 1] / FP[i]

        # rhombus rule for filling triangular Q-D table (e & q)
        for r in range(1, M+1):            
            # start with e, column 1, 0:2*M-2
            mr = 2 * (M - r) + 1
            e[0:mr, r] = q[1:mr + 1, r - 1] - q[0:mr, r - 1] + e[1:mr + 1, r - 1]
            if not r == M:
                rq = r + 1
                mr = 2 * (M - rq) + 1 + 2
                for i in range(mr):
                    q[i, rq - 1] = q[i + 1, rq - 2] * e[i + 1, rq - 1] / e[i, rq - 1]

        # build up continued fraction coefficients (d)
        d[0] = FP[0] / 2
        for r in range(1, M + 1):
            d[2 * r - 1] = -q[0, r - 1]  # even terms
            d[2 * r]     = -e[0, r]      # odd terms

        # seed A and B for recurrence
        A[0] = 0.0 + 0.0j
        A[1] = d[0]
        B[0:2] = 1.0 + 0.0j

        results = np.zeros(len(times), dtype=np.float64)

        for it, t in enumerate(times):
            # base of the power series
            z = np.exp(1j * np.pi * t / T)

            # coefficients of Pade approximation (A & B)
            # using recurrence for all but last term
            for i in range(1, 2 * M):
                A[i + 1] = A[i] + d[i] * A[i - 1] * z
                B[i + 1] = B[i] + d[i] * B[i - 1] * z

            # "improved remainder" to continued fraction
            brem = (1 + (d[2 * M - 1] - d[2 * M]) * z) / 2
            # powm1(x,y) computes x^y - 1 more accurately near zero
            # rem = brem * self.ctx.powm1(1 + d[2 * M] * z / brem, 1 / 2)
            rem = brem * (np.sqrt(1 + d[2 * M] * z / brem) - 1)

            # last term of recurrence using new remainder
            A[NP] = A[2 * M] + rem * A[2 * M - 1]
            B[NP] = B[2 * M] + rem * B[2 * M - 1]

            # diagonal Pade approximation
            # F = A / B represents accelerated trapezoid rule
            results[it] = np.exp(gamma * t) / T * (A[NP] / B[NP]).real

        return results

# ****************************************

class Cohen(InverseLaplaceTransform):
    r"""The Cohen algorithm accelerates the convergence of the nearly
    alternating series resulting from the application of the trapezoidal
    rule to the Bromwich contour inversion integral.
    
    .. math ::
        p_k = \frac{\gamma}{2 t} + \frac{\pi i k}{t} \qquad 0 \le k < M
    
    where
    
    .. math ::
        \gamma = \frac{2}{3} (d + \log(10) + \log(2 t)),
    
    `d = \mathrm{dps\_goal}`, which is chosen based on the desired
    accuracy using the method developed in [1] to improve numerical
    stability. The Cohen algorithm shows robustness similar to the de Hoog
    et al. algorithm, but it is faster than the fixed Talbot algorithm.
    
    **Optional arguments**
    
    *degree*
        integer order of the approximation (M = number of terms)
    
    *alpha*
        abscissa for `p_0` (controls the discretization error)
    The working precision will be increased according to a rule of
    thumb. If 'degree' is not specified, the working precision and
    degree are chosen to hopefully achieve the dps of the calling
    context. If 'degree' is specified, the working precision is
    chosen to achieve maximum resulting precision for the
    specified degree.
    
    **References**
    1. P. Glasserman, J. Ruiz-Mata (2006). Computing the credit loss
    distribution in the Gaussian copula model: a comparison of methods.
    *Journal of Credit Risk* 2(4):33-66, 10.21314/JCR.2006.057
    """

    def calc_time_domain_solution(self, fp, times, **kwargs):
        r"""Calculate time-domain solution for Cohen algorithm.

        The accelerated nearly alternating series is:

        .. math ::

            f(t, M) = \frac{e^{\gamma / 2}}{t} \left[\frac{1}{2}
            \Re\left(\bar{f}\left(\frac{\gamma}{2t}\right) \right) -
            \sum_{k=0}^{M-1}\frac{c_{M,k}}{d_M}\Re\left(\bar{f}
            \left(\frac{\gamma + 2(k+1) \pi i}{2t}\right)\right)\right],

        where coefficients `\frac{c_{M, k}}{d_M}` are described in [1].

        1. H. Cohen, F. Rodriguez Villegas, D. Zagier (2000). Convergence
        acceleration of alternating series. *Experiment. Math* 9(1):3-12

        """
        if np.isscalar(times):
            times = [times]

        N = kwargs.get('degree', 22)
        M = N + 1

        def gamma(t):
            dps = 15
            return 2 / 3 * (dps * np.log(10) + np.log(2 * t))        
        
        A = np.zeros(M)
        d = (3 + np.sqrt(8)) ** N
        d = (d + 1 / d) / 2
        b = -1
        c = -d
        s = 0

        results = np.zeros(len(times))

        for it, t in enumerate(times):            
            p = gamma(t) / (2 * t) + 1j * np.pi * np.arange(
                M, dtype=np.complex128) / t    
            
            for m in range(M):
                A[m] = fp(p[m]).real

            for k in range(N):
                c = b - c
                s = s + c * A[k + 1]
                b = 2 * (k + N) * (k - N) * b / ((2 * k + 1) * (k + 1))

            results[it] = np.exp(gamma(t) / 2) / t * (A[0] / 2 - s / d)

        return results

# ****************************************

class LaplaceTransformInversionMethods(object):
    def __init__(ctx, *args, **kwargs):
        ctx._fixed_talbot = FixedTalbot(ctx)
        ctx._stehfest = Stehfest(ctx)
        ctx._de_hoog = deHoog(ctx)
        ctx._cohen = Cohen(ctx)

    def invertlaplace(ctx, f, t, **kwargs):
        r"""Computes the numerical inverse Laplace transform for a
        Laplace-space function at a given time.  The function being
        evaluated is assumed to be a real-valued function of time.

        The user must supply a Laplace-space function `\bar{f}(p)`,
        and a desired time at which to estimate the time-domain
        solution `f(t)`.

        A few basic examples of Laplace-space functions with known
        inverses (see references [1,2]) :

        .. math ::

            \mathcal{L}\left\lbrace f(t) \right\rbrace=\bar{f}(p)

        .. math ::

            \mathcal{L}^{-1}\left\lbrace \bar{f}(p) \right\rbrace = f(t)

        .. math ::

            \bar{f}(p) = \frac{1}{(p+1)^2}

        .. math ::

            f(t) = t e^{-t}

        >>> from mpmath import *
        >>> mp.dps = 15; mp.pretty = True
        >>> tt = [0.001, 0.01, 0.1, 1, 10]
        >>> fp = lambda p: 1/(p+1)**2
        >>> ft = lambda t: t*exp(-t)
        >>> ft(tt[0]),ft(tt[0])-invertlaplace(fp,tt[0],method='talbot')
        (0.000999000499833375, 8.57923043561212e-20)
        >>> ft(tt[1]),ft(tt[1])-invertlaplace(fp,tt[1],method='talbot')
        (0.00990049833749168, 3.27007646698047e-19)
        >>> ft(tt[2]),ft(tt[2])-invertlaplace(fp,tt[2],method='talbot')
        (0.090483741803596, -1.75215800052168e-18)
        >>> ft(tt[3]),ft(tt[3])-invertlaplace(fp,tt[3],method='talbot')
        (0.367879441171442, 1.2428864009344e-17)
        >>> ft(tt[4]),ft(tt[4])-invertlaplace(fp,tt[4],method='talbot')
        (0.000453999297624849, 4.04513489306658e-20)

        The methods also work for higher precision:

        >>> mp.dps = 100; mp.pretty = True
        >>> nstr(ft(tt[0]),15),nstr(ft(tt[0])-invertlaplace(fp,tt[0],method='talbot'),15)
        ('0.000999000499833375', '-4.96868310693356e-105')
        >>> nstr(ft(tt[1]),15),nstr(ft(tt[1])-invertlaplace(fp,tt[1],method='talbot'),15)
        ('0.00990049833749168', '1.23032291513122e-104')

        .. math ::

            \bar{f}(p) = \frac{1}{p^2+1}

        .. math ::

            f(t) = \mathrm{J}_0(t)

        >>> mp.dps = 15; mp.pretty = True
        >>> fp = lambda p: 1/sqrt(p*p + 1)
        >>> ft = lambda t: besselj(0,t)
        >>> ft(tt[0]),ft(tt[0])-invertlaplace(fp,tt[0],method='dehoog')
        (0.999999750000016, -6.09717765032273e-18)
        >>> ft(tt[1]),ft(tt[1])-invertlaplace(fp,tt[1],method='dehoog')
        (0.99997500015625, -5.61756281076169e-17)

        .. math ::

            \bar{f}(p) = \frac{\log p}{p}

        .. math ::

            f(t) = -\gamma -\log t

        >>> mp.dps = 15; mp.pretty = True
        >>> fp = lambda p: log(p)/p
        >>> ft = lambda t: -euler-log(t)
        >>> ft(tt[0]),ft(tt[0])-invertlaplace(fp,tt[0],method='stehfest')
        (6.3305396140806, -1.92126634837863e-16)
        >>> ft(tt[1]),ft(tt[1])-invertlaplace(fp,tt[1],method='stehfest')
        (4.02795452108656, -4.81486093200704e-16)

        **Options**

        :func:`~mpmath.invertlaplace` recognizes the following optional
        keywords valid for all methods:

        *method*
            Chooses numerical inverse Laplace transform algorithm
            (described below).
        *degree*
            Number of terms used in the approximation

        **Algorithms**

        Mpmath implements four numerical inverse Laplace transform
        algorithms, attributed to: Talbot, Stehfest, and de Hoog,
        Knight and Stokes. These can be selected by using
        *method='talbot'*, *method='stehfest'*, *method='dehoog'* or
        *method='cohen'* or by passing the classes *method=FixedTalbot*,
        *method=Stehfest*, *method=deHoog*, or *method=Cohen*. The functions
        :func:`~mpmath.invlaptalbot`, :func:`~mpmath.invlapstehfest`,
        :func:`~mpmath.invlapdehoog`, and :func:`~mpmath.invlapcohen`
        are also available as shortcuts.

        All four algorithms implement a heuristic balance between the
        requested precision and the precision used internally for the
        calculations. This has been tuned for a typical exponentially
        decaying function and precision up to few hundred decimal
        digits.

        The Laplace transform converts the variable time (i.e., along
        a line) into a parameter given by the right half of the
        complex `p`-plane.  Singularities, poles, and branch cuts in
        the complex `p`-plane contain all the information regarding
        the time behavior of the corresponding function. Any numerical
        method must therefore sample `p`-plane "close enough" to the
        singularities to accurately characterize them, while not
        getting too close to have catastrophic cancellation, overflow,
        or underflow issues. Most significantly, if one or more of the
        singularities in the `p`-plane is not on the left side of the
        Bromwich contour, its effects will be left out of the computed
        solution, and the answer will be completely wrong.

        *Talbot*

        The fixed Talbot method is high accuracy and fast, but the
        method can catastrophically fail for certain classes of time-domain
        behavior, including a Heaviside step function for positive
        time (e.g., `H(t-2)`), or some oscillatory behaviors. The
        Talbot method usually has adjustable parameters, but the
        "fixed" variety implemented here does not. This method
        deforms the Bromwich integral contour in the shape of a
        parabola towards `-\infty`, which leads to problems
        when the solution has a decaying exponential in it (e.g., a
        Heaviside step function is equivalent to multiplying by a
        decaying exponential in Laplace space).

        *Stehfest*

        The Stehfest algorithm only uses abscissa along the real axis
        of the complex `p`-plane to estimate the time-domain
        function. Oscillatory time-domain functions have poles away
        from the real axis, so this method does not work well with
        oscillatory functions, especially high-frequency ones. This
        method also depends on summation of terms in a series that
        grows very large, and will have catastrophic cancellation
        during summation if the working precision is too low.

        *de Hoog et al.*

        The de Hoog, Knight, and Stokes method is essentially a
        Fourier-series quadrature-type approximation to the Bromwich
        contour integral, with non-linear series acceleration and an
        analytical expression for the remainder term. This method is
        typically one of the most robust. This method also involves the
        greatest amount of overhead, so it is typically the slowest of the
        four methods at high precision.

        *Cohen*

        The Cohen method is a trapezoidal rule approximation to the Bromwich
        contour integral, with linear acceleration for alternating
        series. This method is as robust as the de Hoog et al method and the
        fastest of the four methods at high precision, and is therefore the
        default method.

        **Singularities**

        All numerical inverse Laplace transform methods have problems
        at large time when the Laplace-space function has poles,
        singularities, or branch cuts to the right of the origin in
        the complex plane. For simple poles in `\bar{f}(p)` at the
        `p`-plane origin, the time function is constant in time (e.g.,
        `\mathcal{L}\left\lbrace 1 \right\rbrace=1/p` has a pole at
        `p=0`). A pole in `\bar{f}(p)` to the left of the origin is a
        decreasing function of time (e.g., `\mathcal{L}\left\lbrace
        e^{-t/2} \right\rbrace=1/(p+1/2)` has a pole at `p=-1/2`), and
        a pole to the right of the origin leads to an increasing
        function in time (e.g., `\mathcal{L}\left\lbrace t e^{t/4}
        \right\rbrace = 1/(p-1/4)^2` has a pole at `p=1/4`).  When
        singularities occur off the real `p` axis, the time-domain
        function is oscillatory. For example `\mathcal{L}\left\lbrace
        \mathrm{J}_0(t) \right\rbrace=1/\sqrt{p^2+1}` has a branch cut
        starting at `p=j=\sqrt{-1}` and is a decaying oscillatory
        function, This range of behaviors is illustrated in Duffy [3]
        Figure 4.10.4, p. 228.

        In general as `p \rightarrow \infty` `t \rightarrow 0` and
        vice-versa. All numerical inverse Laplace transform methods
        require their abscissa to shift closer to the origin for
        larger times. If the abscissa shift left of the rightmost
        singularity in the Laplace domain, the answer will be
        completely wrong (the effect of singularities to the right of
        the Bromwich contour are not included in the results).

        For example, the following exponentially growing function has
        a pole at `p=3`:

        .. math ::

            \bar{f}(p)=\frac{1}{p^2-9}

        .. math ::

            f(t)=\frac{1}{3}\sinh 3t

        >>> mp.dps = 15; mp.pretty = True
        >>> fp = lambda p: 1/(p*p-9)
        >>> ft = lambda t: sinh(3*t)/3
        >>> tt = [0.01,0.1,1.0,10.0]
        >>> ft(tt[0]),invertlaplace(fp,tt[0],method='talbot')
        (0.0100015000675014, 0.0100015000675014)
        >>> ft(tt[1]),invertlaplace(fp,tt[1],method='talbot')
        (0.101506764482381, 0.101506764482381)
        >>> ft(tt[2]),invertlaplace(fp,tt[2],method='talbot')
        (3.33929164246997, 3.33929164246997)
        >>> ft(tt[3]),invertlaplace(fp,tt[3],method='talbot')
        (1781079096920.74, -1.61331069624091e-14)

        **References**

        1. [DLMF]_ section 1.14 (http://dlmf.nist.gov/1.14T4)
        2. Cohen, A.M. (2007). Numerical Methods for Laplace Transform
           Inversion, Springer.
        3. Duffy, D.G. (1998). Advanced Engineering Mathematics, CRC Press.

        **Numerical Inverse Laplace Transform Reviews**

        1. Bellman, R., R.E. Kalaba, J.A. Lockett (1966). *Numerical
           inversion of the Laplace transform: Applications to Biology,
           Economics, Engineering, and Physics*. Elsevier.
        2. Davies, B., B. Martin (1979). Numerical inversion of the
           Laplace transform: a survey and comparison of methods. *Journal
           of Computational Physics* 33:1-32,
           http://dx.doi.org/10.1016/0021-9991(79)90025-1
        3. Duffy, D.G. (1993). On the numerical inversion of Laplace
           transforms: Comparison of three new methods on characteristic
           problems from applications. *ACM Transactions on Mathematical
           Software* 19(3):333-359, http://dx.doi.org/10.1145/155743.155788
        4. Kuhlman, K.L., (2013). Review of Inverse Laplace Transform
           Algorithms for Laplace-Space Numerical Approaches, *Numerical
           Algorithms*, 63(2):339-355.
           http://dx.doi.org/10.1007/s11075-012-9625-3

        """

        rule = kwargs.get('method', 'cohen')
        if type(rule) is str:
            lrule = rule.lower()
            if lrule == 'talbot':
                rule = ctx._fixed_talbot
            elif lrule == 'stehfest':
                rule = ctx._stehfest
            elif lrule == 'dehoog':
                rule = ctx._de_hoog
            elif rule == 'cohen':
                rule = ctx._cohen
            else:
                raise ValueError("unknown invlap algorithm: %s" % rule)
        else:
            rule = rule(ctx)

        # determine the vector of Laplace-space parameter
        # needed for the requested method and desired time
        rule.calc_laplace_parameter(t, **kwargs)

        # compute the Laplace-space function evalutations
        # at the required abscissa.
        fp = [f(p) for p in rule.p]

        # compute the time-domain solution from the
        # Laplace-space function evaluations
        return rule.calc_time_domain_solution(fp, t)

    # shortcuts for the above function for specific methods
    def invlaptalbot(ctx, *args, **kwargs):
        kwargs['method'] = 'talbot'
        return ctx.invertlaplace(*args, **kwargs)

    def invlapstehfest(ctx, *args, **kwargs):
        kwargs['method'] = 'stehfest'
        return ctx.invertlaplace(*args, **kwargs)

    def invlapdehoog(ctx,  *args, **kwargs):
        kwargs['method'] = 'dehoog'
        return ctx.invertlaplace(*args, **kwargs)

    def invlapcohen(ctx, *args, **kwargs):
        kwargs['method'] = 'cohen'
        return ctx.invertlaplace(*args, **kwargs)


# ****************************************

if __name__ == '__main__':
    import doctest
    doctest.testmod()
