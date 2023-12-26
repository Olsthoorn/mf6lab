#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:10:00 2020

@author: Theo
"""
import matplotlib.pyplot as plt
import numpy as np
import pdb

def newfig(title='title?', xlabel='xlabel?', ylabel='ylabel?', xscale='linear', yscale='linear',
           xlim=None, ylim=None, size_inches=(12, 11), **kwargs):
    """Return ax for new plot."""
    fig, ax = plt.subplots()
    fig.set_size_inches(size_inches)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.grid()
    return ax


def newfig2(titles=['title1', 'title2'], xlabel='time',
            ylabels=['heads [m]', 'flows [m2/d]'],
            xscale='linear',
            yscales=['linear', 'linear'],
            sharex=True,
            sharey=False,
            xlims=None,
            ylims=None,
            size_inches=(12, 6),
            **kwargs):
    """Return ax[0], ax[1] for new plot."""
    fig, ax = plt.subplots(2, 1, sharex=sharex, sharey=sharey)
    fig.set_size_inches(size_inches)
    for a, title, ylabel, yscale in zip(ax, titles, ylabels, yscales):
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.set_xscale(xscale)
        a.set_yscale(yscale)
        a.grid()
    if xlims is not None:
        ax[0].set_xlim(xlims[0])
        ax[1].set_xlim(xlims[1])
    if ylims is not None:
        ax[0].set_ylim(ylims[0])
        ax[1].set_ylim(ylims[1])
    return ax


def get_cdr(h0=None, zdr=None, cdrain=None, htol=None, rising=None, **kwargs):
    if h0 > zdr + htol:
        cdr = np.inf
    elif h0 < zdr - htol:
        cdr = cdrain
    else:
        if rising:
            cdr = np.inf
        else:
            cdr = cdrain
    return cdr


def headt(ax=None, t=None, t0=None, phih=None, h0=None, Th=None, B=None, **kwargs):
    e = np.exp(-(t - t0) / Th)
    htau = kw['phih'] +  (kw['h0'] - kw['phih']) * e + kw['B'] * (1 -e)
    return htau

def plot_all(ax=None, t=None, ht=None, h0=None, phi=None, phih=None,
             tauhit=None, htauhit=None, rising=None, **kwargs):
    ax.hlines(zdr, t[0], t[-1], label='zdr', lw = 0.5, color='green')
    ax.hlines(h0,  t[0], t[-1], label='h0' , lw=0.5, color='blue')
    ax.plot(t, ht, '.-', label='head', lw=2, color='black')
    if tauhit:
        ax.hlines(phih, t[0], t[ 0] + tauhit, label='phih', lw=2, color='red')
        ax.hlines(phi,  t[0], t[ 0] + tauhit, label='phi' , lw=0.5, color='purple')
        ax.plot(t[0] + tauhit, htauhit, 'ro', label='hit')
    print('tauhit = {}, rising = {}'.format(tauhit, rising))


# Because cdr is inf for h below zdr, the way ch and phih are calculated matter!
def update_values(cdrain=None, c=None, zdr=None, h0=None, phi=None, hLR=None,
           N=None, kD=None, **kwargs):
    """Update values in kwargs to new values in input, in place."""
    cdr = get_cdr(h0=h0, cdrain=cdrain, zdr=zdr, htol=0.001, **kwargs)
    ch  = c / (c / cdr + 1)
    phih = (phi  + zdr * c / cdr) / (c / cdr + 1)
    lamh = np.sqrt(kD * ch)
    Th   = mu * ch
    Lamh = (1 /((b / lamh) *np.cosh(b / lamh)
            / np.sinh(b / lamh) + (b / ch)/ (D / w)))
    B = N * ch - (N * ch - (hLR - phih)) * Lamh

    logarg = (h0 - phih - B) / (zdr - phih - B)
    if logarg > 0:
        tauhit =  Th * np.log(logarg)
        rising = (phih - h0 + B) > 0
    else:
        tauhit = np.nan
        rising = np.nan

    return {'cdr':cdr, 'ch': ch, 'phih': phih, 'Th': Th, 'lamh': lamh,
          'Lamh':Lamh, 'B':B, 'tauhit': tauhit, 'rising':rising}

#%%
if __name__ == '__main__':

    # parameters
    N = 0.002
    h0, phi, hLR, zdr,  = 0., 0.5, 0.5, 1.0
    b, mu, k, D = 75., 0.15, 10., 10.,
    c, cdrain, w = 150, 5, 0.5
    tend = 40

    # Initialize kw dict
    kw = {'N':N, 'h0':h0, 'phi':phi, 'hLR': hLR, 'zdr':zdr,
             'b':b, 'mu':mu, 'k':k, 'D':D, 'c':c, 'cdrain':cdrain, 'kD': k * D}

    ax = newfig('', xlabel='time', ylabel='m', size_inches=(14, 8))

    kw['ax'] = ax
    pdb.set_trace()
    kw.update(update_values(**kw))

    title=('tauhit={:.4g}, rising ={:.4g}, kD={}, c={}, cdr={}, ch={:.4g}'.format(
            kw['tauhit'], kw['rising'], kw['kD'], kw['c'], kw['cdr'], kw['ch']) +
            ', phi={}, phih={:.4g}, N={}'.format(kw['phi'], kw['phih'], kw['N']))

    t  = np.linspace(0, tend, 30)
    ht = headt(t=t, t0=t[0], **kw)
    hthit = headt(t=t[0] + kw['tauhit'], t0=t[0], **kw)
    plot_all(t=t, ht=ht, htauhit=hthit, **kw)

    thit = t[0] + kw['tauhit']
    h0   = headt(tau=t[0] + kw['tauhit'],t0=t[0], **kw)

    #  Second part of time step, after hitting zdr
    t = thit + t

    cdr = np.inf if h0 < zdr else cdrain
    plot_all(ax, t, cdr, h0)


    ax.set_title(title)
    ax.set_xlabel('tau [d]')
    ax.set_ylabel('elevation and head')
    ax.legend()
    plt.show()
