# The purpose of this case Brug220_03

The case concerns a cylindrical building pit of radius R within at t=0 a sudden head change occurs to be fixed thereafter. The idea is to compute the flow into the surrounding infinite aquifer around it over time.

The analytical solution is given by Bruggman (1999) as solution 220.03. But it can also be computed apprximately by convolution, by TTIM, or numerically by Modflow of my own fdm3t.

I have not been able to implement the analytical solution of Bruggeman so that it comes close to the flow computed by TTIM which is close to the approximate convolution method.

However, both Modflow and fdm3t yield the same results which are close to the convolution method.

I stil have to add the TTIM to this example.

But I'm at least convinced that fdm3t and Modflow are correct as is the convolution method.

@TO 20241218
