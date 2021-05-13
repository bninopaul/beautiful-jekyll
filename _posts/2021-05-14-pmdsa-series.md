---
layout: post
title: Stat217 Machine Problem 1
image: /img/hello_world.jpeg
tags: [pmdsa]
---

## Problem 1

(Inference) The following data are an IID sample from a $Cauchy (\theta, 1)$ distribution: 1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21.

a. Graph the log likelihood function.

b. Find the MLE for $\theta$ using the Halley’s method. Try the following starting points: -11, 0, 1.5, 4, 4.7, 7, 38, and 51. Consider also using the mean and median of the data as starting points.Discuss the results.

## Solution to Problem 1

To answer item (a), we need to first derive the log likelihood function of $Cauchy (\theta, 1)$ distribution. The probability density function is given by:
$$f(x) = \frac{1}{\pi(1+(x-\theta)^2)}$$
Thus, the likelihood function will be the product of the pdfs evaluated at different values of x:
$$\mathcal{L}(\theta|x) = \prod_{i=1}^{n} \frac{1}{\pi(1+(x_i-\theta)^2)}$$
$$\mathcal{L}(\theta|x) =\frac{1}{\pi^n}\prod_{i=1}^{n} \frac{1}{(1+(x_i-\theta)^2)}$$
Taking the log, we get:
$$\ell(\theta|x) = log \Big( \frac{1}{\pi^n}\prod_{i=1}^{n} \frac{1}{(1+(x_i-\theta)^2)}\Big)$$
$$\vdots$$
$$\ell(\theta|x) = -n log(\pi) - \sum_{i=1}^{n}log(1+(x_i-\theta)^2)$$
Now, we can plot the log likelihood function.

```{r problem1}
samples <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

# Create the function for the computing log likelihood of Cauchy distribution
compute_log_likelihood_cauchy <- function(theta, x){
  n <- length(x)
  value <- - n*log(pi)-sum(log(1+(x-theta)**2)) 
  return(value)
}

# Set different values of theta from -100 to 100
thetas <- seq(from=-100, to=100, by=0.01)

# Evaluate the log likelihoods for each value of theta
log_likelihood_values <- unlist(lapply(thetas, compute_log_likelihood_cauchy, x=samples))

# Plot the log likelihood values
plot(thetas, log_likelihood_values, type="l", col="blue", xlab=TeX("theta($\\theta$)"), ylab="log likelihood value")
```


The maximum value is located around $\theta=0$ however, there are also noticeable local peaks/cusps around $\theta=50$ and $\theta=-15$. To explore more, let's answer item (b), which is to find the MLE using Halley's method.

```{r fig.show="hold", out.width="33%", fig.height=10}
plot(thetas, log_likelihood_values, type="l", col="blue", xlab=TeX("theta($\\theta$)"), ylab="log likelihood value", xlim = c(-25, -5))
plot(thetas, log_likelihood_values, type="l", col="blue", xlab=TeX("theta($\\theta$)"), ylab="log likelihood value", xlim = c(-5, 5))
plot(thetas, log_likelihood_values, type="l", col="blue", xlab=TeX("theta($\\theta$)"), ylab="log likelihood value", xlim = c(30, 70))
```

Halley's method states that we can find the value of $x$(root) of $f(x)=0$ using the update formula: 

$$x_{n+1} = x_{n} - \frac{2f(x_n)f^{'}(x_n)}{2[f^{'}(x_n)]^2 - f(x_n)f^{''}(x_n)}$$

First, let us set up our $f(x) = 0$ expression, which is  $f(\theta)=\frac{\partial \ell}{\partial \theta} = 0$.

Taking the derivative of log likelihood function we get:
$$\ell(\theta|x) = -n log(\pi) - \sum_{i=1}^{n}log(1+(x_i-\theta)^2)$$

$$\frac{\partial\ell}{\partial\theta} = 0 - \sum_{i=1}^{n} \frac{1}{1+(x_i-\theta)^2}(2)(x_i-\theta)(-1)$$
$$\frac{\partial\ell}{\partial\theta} =2\sum_{i=1}^{n} \frac{(x_i - \theta)}{1+(x_i-\theta)^2}$$
Thus, we need to solve:
$$f(\theta) = \frac{\partial\ell}{\partial\theta} =2\sum_{i=1}^{n} \frac{(x_i - \theta)}{1+(x_i-\theta)^2} = 0$$
First, we need to derive $f'(\theta)$ and $f''(\theta)$.

Computing for $f'(\theta)$:
$$f'(\theta) = 2\sum_{i=1}^{n} \frac{\big(1+(x_i-\theta)^2\big)\big(-1\big) - \big(x_i - \theta\big)\big((2)(x_i - \theta)(-1)\big)}{\big(1+(x_i - \theta)^2\big)^2}$$
$$\vdots$$
$$f'(\theta) = 2\sum_{i=1}^{n} \frac{(x_i-\theta)^2 - 1}{\big(1+(x_i - \theta)^2\big)^2}$$
Computing for $f''(\theta)$:

$$f''(\theta) = 2\sum_{i=1}^{n} \frac{\bigg(\big(1+(x_i - \theta)^2\big)^2\bigg)\bigg((2)(x_i-\theta)(-1)\bigg) - \bigg((x_i-\theta)^2 - 1\bigg)\bigg(2\big(1+(x_i - \theta)^2\big)(2)(x_i - \theta)(-1)\bigg)}{\bigg(\big(1+(x_i - \theta)^2\big)^2\bigg)^2}$$

$$f''(\theta) = 2\sum_{i=1}^{n} \frac{-2(x_i-\theta)\big(1+(x_i -\theta)^2\big)\big[\big(1+(x_i - \theta)^2\big) - 2\big((x_i-\theta)^2 - 1\big)\big]}{\big(1+(x_i - \theta)^2\big)^4}$$
$$f''(\theta) = 2\sum_{i=1}^{n} \frac{-2(x_i-\theta)\big(1+(x_i -\theta)^2\big)\big[3-(x_i-\theta)^2\big]}{\big(1+(x_i - \theta)^2\big)^4}$$
$$f''(\theta) = 2\sum_{i=1}^{n} \frac{-2(x_i-\theta)\big[3-(x_i-\theta)^2\big]}{\big(1+(x_i - \theta)^2\big)^3}$$
Now, we can start implementing in R, we need to create an R function for $f(\theta)$, $f'(\theta)$,  $f''(\theta)$ and the general function for Halley's method.