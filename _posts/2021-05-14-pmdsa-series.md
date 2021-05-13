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


        ```{r}
        # Create function for f(theta)
        compute_f_cauchy_log_likelihood <- function(theta, x){
        numerator = x - theta
        denominator = 1 + (x - theta)**2
        value <- 2*sum(numerator/denominator)
        return(value)
        }

        # Create function for f'(theta)
        compute_fprime_cauchy_log_likelihood <- function(theta, x){
        numerator = (x - theta)**2 - 1
        denominator = (1 + (x - theta)**2)**2
        value <- 2*sum(numerator/denominator)
        return(value)
        }

        # Create function for f''(theta)
        compute_fdoubleprime_cauchy_log_likelihood <- function(theta, x){
        numerator = -2*(x - theta)*(3 - (x - theta)**2)
        denominator = (1 + (x - theta)**2)**3
        value <- 2*sum(numerator/denominator)
        return(value)
        }

        # Implement the main function for Halley's method
        find_root_using_halleys_method <- function(theta0, max_iters, tol, f, f_prime, f_double_prime, x){
        for (i in seq(0, max_iters, 1)){
            f_value <- f(theta0, x)
            f_prime_value <- f_prime(theta0, x)
            f_double_prime_value <- f_double_prime(theta0, x)  
            theta1 <- theta0 - (2*f_value*f_prime_value)/(2*f_prime_value**2 - f_value*f_double_prime_value)
            theta_diff <- abs(theta1-theta0) 
            print(paste("Iter:", i, "Theta:", theta1, "Diff:", theta_diff))
            
            #catch divergent cases
            if (is.infinite(theta1)){
            print("Value has diverged")
            return(NULL)
            }
            if (theta_diff <= tol) {
            return(theta1)
            }
            theta0 <- theta1
        }
        return(theta1)
        }
        ```

Now, let's start running the Halley's mehod algorithm to find the MLE using the different starting points theta = -11, 0, 1.5, 4, 4.7, 7,38, and 51. 

        ```{r}
        #when start point is theta = -11
        find_root_using_halleys_method(theta0 = -11, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

It resulted to Nan.

        ```{r}
        #when start point is theta = 0
        find_root_using_halleys_method(theta0 = 0, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

It converged to $\theta=-0.1922866$.

        ```{r}
        #when start point is theta = 1.5
        find_root_using_halleys_method(theta0 = 1.5, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```
It converged to 1.713587

        ```{r}
        #when start point is theta = 4
        find_root_using_halleys_method(theta0 = 4, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

        ```{r}
        #when start point is theta = 4.7
        find_root_using_halleys_method(theta0 = 4.7, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

        ```{r}
        #when start point is theta = 7
        find_root_using_halleys_method(theta0 = 7, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

        ```{r}
        #when start point is theta = 38
        find_root_using_halleys_method(theta0 = 38, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

        ```{r}
        #when start point is theta = 51
        find_root_using_halleys_method(theta0 = 51, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```
        ```{r}
        #when start point is theta = mean(samples)
        mean_samples <- mean(samples)
        print(paste("Mean of Samples", mean_samples))
        find_root_using_halleys_method(theta0 = mean_samples, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```
        ```{r}
        #when start point is theta = median(samples)
        mean_samples <- median(samples)
        print(paste("Mean of Samples", mean_samples))
        find_root_using_halleys_method(theta0 = mean_samples, max_iters = 1000, tol = 1e-8, f = compute_f_cauchy_log_likelihood, f_prime = compute_fprime_cauchy_log_likelihood, f_double_prime = compute_fdoubleprime_cauchy_log_likelihood, x = samples)
        ```

From above, we see that we get different results on the MLE values for the log likelihood function of Cauchy distribution using Halley's method. If the starting point is -11 or 7, the value of theta diverges to infinity. Otherwise for other starting points, we get different values of the final theta namely: -0.1922866, 1.713587, 2.817472, 41.04085, 56.25336. To verify if these are all valid values of theta MLE, let's plot the $f(\theta)$ and verify if these derived values are actually roots of $f(\theta)= 2\sum_{i=1}^{n} \frac{(x_i - \theta)}{1+(x_i-\theta)^2} = 0$.

        ```{r}
        # plot f(theta) 
        thetas <- seq(from=-100, to=100, by=0.001)
        first_deriv_log_likelihood_values <- unlist(lapply(thetas, compute_f_cauchy_log_likelihood, x=samples))


        plot(thetas, first_deriv_log_likelihood_values, type="l", col="blue", xlab=TeX("theta($\\theta$)"), ylab=TeX("f($\\theta$)"))
        lines(thetas, rep(0, length(thetas)), type="l", col="green",)
        ```

We can see from the plot above that the derived values of $\theta MLE$ are the roots of $f(\theta)$, the intersection of the blue($f(\theta)$) and green line($f(\theta)=0$). We can see that around $\theta=0$ , there are 3 possible roots, which correspond to what we got in the runs: -0.1922866, 1.713587, 2.817472. We also have 4 other roots around 50, and we found 2 of them in our runs: 41.04085, 56.25336. 


## Problem 2

(Algorithm) Twenty values were assumed to be observed from a continuous distribution with mean $\theta$. The values obtained are as follows:


| 102 | 112 | 131 | 107 | 114 | 95  | 133 | 145 | 139 | 117 |
| 93  | 111 | 124 | 122 | 136 | 141 | 119 | 122 | 151 | 143 |


Suppose now that, as in a simulation, we have the option of continually generating additional data values. How many additional observations do you think will we have to generate if we want to be 99% certain that our final estimate of $\theta$ is correct within ±0.5 units? Discuss how you will proceed (include assumptions and/or reasons), present the necessary steps (as in an algorithm),and implement in $R$.

## Solution to Problem 2

The formula for confidence interval for population mean is given by: 

$$\bar{X} \pm Z_{\alpha/2}\frac{\sigma}{\sqrt{n}}$$

provided that data follows a normal distribution with population mean $\theta$ and population standard deviation $\sigma$. However, the given data provided in Problem 2 did not provide a value for $\sigma$. Thus we need to find a way to use the sample standard deviation $s$ instead. Hence, we can assume that the $\frac{X-\theta}{\sigma/sqrt(n)}$ follows a student T distribution instead of the standard normal distribution $Z$, with confidence interval for the population mean $\theta$ given by:

$$\bar{X} \pm t_{\alpha/2}\frac{s}{\sqrt{n}}$$

with degrees of freedom = n-1.

The problem asks how many additional samples we should generate if we want to have a margin of error $\pm 0.5$ with 99% confidence level($\alpha=0.01$). In other words, what is the value of $n$ for which the margin of error is `E=0.5`?

Here is the algorithm:

  1. Starting with $n=20$, we compute for the margin of error $E = t_{\alpha/2}\frac{s}{\sqrt{n}}$.
  2. If $E \le 0.5$: stop and return the value of $n$ else we increment the value of $n$ and evaluate again the margin of error. Also, in this algorithm, we assume that the sample standard deviation does not change as we increase n.
  3. We continue this incremental process and stop until n satisfy the  $E \le 0.5$. 


        ```{r}
        values <- c(102, 112, 131, 107, 114, 95, 133, 145, 139, 117, 93, 111, 124, 122, 136, 141, 119, 122, 151, 143)
        sample_std_dev = sd(values)
        paste("Sample Standard Deviation: ", sample_std_dev)

        error_thresh <- 0.5
        alpha = 0.01 #since 99% confidence level
        n = length(values)
        error = qt(1-alpha/2, n-1)*(sample_std_dev/sqrt(n))

        while (error > error_thresh){
        n <- n + 1
        error <- qt(1-alpha/2, n-1)*(sample_std_dev/sqrt(n))
        if (n%%1000 == 0){
            print(paste("Sample size n:", n, "|Margin of error:", error))
        }
        }
        # Final value of sample size
        print(paste("Sample size n:", n, "|Margin of error:", error))
        ```

Thus, it requires to have 7495 samples to reach a marginal error of less than or equal to 0.5. Therefore, we need to generate 7495-20 = 7475 more additional samples to make sure that we are 99% confident that our final estimate for population mean $\theta$ is correct within ±0.5 units.


## Problem 3

Data Generating Process or DGP) Other than MCAR, MAR, and NMAR classification on
missingness, literature (e.g. Little and Rubin (2014), Statistical Analysis with Missing Data, 2nd ed.) identifies MCAR and MAR as ignorable missingness, and NMAR as non-ignorable missingness (also, missing not at random or MNAR, cf. Im and Kim, 2017 in https://www.sciencedirect.com/science/article/pii/S1226319217300406 ; note Science Direct available only within DilNet or OpenAthens; see also simpler discussions in https://en.wikipedia.org/wiki/Missing_data ). The main difference between ignorable and non-ignorable missing data is on whether missingness depends on the missing values themselves. Non-ignorable missingness essentially suggests that the missing data should not be ignored in any analysis, on the premise that the missingness (of a data) and/or nonresponse (of a subject) are likely to be correlated with the unobserved or missing values. That is, the probability of missingness depends on the (potentially missing) variable itself. The cases with missing data therefore differ from cases with complete data for some reason rather than completely or conditionally random. Under this case, information on the missing variable and missing mechanism is required to obtain a consistent statistical analysis results (Im and Kim, 2017). In contrast, for MCAR or MAR, the missing values do not depend on the unobserved values, but are related to other variables (outside the data for MCAR, other variables in the data for MAR).

a. Suggest and/or recommend a framework (i.e., characterization, requirements, assumptions) and an algorithm (i.e., process flow, step-by-step procedure) for the DGP that illustrates characteristics of a (univariate) data with non-ignorable missing values.
b. Implement an R code for your recommendation in (a). Consider creating a function call with appropriate parameters. The generated data must be one of the function’s output objects.
c. Update your code/implementation in (b) for a DGP of a multivariate data (i.e., several random variables) with non-ignorable missing values.

## Solution to Problem 3

*a. Suggest and/or recommend a framework (i.e., characterization, requirements, assumptions) and an algorithm (i.e., process flow, step-by-step procedure) for the DGP that illustrates characteristics of a (univariate) data with non-ignorable missing values.*

For problem a, to create a DGP that illustrates the characteristics of a univariate data with non-ignorable missing values, we need to describe the data from a distribution perspective on what it means to have a non-ignorable missing values, that is, experiencing MNAR missingness. Though, identifying if data has MNAR missing values requires domain expertise in the field of study, I think we can set some obvious properties, for example:

  1. The missing values happens in gaps or clusters of values of the distribution.
  2. Missing data is usually observed on the extreme values(lowest or highest values)

Now on the algorithm, to simplify things, we will only assume generating a standard normal distributed data. Here are the steps:

  1. Generate the N random samples from the distribution
  2. Require an interval of values for which the data will be missing for MNAR, i.e. x> certain value, and set it to NA.
  3. To make it more realistic, also assign some of the values to be missing via MCAR, that is randomly draw from the drawn samples in step 1, and assign them as MCAR(ignorable missingness).


*b. Implement an R code for your recommendation in (a). Consider creating a function call with appropriate parameters. The generated data must be one of the function’s output objects.*


        ```{r}
        # Create the function for data generating process of data including the missing values
        dgp_standard_normal_univariate <- function(num_samples, mnar_lower_bnd=NULL, mnar_upper_bnd=NULL, mcar_prob=0){

        #draw samples from the distribution
        data_df <- data.frame(data=rnorm(n=num_samples))
        
        #define the MNAR missing data
        if (is.null(mnar_lower_bnd) & is.null(mnar_upper_bnd)){
            stop("At least either mnar_lower_bnd or mnar_upper_bnd should be provided")
        }else if (is.null(mnar_lower_bnd)){
            data_df$is_mnar <- data_df$data<=mnar_upper_bnd
        }else if (is.null(mnar_upper_bnd)){
            data_df$is_mnar <- data_df$data>=mnar_lower_bnd
        }else{
            data_df$is_mnar <- (data_df$data>=mnar_lower_bnd) & (data_df$data<=mnar_upper_bnd)
        }
        
        #define the MCAR missing data
        mcar_indices = sample(which(!data_df$is_mnar), size=as.integer(mcar_prob*num_samples))
        data_df$is_mcar <- 1:num_samples %in% mcar_indices
        
        #create a stacked format of the data
        obs_bool_index <- (!data_df$is_mcar)&(!data_df$is_mnar)
        mcar_bool_index <- data_df$is_mcar
        mnar_bool_index <- data_df$is_mnar
        
        stacked_data_df <- data.frame(data = c(data_df$data[obs_bool_index], 
                                                data_df$data[mcar_bool_index], 
                                                data_df$data[mnar_bool_index]),
                                        type = c(rep("Observed data", sum(obs_bool_index)),
                                                rep("MCAR data", sum(mcar_bool_index)),
                                                rep("MNAR data", sum(mnar_bool_index))))
        
        return(stacked_data_df)
        }

        #create a function for visualization
        visualize_dgp<-function(stacked_data_df){
        num_samples <- nrow(stacked_data_df)
        
        # visualize the data
        plot1 <- ggplot(data=stacked_data_df, aes(x=data, fill=rep("Complete Data", num_samples))) +
            geom_histogram(colour = 'white', binwidth = 0.15) +
            labs(fill="Data Type") +
            scale_fill_manual(values="gray")
        
        #get the y range to be used to standardize the y limits for all the plots
        yrange <- ggplot_build(plot1)$layout$panel_params[[1]]$y.range
        xrange <- ggplot_build(plot1)$layout$panel_params[[1]]$x.range
        plot1 <- plot1 +
            coord_cartesian(ylim=c(yrange[1], yrange[2]), xlim=c(xrange[1], xrange[2]))
        
        
        plot2 <- ggplot(data=subset(stacked_data_df, type=="Observed data"), aes(x=data, fill=type)) +
            geom_histogram(colour = 'white', binwidth = 0.15) +
            labs(fill="Data Type") +
            scale_fill_manual(values= "#E7B800") +
            coord_cartesian(ylim=c(yrange[1], yrange[2]), xlim=c(xrange[1], xrange[2]))
        
        plot3 <- ggplot(data=stacked_data_df, aes(x=data, fill=type)) +
            geom_histogram(colour = 'white', binwidth = 0.15) +
            labs(fill="Data Type")  +
            scale_fill_manual(values=c("#00AFBB", "#FC4E07", "#E7B800")) +
            coord_cartesian(ylim=c(yrange[1], yrange[2]), xlim=c(xrange[1], xrange[2]))
        
        print(plot1)
        print(plot2)
        print(plot3)
        }
        ```


Now let's take a look at some examples. Suppose, we want to generate N=50000 samples from standard normal distribution with MNAR missingness for values between 1 and 1.5 and MCAR missingness with size probability of 10%.


        ```{r fig.show="hold", out.width="33%", fig.height=5}
        set.seed(1234)
        num_samples = 50000
        mnar_lower_bnd = 1
        mnar_upper_bnd = 1.5
        mcar_prob = 0.10

        #generate the data
        stacked_data_df <- dgp_standard_normal_univariate(num_samples=num_samples, mnar_lower_bnd=mnar_lower_bnd, mnar_upper_bnd=mnar_upper_bnd, mcar_prob=mcar_prob)

        #show some sample values
        head(stacked_data_df)

        visualize_dgp(stacked_data_df)
        ```

The leftmost plot shows the complete data. The center figure shows only the observed data when MNAR and MCAR missingness was introduced. The rightmost plot shows the breakdown of the observed data in yellow color, the MCAR data in cyan and the MNAR data in bright red color.

Let's show another example where there is MNAR data for values less than -1 and 15% MCAR data with the same number of samples N=50000


        ```{r fig.show="hold", out.width="33%", fig.height=5}
        set.seed(1234)
        num_samples = 50000
        mnar_upper_bnd = -1
        mcar_prob = 0.15

        #generate the data
        stacked_data_df <- dgp_standard_normal_univariate(num_samples=num_samples, mnar_upper_bnd=mnar_upper_bnd, mcar_prob=mcar_prob)

        visualize_dgp(stacked_data_df)
        ```


*c. Update your code/implementation in (b) for a DGP of a multivariate data (i.e., several random variables) with non-ignorable missing values.*

        ```{r}
        # Create the function for data generating process of data including the missing values
        dgp_standard_normal_multivariate <- function(num_samples, num_dimensions, mnar_lower_bnds, mnar_upper_bnds, mcar_prob){
        # generate the data 
        data <- lapply(1:num_samples, function(x) rnorm(num_dimensions))
        data <- do.call(rbind, data)
        # define the MNAR missing data
        mnar_upper_bnd_matrix <- matrix(mnar_upper_bnds, nrow=num_samples, 
                                        ncol=length(mnar_upper_bnds), byrow=TRUE)
        mnar_lower_bnd_matrix <- matrix(mnar_lower_bnds, nrow=num_samples, 
                                        ncol=length(mnar_lower_bnds), byrow=TRUE)
        # indicator matrix for MNAR data
        is_mnar_matrix <- (data<=mnar_upper_bnd_matrix) & (data>=mnar_lower_bnd_matrix)
        # define the MCAR missing data
        valid_indices_to_sample <- which(!is_mnar_matrix, arr.ind = T)
        chosen_indices <- sample(1:nrow(valid_indices_to_sample), size=as.integer(mcar_prob*num_samples*num_dimensions))
        #indicator matrix for MCAR data
        is_mcar_matrix <- matrix(F, nrow=num_samples, ncol=num_dimensions)
        for (i in chosen_indices){
            row_col <- valid_indices_to_sample[i,]
            is_mcar_matrix[row_col[1], row_col[2]] <- T
        }
        return(list(data, is_mnar_matrix, is_mcar_matrix))
        }
        ```

Now, let's test the DGP for multivariate standard normal distribution with dimension=3. On dimension 1 MNAR is observed on $x\ge1.5$, on dimension 2, in the interval $[0, 1]$ and on the third dimension on $x\le-2$. Also, an MCAR with size propbability of 10%

        ```{r}
        set.seed(1234)
        num_samples = 50000
        num_dimensions = 3
        mnar_upper_bnds <- c(Inf, 1, -2)
        mnar_lower_bnds <- c(1.5, 0, -Inf)
        mcar_prob = 0.1

        output <- dgp_standard_normal_multivariate(num_samples, num_dimensions, mnar_lower_bnds, mnar_upper_bnds, mcar_prob)

        data <- output[[1]]
        is_mnar_matrix <- output[[2]]
        is_mcar_matrix <- output[[3]]

        #the final data generated with missing data
        indicator_matrix <- matrix(NA, nrow=num_samples, ncol=num_dimensions)
        indicator_matrix[!(is_mnar_matrix | is_mcar_matrix)] = 1
        final_data <- data*indicator_matrix
        head(final_data)
        ```

        ```{r fig.show="hold", out.width="33%", fig.height=5}
        #let's visualize each dimension and validate 
        #create a stacked format of the data
        #on dimension 1
        dimension <- 1
        obs_bool_index <- (!is_mcar_matrix[,dimension])&(!is_mnar_matrix[,dimension])
        mcar_bool_index <- is_mcar_matrix[,dimension]
        mnar_bool_index <- is_mnar_matrix[,dimension]

        stacked_data_df_dimension1 <- data.frame(data = c(data[,dimension][obs_bool_index], 
                                            data[,dimension][mcar_bool_index], 
                                            data[,dimension][mnar_bool_index]),
                                    type = c(rep("Observed data", sum(obs_bool_index)),
                                            rep("MCAR data", sum(mcar_bool_index)),
                                            rep("MNAR data", sum(mnar_bool_index))))
        head(stacked_data_df_dimension1)

        # visualize the data
        visualize_dgp(stacked_data_df_dimension1)
        ```

        ```{r fig.show="hold", out.width="33%", fig.height=5}
        #let's visualize each dimension and validate 
        #create a stacked format of the data
        #on dimension 2
        dimension <- 2
        obs_bool_index <- (!is_mcar_matrix[,dimension])&(!is_mnar_matrix[,dimension])
        mcar_bool_index <- is_mcar_matrix[,dimension]
        mnar_bool_index <- is_mnar_matrix[,dimension]

        stacked_data_df_dimension2 <- data.frame(data = c(data[,dimension][obs_bool_index], 
                                            data[,dimension][mcar_bool_index], 
                                            data[,dimension][mnar_bool_index]),
                                    type = c(rep("Observed data", sum(obs_bool_index)),
                                            rep("MCAR data", sum(mcar_bool_index)),
                                            rep("MNAR data", sum(mnar_bool_index))))
        head(stacked_data_df_dimension2)

        # visualize the data
        visualize_dgp(stacked_data_df_dimension2)
        ```

        ```{r fig.show="hold", out.width="33%", fig.height=5}
        #let's visualize each dimension and validate 
        #create a stacked format of the data
        #on dimension 3
        dimension <- 3
        obs_bool_index <- (!is_mcar_matrix[,dimension])&(!is_mnar_matrix[,dimension])
        mcar_bool_index <- is_mcar_matrix[,dimension]
        mnar_bool_index <- is_mnar_matrix[,dimension]

        stacked_data_df_dimension3 <- data.frame(data = c(data[,dimension][obs_bool_index], 
                                            data[,dimension][mcar_bool_index], 
                                            data[,dimension][mnar_bool_index]),
                                    type = c(rep("Observed data", sum(obs_bool_index)),
                                            rep("MCAR data", sum(mcar_bool_index)),
                                            rep("MNAR data", sum(mnar_bool_index))))
        head(stacked_data_df_dimension3)

        # visualize the data
        visualize_dgp(stacked_data_df_dimension3)
        ```
