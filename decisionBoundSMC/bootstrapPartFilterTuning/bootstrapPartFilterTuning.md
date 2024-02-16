---
title: Initial fine-tuning of bootstrap particle filter
author: Cas Ludwig
date: 16-02-2024
---

# Preliminary runs on toy data

In this toy dataset, the agent stepped through 9 different parameter combinations (3 levels for $\alpha \times$ 3 levels for $\beta$, keeping $\eta = 2.5$​ constant). The dynamics of these parameter changes were non-sensical, but served the development of the bootstrap particle filter. In this document, I have summarised the results of some tuning of this filter.

## Tuning of the bootstrap filter

- In all cases, we kept the prior means constant at $\mu_0 = \lbrace 5, 0, 1 \rbrace$.

- We varied the variances of the multivariate Gaussian prior, $\sigma^2_0$ and of the multivariate Gaussian transition distribution, $\sigma^2_q$​. 

- We varied the number of particles $n \in \lbrace 2000, 5000 \rbrace$. 

In total, there are 10 runs (ok, 12) and the outputs are shown as 'prelimRunX.png', where X goes from 1 - 10. It does not help that these figures do not have exactly the same fixed y-axes, but we can nevertheless draw some conclusions from them.

In the table below, a triplet of parameters always relate to $\lbrace \alpha, \beta, \eta \rbrace$, respectively. Note that we show the standard deviations, rather than the variances. 

| Run  | $n$  | $\sigma_0$                  | $\sigma_q$                        | Comments                                                     |
| ---- | ---- | --------------------------- | --------------------------------- | ------------------------------------------------------------ |
| 1    | 2000 | $\lbrace 5, 0.4, 3 \rbrace$ | $\lbrace 2, 0.2, 0.5 \rbrace$     | Parameter $\alpha$ tracks well, but estimates are quite variable (CIs around 10). Parameter $\beta$ is largely flat, but centred around its average value and therefore actually fairly close. Again, variability is large: CIs generally include 0, even though the gradient was consistently negative. |
| 2    | 5000 | $\lbrace 5, 0.4, 3 \rbrace$ | $\lbrace 2, 0.2, 0.5 \rbrace$     | increasing the number of particles does not seem to do much. |
| 3    | 2000 | $\lbrace 3, 0.2, 1 \rbrace$ | $\lbrace 1, 0.1, 0.1 \rbrace$     | Parameter $\alpha$ again tracks quite well, although it does not quite reach the higher intercepts. This is probably the price of having a tighter transition function: the particles are not allowed to take very large steps, so it takes more trials to adapt to a sudden jump. Parameter $\beta$ seems to track a little better compared to Runs 1 and 2, but the variability is still large and the CIs generally include 0. Overall, the estimates appear somewhat less variable than in Runs 1 and 2; this is particularly noticeable in the estimates for $\eta$. The $\eta$ estimates themselves also appear more stable (albeit consistently underestimating the true value). |
| 4    | 5000 | $\lbrace 3, 0.2, 1 \rbrace$ | $\lbrace 1, 0.1, 0.1 \rbrace$     | Tracking seems slightly better for all three parameters, but this improvement is most noticeable for $\beta$. CIs about as tight as in Run 3, so better than in Runs 1 - 2. |
| 5    | 2000 | $\lbrace 3, 0.2, 1 \rbrace$ | $\lbrace 0.5, 0.05, 0.05 \rbrace$ | Tracking is clearly off here, apart from $\eta$, which is of course the only parameter that does not change. CIs are pretty tight, but there are now many timepoints where they do not include the true values. |

It appears that the variance of the transition distribution is important for tracking: set it too small and the particles just cannot jump to the region of the space where the likelihood is high(er) or, at least, they can only get there with small steps (i.e. slowly). However, small variance in the transition distribution seems to result in tighter estimates and for parameters that do not change over time ($\eta$), estimation is improved. So the $\sigma_q$​ hyperparameters control the trade-off between posterior variance and tracking performance (at least, for non-stationary parameters). The number of particles doesn't seem to have a huge effect, but in one case (Run 4) the higher number of particles seemed to give slightly better tracking.

Note that Runs 1-2 on the one hand and Runs 3-4 on the other hand confound the prior variance with the transition variance. Therefore, the next two runs decouple these sets of parameters.

| Run  | $n$  | $\sigma_0$                  | $\sigma_q$                    | Comments                                                     |
| ---- | ---- | --------------------------- | ----------------------------- | ------------------------------------------------------------ |
| 6    | 2000 | $\lbrace 5, 0.4, 3 \rbrace$ | $\lbrace 1, 0.1, 0.1 \rbrace$ | Compare with Run 1 to get a sense of the role of the transition function. Tracking seems pretty similar as in Run 1. It looks like $\eta$ is more stable around its true value. CIs are generally a bit tighter, especially later in time (particularly noticeable for $\eta$). Compare with Run 3 to get a sense of the role of the prior variance. Tracking looks very similar, as do the CIs. Again, the only noticeable improvement is for $\eta$. However, that might just be the kind of variation with repeated runs without any changes in the hyperparameters. I'll address this below. |
| 7    | 2000 | $\lbrace 3, 0.2, 1 \rbrace$ | $\lbrace 2, 0.2, 0.5 \rbrace$ | Comparison with Run 1 gives a sense of the influence of the prior variance. There are no noticeable differences between these runs. Comparison with Run 3 gives insight into what happens when we increase the transition variance. It appears that our estimates become more variable, but seem to track a bit better. |

In summary, the impressions formed above are confirmed. Manipulating the prior variance does not do a huge amount. I guess that we just need sufficient variance so that at least some part of the prior covers the high likelihood region of the parameter space. Then particles will be drawn to that region quickly enough. In our case, the prior variance was probably sufficiently high in all these runs, so it makes sense that it did not affect the results very much. Indeed, when we shrink the prior variance right down (Run 8: $\sigma_0 = \lbrace 0.5, 0.01, 0.1 \rbrace$; $\sigma_q = \lbrace 1, 0.1, 0.1 \rbrace$), we see that it takes longer for the estimates to converge to the true values, which is problematic if the prior means are some way off (e.g. as is the case for $\eta$). 

Finally, I wanted to gain some insight into the extent to which the estimates and their CIs vary across repeated runs with the same hyperparameters. Based on the results so far, I decided to set the prior quite wide, and to set the variance of the transitions a little bit wider than in my very original fits (corresponding to Run 3 here). So for Runs 9a and 9b, we have: $\sigma_0 = \lbrace 5, 0.4, 0.3 \rbrace$; $\sigma_q = \lbrace 1.5, 0.15, 0.25 \rbrace$, i.e. in between the values for Runs 1 and 3. These two runs were identical. In addition to the consistency between Runs 9a and 9b, we're interested in how the new set of parameters for $\sigma_q$ trades off tracking performance with posterior variance.

Starting with the latter, comparison of Run(s) 9 with Runs 1 and 6 provides an assessment of the changing transition variance. There does not seem to be a big difference between Runs 9 and 1, suggesting that a further increase in transition variance (from Run 9 to Run 1) does not improve tracking. Compared with Run 6, tracking of the time-varying parameters seems a little improved, in line with the comparison between Runs 1 and 6. However, again the posterior variance is greater in Run 9 and the "tracking" of the stationary parameter $\eta$ is worse.

Comparison between Runs 9a and 9b suggests that the estimates and their CIs are highly replicable. There are some apparent differences, but they are really very minor. One question is whether this type of variation is eliminated by increasing the number of particles. So I repeated these fits with $n = 5000$ (Runs 10a and 10b). Again, there is some very slight variation, but nothing that detracts from the quality of the tracking performance and recovery of the ground-truth parameters. Perhaps the variation is slightly smaller here than in Run 9, but there really is not much in it. Especially if computation time is an issue, there seems little reason to have $n > 2000$. 

