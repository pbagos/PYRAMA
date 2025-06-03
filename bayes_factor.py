import numpy as np
import math
from scipy.stats import norm


# Alternative Bayes Factor
def bayes_factor(or_val, se):
    target_prob = 0.95  # with this probability
    b = np.log(or_val)
    tau_seq = np.arange(0.0000001, 1.0001, 0.001)  # grid to evaluate tau
    # Calculate the CDF using scipy.stats.norm.cdf
    x = 1 - norm.cdf(b, loc=0, scale=tau_seq)
    tau = tau_seq[np.argmin(np.abs(x - (1 - target_prob) / 2))]  # which is closest to target prob?
    # # Calculate the standard deviation of the normal distribution
    # standard_deviation = math.sqrt(tau ** 2 + se ** 2)
    #
    # # Calculate the log density of the normal distribution
    # log_density = norm.logpdf(b, loc=0, scale=standard_deviation) - norm.logpdf(b, loc=0, scale=se)
    #
    # bf = np.exp(log_density)
    bf = np.sqrt(se ** 2 / (tau ** 2 + se ** 2)) * np.exp(
        0.5 * (or_val ** 2 / se ** 2) * (tau ** 2 / (tau ** 2 + se ** 2)))
    return bf

# # Example usage
# b_est = 5  # Replace this with the value of b_est in your code
# se = 10# Replace this with the value of se in your code
#
# result_bf = bayes_factor(b_est, se)
# print(result_bf)
