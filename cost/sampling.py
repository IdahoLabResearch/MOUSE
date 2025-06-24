import numpy as np

def create_lognormal_sampler(low_cost, high_cost, class3_cost):
    # Calculate the natural logarithms of the given costs
    ln_low_cost = np.log(low_cost)
    ln_high_cost = np.log(high_cost)
    ln_class3_cost = np.log(class3_cost)

    # Calculate the mean (mu) and standard deviation (sigma) of the logarithms
    mu = np.mean([ln_low_cost, ln_high_cost, ln_class3_cost])
    sigma = np.std([ln_low_cost, ln_high_cost, ln_class3_cost], ddof=0) # Population std deviation

    # Define the sampler function
    def sampler():
        return np.random.lognormal(mean=mu, sigma=sigma)

    return sampler()

def truncated_normal_sample(mean, std, lower_bound, upper_bound):
    while True:
        sample = np.random.normal(mean, std)
        if lower_bound <= sample <= upper_bound:
            return sample

def uniform_sample(low, high):
    return np.random.uniform(low, high)

def sampler(distribution, **kwargs):
    if distribution == "Lognormal":
        return create_lognormal_sampler(kwargs['low_cost'], kwargs['high_cost'], kwargs['class3_cost'])
    elif distribution == "Truncated Normal":
        return truncated_normal_sample(kwargs['mean'], kwargs['std'], kwargs['lower_bound'], kwargs['upper_bound'])
    elif distribution == "Uniform":
        return uniform_sample(kwargs['low'], kwargs['high'])    
    else:
        raise ValueError("Unavailable Distribution")
