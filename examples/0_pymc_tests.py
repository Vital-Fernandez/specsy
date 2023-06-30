import pymc as pm
import arviz as az
from matplotlib import pyplot as plt

# Assume 10 trials and 5 successes out of those trials
# Change these numbers to see how the posterior plot changes
trials = 10
successes = 5

print(pm.__version__)
print(az.__version__)

# Set up model context
with pm.Model() as coin_flip_model:
    # Probability p of success we want to estimate
    # and assign Beta prior
    p = pm.Beta("p", alpha=1, beta=1)

    # Define likelihood
    obs = pm.Binomial("obs", p=p, n=trials, observed=successes)

    # Hit Inference Button
    idata = pm.sample(cores=1, chains=4)

az.plot_posterior(idata, show=True)

plt.show()


