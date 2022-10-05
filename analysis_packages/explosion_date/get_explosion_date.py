import utils.format_data

def fetch_table_element(colname, table):
    '''
    to avoid .data vs .data.data nonsense
    '''
    if type(colname) == str:
        if type(table[colname].data.data) == memoryview:
            dat_ = table[colname].data
        else:
            dat_ = table[colname].data.data
    elif type(colname) == list:
        dat_ = []
        for col in colname:
            dat_.append(fetch_table_element(col, table))
    return dat_

def read_in_lc(lc_dir):
    '''
    Read in target lightcurve from directory.
    '''
    tab = Table.read(f'{lc_dir}/{ZTF_name}.txt', format ='ascii')
    jd, filt_, mag, sigma_mag, limmag, instrument = fetch_table_element(['jdobs', 'filter', 'magpsf', 'sigmamagpsf', 'limmag', 'instrument'], tab)
    m = (filt_ == 'g')  & (instrument == 'P48+ZTF') & (mag > 16)
    jd_, filt_, mag_, sigma_mag_, limmag_, instrument_  = jd[m], filt_[m], mag[m], sigma_mag[m], limmag[m], instrument[m]
    
    return jd_, filt_, mag_, sigma_mag_, limmag_, instrument_

def calibrate_target(mag, sigma_mag, limmag):
    '''
    Calibrate to eliminate mag == 99.
    '''
    is_upper_limit_ = mag == 99
    mag[is_upper_limit_] = limmag[is_upper_limit_]
    flux_obs_ = 10**(-mag/2.5)/np.max(10**(-mag[~is_upper_limit_]/2.5))
    flux_err_ = sigma_mag/1.09*flux_obs_
    var_mag_ = (sigma_mag/1.09)**2
    var_ = var_mag_/1.09*flux_obs_
    
    return is_upper_limit_, flux_obs_, flux_err_, var_mag_, var_

# For fitting

def scale_target_times(jd, is_upper_limit):
    t_d1 = np.min(jd[~is_upper_limit])
    time_zero_ = np.min(jd[~is_upper_limit])-15
    scaled_time_ = jd - time_zero_
    scaled_t_d1_ = t_d1 - time_zero_
    t_end = 20 #arbitrary choice for end date
    times_to_fit_ = scaled_time_ < t_end
    
    return time_zero_, scaled_time_, scaled_t_d1_, times_to_fit_

def calc_flux(time, fm, t0, te):
    flx = fm*(1 - np.exp(-(time - t0)/te))
    flx[time < t0] = 0
    return flx


# In[20]:


'''
Read in target and calibrate it to only get directly pre-/post-explosion data.
Convert to flux.
'''

def fetch_table_element(colname, table):
    '''
    to avoid .data vs .data.data nonsense
    '''
    if type(colname) == str:
        if type(table[colname].data.data) == memoryview:
            dat_ = table[colname].data
        else:
            dat_ = table[colname].data.data
    elif type(colname) == list:
        dat_ = []
        for col in colname:
            dat_.append(fetch_table_element(col, table))
    return dat_

def read_in_lc(lc_dir):
    '''
    Read in target lightcurve from directory.
    '''
    tab = Table.read(f'{lc_dir}/{ZTF_name}.txt', format ='ascii')
    jd, filt_, mag, sigma_mag, limmag, instrument = fetch_table_element(['jdobs', 'filter', 'magpsf', 'sigmamagpsf', 'limmag', 'instrument'], tab)
    m = (filt_ == 'g')  & (instrument == 'P48+ZTF') & (mag > 16)
    jd_, filt_, mag_, sigma_mag_, limmag_, instrument_  = jd[m], filt_[m], mag[m], sigma_mag[m], limmag[m], instrument[m]
    
    return jd_, filt_, mag_, sigma_mag_, limmag_, instrument_

def calibrate_target(mag, sigma_mag, limmag):
    '''
    Calibrate to eliminate mag == 99.
    '''
    is_upper_limit_ = mag == 99
    mag[is_upper_limit_] = limmag[is_upper_limit_]
    flux_obs_ = 10**(-mag/2.5)/np.max(10**(-mag[~is_upper_limit_]/2.5))
    flux_err_ = sigma_mag/1.09*flux_obs_
    var_mag_ = (sigma_mag/1.09)**2
    var_ = var_mag_/1.09*flux_obs_
    
    return is_upper_limit_, flux_obs_, flux_err_, var_mag_, var_

# For fitting

def scale_target_times(jd, is_upper_limit):
    t_d1 = np.min(jd[~is_upper_limit])
    time_zero_ = np.min(jd[~is_upper_limit])-15
    scaled_time_ = jd - time_zero_
    scaled_t_d1_ = t_d1 - time_zero_
    t_end = 20 #arbitrary choice for end date
    times_to_fit_ = scaled_time_ < t_end
    
    return time_zero_, scaled_time_, scaled_t_d1_, times_to_fit_

def calc_flux(time, fm, t0, te):
    flx = fm*(1 - np.exp(-(time - t0)/te))
    flx[time < t0] = 0
    return flx


# In[27]:


'''
Functions for MCMC
'''

def log_likelihood(theta, time_values, f_obs, f_err, var, is_upper_limit):
    '''
    include intrinsic var
    '''
    fm, t0, te = theta
    flux_pred = calc_flux(time = time_values, fm = fm, t0 = t0, te = te)
    lnL_array = np.log(2*np.pi*np.sqrt(f_err**2 + var)**2) + (flux_pred - f_obs)**2/np.sqrt(f_err**2 + var)**2 
    lnL_array[is_upper_limit & (flux_pred < f_obs)] =0 
    lnL = -0.5*np.sum(lnL_array )
    return lnL

def log_prior(theta):
    fm, t0, te = theta
    if not ((fm > 0.5) and (fm < 1.5)):
        return -np.inf
    if not ((t0 > (scaled_t_d1 - 10)) and (t0 < scaled_t_d1)):
        return -np.inf
    if not ((te > 0) and (te < 10)):
        return -np.inf
    return 0


def log_posterior(theta, time_values, f_obs, f_err, var, is_upper_limit):
    '''
    with intrinsic var
    '''
    if np.isfinite(log_prior(theta = theta)):
        ln_post = log_prior(theta = theta) + log_likelihood(theta = theta, time_values = time_values, 
                                    f_obs = f_obs, f_err = f_err, var = var, is_upper_limit = is_upper_limit)
    else:
        ln_post = -np.inf
    return ln_post

    
def dummy_log_posterior(theta):
    '''
    with intrinsic var
    '''
    lnP = log_posterior(theta = theta, time_values = scaled_time[times_to_fit], f_obs = flux_obs[times_to_fit], 
                        f_err = flux_err[times_to_fit], var = var[times_to_fit], is_upper_limit = is_upper_limit[times_to_fit])
    return lnP


'''
Get lightcurve target info for MCMC
'''

jd, filt, mag, sigma_mag, limmag, instrument = read_in_lc(target_lc_direc)
is_upper_limit, flux_obs, flux_err, var_mag, var = calibrate_target(mag, sigma_mag, limmag)
time_zero, scaled_time, scaled_t_d1, times_to_fit = scale_target_times(jd, is_upper_limit)


# In[28]:


'''
MCMC
'''

nwalkers = 64
ndim = 3
ball_p0 = 0.2*np.random.randn(3, 64).T + np.array([1, 12.5, 3])

nburn = 200
nstep = 1000

sampler = emcee.EnsembleSampler(nwalkers, ndim, dummy_log_posterior)
pos, prob, state = sampler.run_mcmc(ball_p0, nburn)
sampler.reset()

for i, result in enumerate(sampler.sample(pos, iterations=nstep)):
    if (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nstep))


# In[29]:


'''
Plotting t0
'''
f = plt.figure(figsize= (15,10))

gs = GridSpec(10,1)
ax = f.add_subplot(gs[:8,0])
ax2 = f.add_subplot(gs[8:,0])
gs.update(hspace=0)

ax.errorbar(scaled_time[~is_upper_limit], flux_obs[~is_upper_limit],  yerr = flux_err[~is_upper_limit], fmt='r.', 
            label = r"$\sigma_{\mathrm{comb}}$")
ax.errorbar(scaled_time[(~is_upper_limit) & times_to_fit], flux_obs[(~is_upper_limit) & times_to_fit],  fmt='.', 
            color = 'lightsalmon', marker = 'o', mec = 'k', ms = 10)
ax.errorbar(scaled_time[is_upper_limit], flux_obs[is_upper_limit], fmt='k.', label = 'non-detection')


time_grid = np.linspace(0, 30, 1000)

p0 = np.median(sampler.flatchain, axis=0)

Nsamples = 1000
samples = sampler.flatchain[np.random.randint(0, len(sampler.flatchain), Nsamples)]


med, upper, lower = [], [], []
med1sig, upper1sig, lower1sig = [], [], []

pred_fluxes = []
for sample in samples:
    pred_fluxes.append(calc_flux(time_grid, sample[0], sample[1], sample[2]))
    
pred_fluxes = np.array(pred_fluxes)

for i, t_ in enumerate(time_grid):
    med_, upper_, lower_ = np.percentile(pred_fluxes.T[i], [50, 97.5, 2.5])
    med1sig_, upper1sig_, lower1sig_ = np.percentile(pred_fluxes.T[i], [50, 84.1, 15.9])
    med.append(med_)
    upper.append(upper_)
    lower.append(lower_)
    med1sig.append(med1sig_)
    upper1sig.append(upper1sig_)
    lower1sig.append(lower1sig_) 

ax.set_title(('%s'%(ZTF_name)), loc='left', color = 'red', fontsize = '30')


ax.set_ylabel('Normalized flux', fontsize = '20')
ax.set_xlabel('Time [d]', fontsize = '20')

ax.fill_between(time_grid,  lower, upper, alpha = 0.5, color = 'lightsteelblue', label = '95%')
ax.fill_between(time_grid,  lower1sig, upper1sig, alpha = 0.2, color = 'navy', label = '68%')
ax.plot(time_grid, med, alpha = 1, color = 'k', label = 'curve fit')

ax.legend(loc = 2, prop={'size': 17})

'''
extract values for plotting purposes
'''

samples[:, 2] = np.exp(samples[:, 2])
# sampler.flatchain[:, 2] = np.exp(sampler.flatchain[:, 2])
fm_mcmc, t0_mcmc, te_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(sampler.flatchain, [16, 50, 84],
                                                axis=0)))

t0_rd = round(t0_mcmc[0], 2)
t0_exp_rd = round(t0_mcmc[1], 2)
t0_sub_rd = round(t0_mcmc[2], 2)

ax.set_title(f'$t_0 = {t0_rd}^{{{t0_exp_rd}}}_{{{t0_sub_rd}}}$', fontsize = 30)

time_zero_rd = round(time_zero, 2)
set_time = str(time_zero_rd)
ax.set_title(f'MJD$= t + {set_time[2:]}$', loc = 'right', fontsize = '20')

# tick info
ax.xaxis.set_ticklabels([])
ax.set_xticks([])
ax2.yaxis.set_ticklabels([])
ax2.set_yticks([])
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.tick_params(direction='in', length = 6)
ax2.tick_params(direction='in', length = 6)

#PDF

ax2.set_ylabel('$P(t_{0})$', fontsize = '20')
ax2.set_xlabel('Time [d]', fontsize = '20')

med_hist, lower_hist, upper_hist = np.percentile(sampler.flatchain.T[1], [50, 2.5, 97.5])
med1sig_hist, lower1sig_hist, upper1sig_hist = np.percentile(sampler.flatchain.T[1], [50, 15.9, 84.1])

lower_hist_t = sampler.flatchain.T[1].flat[np.abs(sampler.flatchain.T[1] - lower_hist).argmin()]
upper_hist_t = sampler.flatchain.T[1].flat[np.abs(sampler.flatchain.T[1] - upper_hist).argmin()]
lower1sig_hist_t = sampler.flatchain.T[1].flat[np.abs(sampler.flatchain.T[1] - lower1sig_hist).argmin()]
upper1sig_hist_t = sampler.flatchain.T[1].flat[np.abs(sampler.flatchain.T[1] - upper1sig_hist).argmin()]

N, bins, _ = ax2.hist(sampler.flatchain[:,1], bins=100, density=True, 
                   histtype='step', lw=2, color = 'black')

bincenters   = 0.5*(bins[1:]+bins[:-1])  

ax2.fill_between(bincenters, 0, N, interpolate=True,
                where=((bincenters>=lower_hist_t) &
                       (bincenters<=upper_hist_t)), alpha = 0.5, color = 'lightsteelblue')

ax2.fill_between(bincenters, 0, N, interpolate=True,
                where=((bincenters>=lower1sig_hist_t) &
                       (bincenters<=upper1sig_hist_t)), alpha = 0.2, color = 'navy')


ax2.set_ylabel('$P(t_{0})$', fontsize = '20')
ax2.set_xlabel('Time [d]', fontsize = '20')

    
ax.set_xlim(0, 30)
ax2.set_xlim(0, 30)
ax.set_ylim(0, 1.2)
ax2.set_ylim(1.5,0)

plt.xticks(fontsize ='20')
plt.yticks(fontsize ='20')

# f.savefig(f'/global/cscratch1/sd/sarafina/ztf_targets_IIp/lc_fits/{ZTF_name}_explosion_date')


# In[30]:


'''
Rejection sampling from hist
Random sampling from histogram of data by calculating its cumulative distribution, extracting a random number from a uniform distribution between 0 and 1 
and finally find the index of the cumulative distribution where the random number should be inserted to keep the cumulative distribution sorted. 
The sampled number from the histogram will be the histogram bin corresponding to this histogram.
'''

hist, _ = np.histogram(sampler.flatchain[:,1], bins=100)
bin_midpoints = (bins[:-1] + bins[1:])/2
cdf = np.cumsum(hist)
cdf = cdf / cdf[-1]

values = np.random.rand(1000)
value_bins = np.searchsorted(cdf, values)
random_from_cdf = bin_midpoints[value_bins]
explosion_dates = list(random_from_cdf + time_zero - 2400000)