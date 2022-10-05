"""
FOR TOMORROW:

STEP 3. 
This spectrum was then de-redshifted and the restframe V-I color of the SN was calculated 
De-redshift the spetrum and calculate restframe V-I color

- The V and I apparent magnitudes 
- M_I (absolute) 
- mu = I - M_I 
--> - From M_I, subtract off I from sncosmo

-> SETUP MCMC FOR UNCERTAINTIES (SEE BELOW)
-> CALCULATE M_I_ABSOLUTE FROM V_FE_50 & V and I colors above
M_I = -4.6 log10(vFe/5000)+0.7 ((V-I)-0.53) - 17.43
SEE BELOW FUNCTION, JUST NEED TO PLUG IN V_ AND I_
-> CALCULATE MU = I - M_I

QUESTIONS FOR PETER:
- HOW TO FOLD IN MULTIPLE EXPLOSION DATES FOR V50 CALCULATION (EASY CODE QUESTION)?
"""


# In[ ]:


#calculate v(t)

bounds = ((0,-0.481),(10000, -0.447))
p0_powerlaw = [5000, -0.464] #km/sec and slope
alpha = 5.81
M_I0 = -17.52
V_I_0 = 0.53
V_I = V_rest - I_rest
p0 = [5000, -0.464] #km/sec and slope

def v_powerlaw(t, v, m):
    """
    Nugent17
    t = independent variable
    v_50 = velocity at day=50
    m = slope
    """
    v50 = v * (t/50)**(m)
    return(v50)


"""
NEXT: need to double check whether this should be days_post_explosion or explosion_date, because this is throwing it off. 
Think it needs to be days_post_explosion

Need to check what to do about uncertainties
-> Fit to g', r', i' each have an uncertainty 
---> Do another MCMC where the uncertainties are g'r'i', and then you feed those in, dispersed by their uncertainties, to sncosmo to get new 
---> V, I
---> Run 1000 MCMCs to see what the uncertainties are on each band which will impact V and I uncertainties 

-> Fit to velocity has an uncertainty (from the covar) - might be small compared to just the method itself 
-> Question is how to fold all of these uncertainties in! 
---> Easiest way: carry them along, take derivative of M_I
---> d(M_I)/d(v) where v is the velocity. 
---> d(M_I)/d(V) where V is mag, and also for I
---> and just from that equation, can give final uncertainty on velocity and all mags is this value
"""

# popts = [] #velocities
# covars = [] #variances (sig^2)

# for day in days_post_explosion:
#     popt, covar = curve_fit(v_powerlaw, day, vels_fit, p0_powerlaw)#, sigma = final_errs, absolute_sigma = True, bounds = bounds)
#     popts.append(popt)
#     covars.append(covar)

# time_arr = np.linspace(np.min(explosion_dates), np.max(explosion_dates), num=1000)

# calculated_vels_arr = v_powerlaw(time_arr, popt[0], popt[1])


# In[ ]:


v50, covar = curve_fit(v_powerlaw, days_post_explosion, final_vels, p0_powerlaw)
print(f'v50 = {v50[0]} km/sec')
time_arr = np.linspace(np.min(days_post_explosion), np.max(days_post_explosion), 10)
calculated_vels_arr = v_powerlaw(time_arr, v50[0], v50[1])
v50_sig = np.sqrt(covar[0][0])


# In[ ]:


days_post_explosion_all_lst = []
"""
For each target night, calculate days post explosion from the 1000 explosion dates calculated above.
Since we have 4 target nights, this is a list of 4000 days post explosion.
"""
for mjd in mjds:
    mjd_lst = []
    for explosion_date in explosion_dates:
        day_post_explosion = calc_day_post_expl(mjd, explosion_date)
        mjd_lst.append(day_post_explosion)
    days_post_explosion_all_lst.append(mjd_lst)


# In[ ]:


plt.figure(figsize = (10, 6))
for i in range(len(df_mjds.mjds.values)):
    plt.scatter(days_post_explosion[i], final_vels[i], label = df_mjds.mjds.values[i], zorder = 6, edgecolor = 'k', s = 80)
    plt.errorbar(days_post_explosion[i], final_vels[i], final_sigs[i], marker = 'o', ls = '')

plt.plot(time_arr, calculated_vels_arr)
plt.scatter(50, v50[0], marker = 'o', zorder = 6, s = 80, edgecolor = 'k', label = f'v_50 = {round(v50[0], 2)}')
plt.errorbar(50, v50[0], yerr = v50_sig, fmt='r.')
plt.ylabel('velocities (km/s)')
plt.xlabel('days post explosion')
plt.grid()
plt.minorticks_on()
plt.title('%s'%(ZTF_name[:12]))
plt.legend(loc = 0)
# plt.savefig('/global/cscratch1/sd/sarafina/ztf_targets_IIp/lc_fits/ZTF19abqhobb_vels.png')


# In[ ]:


"""
From Graham et al:
They find an intrinsic 𝐼-band magnitude
of −17.84 ± 0.14 mag at 50 days, a distance modulus of 𝜇 = 34.95 ±
0.26 mag and a distance of 𝐷 = 97.6 ± 12 Mpc for SN 2019nvm.
"""

def calculate_M_I(v):
    """
    Nugent17
    t = independent variable
    v_50 = velocity at day=50
    m = slope
    (V-I)0 = 0.53 from paper #base color
    """
    M_I = -alpha * np.log10(v/5000) - 1.36*((V_I) - 0.53) + M_I0
    return M_I

def calculate_mu(I, M_I):
    mu = I - M_I
    return mu

M_I = calculate_M_I(v50[0])
print(f'I-band Mag: {round(M_I, 2)}')
mu = calculate_mu(I_rest, M_I)
print(f'Distance modulus to {ZTF_name}: {round(mu, 3)}')