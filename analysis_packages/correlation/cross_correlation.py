"""
CROSS-CORRELATION & LAG CALCULATIONS
"""
wfe = 5169.0 #A
whb = 4861.0 #A

def calculate_lag(target_spec1d, template_spec1d):
    """
    Compute cross-correlation and extract max correlation/lag
    De-redshift spectral axis by dividing by associated lag
    Creates dictionary of lags, redshifts, and de-redshifted wavelengths. 
    
    Parameters
    ----------
    p_obs : Target spectrum container for 1D spectral data
    template : list
         dict of fluxes and wavelengths for formatted templates
    #values : list
        #Table of correlation and lag values
        
    Returns
    -------
    df_vals_target : dict
        Dictionary of lags, z, and z_blue for each target
    """   
    #calculate lag
    
    corr, lag = correlation.template_correlate(target_spec1d, template_spec1d) #, apodization_window = 0.5, resample = True) #compute cross-correlation

    df = pd.DataFrame({'corr':corr, 'lag':lag})
    argmax = df['corr'].idxmax()
    max_lag = df.lag[argmax] # find max lag
    
    return max_lag

def delag_template(max_lag, template_spec1d):
    
    z = max_lag / const.c.to('km/s') # compute redshift from lag
    delag_const = 1 - z.value
    spectral_axis_template = template_spec1d.spectral_axis.value # define template spectral axis
    delag_specaxis = spectral_axis_template / delag_const # de-redshift template spectral axis 
    
    #de-redshifted template
    template_lag_spec1d = Spectrum1D(spectral_axis=delag_specaxis * u.AA, flux=template_spec1d.flux, uncertainty=template_spec1d.uncertainty)
    
    return template_lag_spec1d


# In[24]:


"""
CHI2 FUNCTIONS
"""
_KMS = u.Unit('km/s')


def _normalize_for_template_matching(observed_spectrum, template_spectrum, stddev=None):
    
    if stddev is None:
        stddev = _uncertainty_to_standard_deviation(observed_spectrum.uncertainty)
        
    num = np.sum((observed_spectrum.flux*template_spectrum.flux) / (stddev**2))
    denom = np.sum((template_spectrum.flux / stddev)**2)

    return num/denom

def _uncertainty_to_standard_deviation(uncertainty):
    """
    Convenience function to convert other uncertainty types to standard deviation,
    for consistency in calculations elsewhere.

    Parameters
    ----------
    uncertainty : :class:`~astropy.nddata.NDUncertainty`
        The input uncertainty

    Returns
    -------
    :class:`~numpy.ndarray`
        The array of standard deviation values.

    """
    if uncertainty is not None:
        if isinstance(uncertainty, StdDevUncertainty):
            stddev = uncertainty.array
        elif isinstance(uncertainty, VarianceUncertainty):
            stddev = np.sqrt(uncertainty.array)
        elif isinstance(uncertainty, InverseVariance):
            stddev = 1 / np.sqrt(uncertainty.array)

        return stddev

def _chi_square_for_templates(observed_spectrum, template_spectrum, resample_method):
    """
    Resample the template spectrum to match the wavelength of the observed
    spectrum. Then, calculate chi2 on the flux of the two spectra.

    Parameters
    ----------
    observed_spectrum : :class:`~specutils.Spectrum1D`
        The observed spectrum.
    template_spectrum : :class:`~specutils.Spectrum1D`
        The template spectrum, which will be resampled to match the wavelength
        of the observed spectrum.

    Returns
    -------
    normalized_template_spectrum : :class:`~specutils.Spectrum1D`
        The normalized spectrum template.
    chi2 : `float`
        The chi2 of the flux of the observed spectrum and the flux of the
        normalized template spectrum.
    """
    # Resample template
    fluxc_resample = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')
    template_obswavelength = fluxc_resample(template_spectrum, observed_spectrum.spectral_axis) #this resamples the template's uncertainty too

    # Convert the uncertainty to standard deviation if needed
    stddev_obs = _uncertainty_to_standard_deviation(observed_spectrum.uncertainty)
    stddev_temp = _uncertainty_to_standard_deviation(template_obswavelength.uncertainty)
    
    # Normalize spectra
    normalization = _normalize_for_template_matching(observed_spectrum, template_obswavelength, stddev_obs)

    # Numerator
    num_right = normalization * template_obswavelength.flux
    num = observed_spectrum.flux - num_right

    # Denominator
    # denom = stddev * observed_spectrum.flux.unit
    denom_sq = ((stddev_obs * observed_spectrum.flux.unit)**2 + (stddev_temp * template_obswavelength.flux.unit)**2) #denom including template uncertainty

    # Get chi square
    result = num**2/denom_sq
    chi2 = np.sum(result.value)
    
    # Create normalized template spectrum, which will be returned with corresponding chi2
    normalized_template_spectrum = Spectrum1D(spectral_axis=template_spectrum.spectral_axis, flux=template_spectrum.flux*normalization)

    return normalized_template_spectrum, chi2


def _reduced_chi2(normalized_template_spectrum, chi2):
    """
    redchi2 = chi2/(len(data)-dof)
    here, dof = 1
    
    Parameters
    ----------
    template :class:`~specutils.Spectrum1D` or :class:`~specutils.SpectrumCollection` or `list`
    chi2 : `float`
        The chi2 of the flux of the observed_spectrum and the flux of the normalized template spectrum.

    Returns
    -------
    reuced_chi2 : `float`
        The reduced chi2
    """
    n_lam = len(normalized_template_spectrum.spectral_axis)
    reduced_chi2 = chi2/(n_lam - 1)
    
    return reduced_chi2

def template_match(target, temp_name, observed_spectrum, spectral_templates, resample_method="linear_interpolated", redshift=None):

    final_redshift = None

    if hasattr(spectral_templates, 'flux') and len(spectral_templates.flux.shape) == 1:
        
        normalized_spectral_template, chi2 = _chi_square_for_templates(observed_spectrum, spectral_templates, resample_method)
        reduced_chi2 = _reduced_chi2(normalized_spectral_template, chi2)
        
        # print(f'{target} x {temp_name} has {reduced_chi2}')
        
        return pd.Series([normalized_spectral_template, reduced_chi2]) 
        """
PUT the Fe & Hbeta templates together in one list.
"""
merged_temps = spec1d_temps_dict_fe + spec1d_temps_dict_hbeta

template_data_keys_and_vals = {'template': temp_names, 'template_data_orig': merged_temps}
template_data = pd.DataFrame(template_data_keys_and_vals)

template_data_w_vels = pd.merge(template_data, vels)

"""
Classify templates as hbeta or fe
"""

merged_temps = pd.merge(template_data_w_vels, fe_vels, on = ['template'], how = 'left', indicator = 'Fe')
merged_temps.drop('velocity_y', inplace = True, axis = 1)
merged_temps['Fe'] = np.where(merged_temps.Fe == 'both', True, False)


# In[40]:


colnames = ['target', 'template', 'target_spectrum', 'Fe', 'lag', 'lag_shifted_temp']#, 'red_chi2']
lst =[]
for row_num, vals in df_mjds.iterrows():
    for row_num_temp, vals_temp in merged_temps.iterrows():
        lag = calculate_lag(vals.mjd_spec1d_obj, vals_temp.template_data_orig)
        lag_shifted_temp = delag_template(lag, vals_temp['template_data_orig'])
        
        #compute chi2 between target and lag shifted template
        lst.append([vals.mjds, vals_temp.template, vals.mjd_spec1d_obj, vals_temp.Fe, lag, lag_shifted_temp])

df_lags = pd.DataFrame(lst, columns = colnames)


# In[41]:


"""
does the template normalization in chi2 calculation have to take into account the template uncertainty? 
"""

df_lags[['normalized_spectral_template', 'red_chi2']]= df_lags.apply(lambda x: template_match(x['target'], x['template'], x['target_spectrum'], x['lag_shifted_temp']), axis = 1)


# In[42]:


df_mjds.sort_values(by = ['mjds'])


# In[49]:


# print(f'day post explosion: {round(days_post_explosion[1],1)}:\n')
# df_lags[df_lags['target'] == 58665.0].sort_values(by=['red_chi2']).head(3)


# In[50]:


# plt.figure(figsize = (10, 5))
# plt.plot(df_lags.lag_shifted_temp[25].spectral_axis, df_lags.lag_shifted_temp[25].flux, label = f'best fit {df_lags.template[0]}')
# plt.plot(df_lags.target_spectrum[25].spectral_axis,df_lags.target_spectrum[25].flux, label = f'target {df_lags.target[0]}')
# plt.plot(df_lags.lag_shifted_temp[0].spectral_axis,df_lags.lag_shifted_temp[0].flux, label = f'3rd {df_lags.template[0]}')

# plt.plot(df_mjds.mjd_spec1d_obj[2].spectral_axis, df_mjds.mjd_spec1d_obj[2].flux, label = f'{df_mjds.mjds[2]}')
# plt.legend()


# In[52]:


# print(f'day post explosion: {round(days_post_explosion[1],1)}:\n')
# df_lags[df_lags['target'] == df_mjds.mjds[1]].sort_values(by=['red_chi2']).head(3)


# In[ ]:


# print(f'day post explosion: {round(days_post_explosion[1],1)}:\n')
# df_lags[df_lags['target'] == 58789.0].sort_values(by=['red_chi2']).head(3)


# In[54]:


def correct_vels(velocity, lag):
    lag_shifted_vel = velocity - lag
    return lag_shifted_vel

unique_names = df_lags.target.unique()
dfLags_dict = {elem : pd.DataFrame for elem in unique_names}

for key in dfLags_dict.keys():
    dfLags_dict[key] = df_lags[:][df_lags.target == key]

df_lags['velocity'] = pd.concat([merged_temps['velocity_x']] * 4, ignore_index = True)

for vals, target in dfLags_dict.items():
    lag_shifted_vels = correct_vels(df_lags["velocity"], target['lag'])
    df_lags[f'{vals}_lagshift_vel'] = lag_shifted_vels


# In[55]:


ratio_lst = [1.09, 1.09, 1.10, 1.11, 1.12, 1.12, 1.13, 1.14, 1.14, 1.15, 1.16, 1.16, 1.17, 1.17, 1.17, 1.18, 1.18, 1.18, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19]
ratio_sig_lst = [0.09, 0.09, 0.08, 0.07, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.07, 0.07, 0.08, 0.08]
epoch_lst = list(range(5, 31))


# In[56]:


cols_vhb = ['epoch','ratio','ratio_sig']
vhb_ratio_keys_and_vals = {'epoch': epoch_lst, 'ratio': ratio_lst, 'ratio_sig': ratio_sig_lst}
vhb_ratio_df = pd.DataFrame(vhb_ratio_keys_and_vals) 


# In[57]:


df_mjds['days_post_explosion'] = days_post_explosion
df_mjds['Fe'] = ['False' if x < 30.0 else 'True' for x in df_mjds['days_post_explosion']]


# In[58]:


df_mjds['epoch'] = [round(x) if x < 30.0 else int(30.0) for x in df_mjds['days_post_explosion']]


# In[59]:


df_mjds


# In[62]:


for idx, row in df_mjds.iterrows():
    epoch = row['epoch']
    ratio_match = vhb_ratio_df[vhb_ratio_df['epoch'] == epoch]['ratio']
    sig_match = vhb_ratio_df[vhb_ratio_df['epoch'] == epoch]['ratio_sig']
    mjd = row.loc['mjds']
    rowIndex = df_mjds.index[df_mjds['mjds'] == mjd]
    df_mjds.loc[rowIndex, 'vhb_ratio'] = ratio_match.values
    df_mjds.loc[rowIndex, 'vhb_ratio_sig'] = sig_match.values


# In[63]:


for col_name, template in df_lags.iteritems():
    if col_name in df_lags.filter(like='lagshift'):
        col_name_strip = col_name.strip('_lagshift_vel')
        for idx, row in df_mjds.iterrows():
            if col_name_strip == str(row['mjds']):
                final_vel = df_lags[col_name]/row['vhb_ratio']
                final_vel_sig = row['vhb_ratio_sig']
                mjd_night = col_name.replace('_lagshift_vel','')
                df_lags.loc[df_lags['Fe'] == True, f'{mjd_night}_final_vel'] = df_lags[col_name]
                df_lags.loc[df_lags['Fe'] == False, f'{mjd_night}_final_vel'] = final_vel
                df_lags[f'{mjd_night}_final_vel_sig'] = final_vel_sig


# In[64]:


"""
CALCULATE FINAL VELOCITIES
"""
def calculate_avg_final_vel(target):
    """
    Compute final velocity by averaging the final 3 velocities of the lowest reduced chi2s
    Parameters
    ----------
    target : str
        target name    
    Returns
    -------
    target : str
        target name
    vf : float
        averaged final velocity
    """   
    v1, v2, v3 = df_lags[df_lags['target'] == target].sort_values(by=['red_chi2'])[f'{target}_final_vel'].head(3)
    lowest_vel_lst = [v1, v2, v3]
    vf = df_lags[df_lags['target'] == target].sort_values(by=['red_chi2'])[f'{target}_final_vel'].head(3).mean()
    return target, vf, lowest_vel_lst

"""
CALCULATE UNCERTAINTIES

Generic baseline uncertainty minimum: 250 km/sec 
If the following uncertainty is larger, then put that
Got the best 3 chi2. Take the highest velocity - lowest velocity out of those 3. If that is > 250, use that. If < 250, use 250. 

Explanation of 250 km/sec:
Typical speed of spiral arms is 250 km/sec. Can't tell which arm the SN is in because of line of sight, so using 250 km/sec for a typical unceratinty
When we measure galaxy redshift, you put a fiber on the core and measure it. that's how we get the galaxy's redshift
Conservatively say that it's 250 km/sec 
This should show up naturally in the dispersion fo the top 3 fits. 
"""

def calculate_uncertainty(lowest_vel_lst):
    vmax = max(lowest_vel_lst)
    vmin = min(lowest_vel_lst)
    v_diff = vmax - vmin 
    if v_diff > 250.0:
        return v_diff
    else:
        return 250.0               
                          
d_vf = {}
for idx, vals in df_mjds.iterrows():
    target_name, vf, lowest_vels = calculate_avg_final_vel(vals.mjds)
    vsig = calculate_uncertainty(lowest_vels)
    if target_name not in d_vf:
        d_vf[target_name] = []
    pair = vf, vsig
    d_vf[target_name].append(pair)

d_vf_lst = [[k, *v] for k, lst in d_vf.items() for v in lst]
df_vels = pd.DataFrame(d_vf_lst, columns = ['mjd', 'avg_vel', 'sig'])
df_final = df_mjds.join(df_vels.drop('mjd', axis = 1))
final_vels = df_final.avg_vel.values
final_sigs = df_final.sig.values
dates_fit = df_final.mjds.values