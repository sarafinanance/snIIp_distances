import specutils
from specutils import analysis, Spectrum1D


targets_spectra_direc = f'/Users/sarafinanance/Desktop/ztf_IIp_targets/spectra/{ZTF_name}'
target_lc_direc = f'/Users/sarafinanance/Desktop/ztf_IIp_targets/lightcurves/{ZTF_name}'
filters_dir = f'/Users/sarafinanance/Desktop/ztf_IIp_targets/{ZTF_name}/color_lcs'
day50_file_dir = f'/Users/sarafinanance/Desktop/ztf_IIp_targets/spectra/ZTF19abqhobb'

files_lc = os.listdir(target_lc_direc)
files_spectra = os.listdir(targets_spectra_direc)

templates_spectra_dir_fe = '/Users/sarafinanance/Desktop/ztf_IIp_targets/templates_spectra/fe'
templates_spectra_dir_hbeta = '/Users/sarafinanance/Desktop/ztf_IIp_targets/templates_spectra/hbeta'

# conversions
fe_file = '/Users/sarafinanance/Desktop/ztf_IIp_targets/conversions/fe.dat'
hbeta_file = '/Users/sarafinanance/Desktop/ztf_IIp_targets/conversions/hbeta.dat'
fe_vels_file = '/Users/sarafinanance/Desktop/ztf_IIp_targets/conversions/fe_vel.dat'
hbeta_vels_file = '/Users/sarafinanance/Desktop/ztf_IIp_targets/conversions/hbeta_vel.dat'
ratios_file = '/Users/sarafinanance/Desktop/ztf_IIp_targets/conversions/ratios.dat'

ratios = pd.read_csv(ratios_file, delim_whitespace=True, names = ['epoch', 'ratio', 'sig_ratio'])
fe_vels = pd.read_csv(fe_vels_file, delim_whitespace=True, names = ['template', 'velocity'])
hbeta_vels = pd.read_csv(hbeta_vels_file, delim_whitespace=True, names = ['template', 'velocity'])
vels = fe_vels.merge(hbeta_vels, how = 'outer')


target_file_names = []
for file in files_spectra:
    # if file.endswith(".txt"):
    target_file_names.append(file)
    

def make_Table_target(file):
    '''
    Create astropy Table of the target night spectrum 
    '''
    try:
        colnames_target = ['wavelength', 'flux', 'sig']
        df_target_ = Table.read(os.path.join(targets_spectra_direc, file), format = 'ascii',  comment='#', names = colnames_target)
        df_target_['sig'] = np.sqrt(df_target_['sig'])
    except:        
        colnames_target = ['wavelength', 'flux']
        df_target_ = Table.read(os.path.join(targets_spectra_direc, file), format = 'ascii',  comment='#', names = colnames_target)
        df_target_['sig'] = df_target_['flux'] * 0.05
        
    target_night_ = file.split('_')[1]   
    mjd_ = Time('{}-{}-{}'.format(target_night_[:4], target_night_[4:6], target_night_[6:])).mjd
    return df_target_, target_night_, mjd_

"""
spectra template imports 
"""

def import_template(file, template_dir):
    colnames_template = ['wavelength', 'flux']
    file_loc = os.path.join(template_dir, file)
    name = os.path.splitext(file)[0]
    df_template = Table.read(file_loc, format = 'ascii', comment = '#', names = colnames_template)
    return df_template, name


# In[23]:


"""
functions for spectra
"""

def format_target(target):
    """Returns formatted target
    
    Parameters
    ----------
    target : Table (astropy)
        Target info [wavelength, flux, sig flux]

    Returns
    -------
    concat_table : Table (astropy)
        Reformatted target by redshifting, constraining wavelength range, normalizing and smoothing
    """
    minwave = 4000.0
    maxwave = 6000.0
    
    target['wavelength'] = target['wavelength'] / (1 + z_input)
    target_min = target[target['wavelength'] > minwave]
    concat_target = target_min[target_min['wavelength'] < maxwave]
    
    p_init = models.Polynomial1D(degree=3)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(p_init, concat_target['wavelength'], concat_target['flux'])
    
    flattened_flux = concat_target['flux']/g(concat_target['wavelength'])
    flattened_sig = concat_target['sig']/g(concat_target['wavelength'])
    
    concat_target['flat_flux'] = flattened_flux
    concat_target['flat_sig'] = flattened_sig

    concat_target.columns['wavelength', 'flat_flux', 'flat_sig']
    
    return concat_target

def format_template(template):
    """Returns formatted template
    
    Parameters
    ----------
    template : Table (astropy)
        Template info [wavelength, flux, sig flux]

    Returns
    -------
    concat_template : df
        Reformatted template by constraining wavelength range, normalizing and smoothing
    """
    minwave = 4000.0
    maxwave = 6000.0
    
    template_min = template[template['wavelength'] > minwave] 
    concat_template = template_min[template_min['wavelength'] < maxwave]
    concat_template = concat_template.copy()
    
    p_init = models.Polynomial1D(degree=3)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(p_init, concat_template['wavelength'], concat_template['flux'])
    
    flattened_flux = concat_template['flux']/g(concat_template['wavelength'])
    flattened_sig = (flattened_flux * 0.05)/g(concat_template['wavelength'])
    
    concat_template['flux'] = flattened_flux
    concat_template['flat_sig'] = flattened_sig
    
    concat_template.columns['wavelength', 'flux', 'flat_sig']
        
    return concat_template


def make_spec1d_target(target):
    """Returns 1D spectrum of target
    
    Parameters
    ----------
    target : Table (astropy)
        Formatted target

    Returns
    -------
    p_obs : Target spectrum container for 1D spectral data
    """    
    spec_unit = u.MJy
    wave = target['wavelength']
    flux = target['flat_flux']
    sig = target['flat_sig']
    
    spec1d_target = Spectrum1D(spectral_axis=wave* u.AA, flux=flux * spec_unit, uncertainty=StdDevUncertainty(sig, unit='Jy'))
    
    g1_fit = fit_generic_continuum(spec1d_target)
    y_continuum_fitted = g1_fit(wave*u.AA)
    
    wht_obs = 1 / sig**2

    dataspec = QTable([wave*u.AA, (spec1d_target.flux - y_continuum_fitted)*spec_unit, wht_obs, sig],
                       names=('wavelength','flux', 'weight', 'uncertainty')) #QTable is a class for heterogeneous tabular data

    dataspec_sub = dataspec[dataspec['weight'] > 0.]

    p_obs = Spectrum1D(spectral_axis=dataspec_sub['wavelength'], flux=dataspec_sub['flux'] * spec_unit, uncertainty=StdDevUncertainty(dataspec_sub['uncertainty']), rest_value = 6000 * u.AA, unit='Jy')
    
    return p_obs

def make_spec1d_temps(template):
    """Returns dictionary of 1D spectra for all templates
    
    Parameters
    ----------
    key : list
        Dictionary keys of names of formatted templates
    values : list
        Dictionary values of flux and wavelngth of formatted templates
    Returns
    -------
    template_spec1d : spec1d container
        dict of spectrum containers for 1D spectral data for all templates
    """     
    spec_unit = u.MJy
    
    wavelength = template['wavelength']
    flux = template['flux']
    sig = template['flat_sig']
    wht_obs = 1 / sig**2
    
    spec1d_template = Spectrum1D(spectral_axis=wavelength* u.AA, flux=flux * u.MJy, uncertainty=StdDevUncertainty(sig, unit='Jy'))
    
    g1_fit = fit_generic_continuum(spec1d_template)
    y_continuum_fitted = g1_fit(wavelength*u.AA)
    
    dataspec = QTable([wavelength*u.AA, (spec1d_template.flux-y_continuum_fitted)*spec_unit, wht_obs, sig],
                       names=('wavelength','flux', 'weight', 'uncertainty')) #QTable is a class for heterogeneous tabular data
    
    dataspec_sub = dataspec[dataspec['weight'] > 0.]
    
    template_spec1d = Spectrum1D(spectral_axis=dataspec_sub['wavelength'], flux=dataspec_sub['flux'] * spec_unit, uncertainty=StdDevUncertainty(dataspec_sub['uncertainty']*np.ones(flux.shape)), unit='Jy')
    
    return template_spec1d 
    
targets_tables = []
target_dates = []
mjds = []

for file in target_file_names:
    """
    Loop through all target spectra to make lists of dfs, dates, and mjds
    """
    df_target, target_date, mjd = make_Table_target(file)
    targets_tables.append(df_target)
    target_dates.append(target_date)
    mjds.append(mjd)
    
def calc_day_post_expl(mjd, explosion_date):
    day_post_explosion = mjd - explosion_date
    return day_post_explosion

days_post_explosion = []
"""
For each target night, calculate days post explosion from the 1000 explosion dates calculated above.
Since we have 4 target nights, this is a list of 4000 days post explosion.
"""
for mjd in mjds:
    for explosion_date in explosion_dates[:1]:
        day_post_explosion = calc_day_post_expl(mjd, explosion_date)
        days_post_explosion.append(day_post_explosion)


# In[32]:


t = dict(zip(mjds, targets_tables))
"""
Make df of dict
"""
s_t = pd.Series(t)
df_t = s_t.to_frame()
df_t.reset_index(inplace=True)
df_t = df_t.rename(columns = {'index': 'mjds', '0' : 'astropy table'})
df_mjds = df_t.rename(columns ={ df_t.columns[1] : 'mjds_orig_spectra'})


# In[33]:


"""
format over pandas df
"""
df_mjds['formatted_mjds'] = df_mjds['mjds_orig_spectra'].map(lambda row: format_target(row))
df_mjds['mjd_spec1d_obj'] = df_mjds['formatted_mjds'].map(lambda row: make_spec1d_target(row))
df_mjds


# In[38]:


"""
Import template spectra
"""

fe_temp_lst = []
hbeta_temp_lst = []

for file in os.listdir(templates_spectra_dir_fe):
    if file.endswith(".dat"):
        df_template_fe = import_template(file, templates_spectra_dir_fe)
        fe_temp_lst.append(df_template_fe)
        
for file in os.listdir(templates_spectra_dir_hbeta):
    if file.endswith(".dat"):
        df_template_hbeta = import_template(file, templates_spectra_dir_hbeta)
        hbeta_temp_lst.append(df_template_hbeta)
        
"""
Format Templates
"""

temp_fe_lst_format = []
temp_hbeta_lst_format = []

#fe temps

for template, name in fe_temp_lst:
    temp_fe_format = format_template(template)
    temp_fe_lst_format.append(temp_fe_format)
    
#hbeta temps


for template, name in hbeta_temp_lst:
    temp_hbeta_format = format_template(template)
    temp_hbeta_lst_format.append(temp_hbeta_format)
    
    
"""
get names of all templates
"""

temp_names_fe = []
temp_names_hbeta = []

for template, name in fe_temp_lst:
    temp_names_fe.append(name)
    
for template, name in hbeta_temp_lst:
    temp_names_hbeta.append(name)
    
temp_names = temp_names_fe + temp_names_hbeta

"""
make spec1d template object
"""

spec1d_temps_dict_fe = []
spec1d_temps_dict_hbeta = []

#fe
for template in temp_fe_lst_format:
    template_1d_fe = make_spec1d_temps(template)
    spec1d_temps_dict_fe.append(template_1d_fe)
    
#hbeta
for template in temp_hbeta_lst_format:
    template_1d_hbeta = make_spec1d_temps(template)
    spec1d_temps_dict_hbeta.append(template_1d_hbeta)