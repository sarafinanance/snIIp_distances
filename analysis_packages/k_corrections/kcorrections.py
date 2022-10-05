
"""
For every t0 [explosion_date] (of which there are 1,000), calculate:
- A velocity at day=50

# To calculate the restframe V-I colors during the plateau, we estimate explosion date and then linear interpolate the g'r'i'z colors at day 50(1+z). To each of the individual filters
# Use the linear fit to state what the M was at day 50 [THIS CODE SHOULD BE HERE SOMEWHERE]
## The light curves are interpolated to 50 days after explosion [Gall]
1. Linearly interpolate restframe V-I colors using g'r'i'z'

# These colors were used to "warp" a day 50 template spectrum, redshifted appropriately, following Nugent02, except they use spline fits to underlying colors to adjust template rather than a reddening law
# This is sncosmo [spectrum]
2. In sncosmo, warp a Day50 template spectrum with these colors using a reddening law

# This spectrum was then de-redshifted and the restframe V-I color of the SN was calculated 
3. De-redshift the spetrum and calculate restframe V-I color

- The V and I apparent magnitudes 
- M_I (absolute) 
- From M_I, subtract off I from sncosmo
- mu = I - M_I 

Plot a histogram of M_Is for all t0s 
"""

"""
Calculate rest-frame V-I colors during the plateau by:
Estimate explosion date
Linear interpolate the g'r'i'z colors at day 50(1+z). 
To each of the individual filters use the linear fit to state what the M was at day 50 
"""

def import_color_lc(file, lc_dir):
    colnames_color_lc = ['flags', 'flux', 'fluxerr',  'id', 'lim_mag', 'mjd', 'zp']
    file_loc = os.path.join(lc_dir, file)
    name = os.path.splitext(file)[0]
    df = pd.read_csv(file_loc, comment='#', delim_whitespace=True, names = colnames_color_lc)
    return df, name

df_g_band, g_filt = import_color_lc(f'{ZTF_name}_g_lightcurve.csv', filters_dir)
df_r_band, r_filt = import_color_lc(f'{ZTF_name}_r_lightcurve.csv', filters_dir)
df_i_band, i_filt = import_color_lc(f'{ZTF_name}_i_lightcurve.csv', filters_dir)

bands = [df_g_band, df_r_band, df_i_band]
filter_names = [g_filt, r_filt, i_filt]
d_bands = dict(zip(filter_names, bands))


# In[ ]:


def f26(f, fsig, zp):
    f26 = f*10**((26-zp)/2.5)
    sig26 = fsig*10**((26-zp)/2.5)
    return f26

def interpolate_lc(band, band_f_conv):
    interpolate_band_obj = scipy.interpolate.interp1d(band.mjd, band_f_conv)
    spec_axis_new = np.arange(min(band.mjd), max(band.mjd), 10)
    band_flux = interpolate_band_obj(spec_axis_new)
    return interpolate_band_obj, spec_axis_new, band_flux

d_flux_calcs = {}
for filt_name, band in d_bands.items():
    band['f_conv'] = band.apply(lambda x: f26(x['flux'], x['fluxerr'], x['zp']), axis = 1)
    band_y_interp, band_specaxis, band_flux = interpolate_lc(band, band.f_conv)
    pair = band_y_interp, band_flux, band_specaxis, band.mjd.values, band.f_conv.values
    if filt_name not in d_flux_calcs:
        d_flux_calcs[filt_name] = []
    d_flux_calcs[filt_name].append(pair)

d_flux_calcs_lst = [[k, *v] for k, lst in d_flux_calcs.items() for v in lst]
df_fluxes = pd.DataFrame(d_flux_calcs_lst, columns = ['filt', 'y_interp_obj', 'interp_fluxes', 'new_specaxis', 'mjds', 'f_conv'])


# In[ ]:


day50 = (explosion_dates[0]+50)#*(1+z_input)
'''
ask about 1+z term-- out of bounds 
'''
plt.figure(figsize = (10, 6))

for i in range(len(df_fluxes.filt.values)):
    plt.plot(df_fluxes.mjds[i], df_fluxes.f_conv[i], 'o', alpha = 0.4)
    plt.plot(df_fluxes.new_specaxis[i], df_fluxes.interp_fluxes[i], '-', alpha = 0.4, label = f'{df_fluxes.filt[i][len(ZTF_name)+1:-11]}')
    plt.plot(day50, df_fluxes.y_interp_obj[i](day50), 'o', markersize = 10, label = f'interp-{df_fluxes.filt[i][len(ZTF_name)+1:-11]}')

plt.axvline(day50)
plt.xlim(58710, 58810)
plt.ylim(300,4000)

plt.grid()
plt.minorticks_on()
plt.legend(loc = 0)
plt.legend()


# In[ ]:


zp = 26

i_flux = iband_y_interp(day50)
g_flux = gband_y_interp(day50)
r_flux = rband_y_interp(day50)

def flux_to_mag(flux):
    mag = -2.5*np.log10(flux) + zp
    return(mag)

app_mag_i = flux_to_mag(i_flux)
app_mag_g = flux_to_mag(g_flux)
app_mag_r = flux_to_mag(r_flux)

print(f'apparent mag_i = {app_mag_i}\napparent mag_g = {app_mag_g}\napparent mag_r = {app_mag_r}')


# In[ ]:


"""
STEP 2.
These colors were used to "warp" a day 50 template spectrum, redshifted appropriately, following Nugent02, except they use spline fits to underlying colors to adjust template rather than a reddening law
This is sncosmo [spectrum]

In sncosmo, warp a Day50 template spectrum with these colors using a reddening law
"""

def import_day50(file, day50_file_dir):
    colnames_day50_spectrum = ['phase', 'wave', 'flux']
    file_loc = os.path.join(day50_file_dir, file)
    df = pd.read_csv(file_loc, comment='#', delim_whitespace=True, names = colnames_day50_spectrum)
    return df

df_day50_spec = import_day50(f'{ZTF_name}_day50_spectrum.dat', day50_file_dir)


# In[ ]:


wave_sncosmo = df_day50_spec.wave.values
flux_sncosmo = df_day50_spec.flux.values
spectrum = sncosmo.Spectrum(wave_sncosmo, flux_sncosmo)
model = sncosmo.Model(source= 'nugent-sn2p')


# In[ ]:


# Find start and end of the input spectrum and call them the first and last spline fit
xsp_1 = spectrum.wave[0]
xsp_end = spectrum.wave[-1]

def flux_to_mag(flux):
    mag = -2.5*np.log10(flux)
    return(mag)

fluxes = spectrum.flux
mags = flux_to_mag(spectrum.flux)
waves = spectrum.wave


# In[ ]:


# get filters

filt_g = sncosmo.get_bandpass('ztfg')
filt_r = sncosmo.get_bandpass('ztfr')
filt_i = sncosmo.get_bandpass('ztfi')

n = 3 #number of filters
nrmflt = filt_r # filter for normalizing
nspline = (2*n)+1

filters = [filt_g, filt_r, filt_i]
wave_g = filt_g.wave
trans_g = filt_g.trans
wave_r = filt_r.wave
trans_r = filt_r.trans
wave_i = filt_i.wave
trans_i = filt_i.trans

#find midpoints:
def find_transmax(filt):
    return np.max(filt.trans)

# midpoints of all filts
g_maxtrans = find_transmax(filt_g)# peak transmission of filter g
r_maxtrans = find_transmax(filt_r)
i_maxtrans = find_transmax(filt_i)

#index of midpts
g_maxtrans_idx = np.where(trans_g == g_maxtrans) # peak transmission index of filter g
r_maxtrans_idx = np.where(trans_r == r_maxtrans)
i_maxtrans_idx = np.where(trans_i == i_maxtrans)

#wave at midpts 
g_wavepeak = wave_g[g_maxtrans_idx]
r_wavepeak = wave_r[r_maxtrans_idx]
i_wavepeak = wave_i[i_maxtrans_idx]

#first and last waves
g_firstwave = wave_g[0]
i_lastwave = wave_i[-1]

g_firstwave_idx = np.where(wave_g == g_firstwave)
i_lastwave_idx = np.where(wave_i == i_lastwave)

# Test to copy fortran code lines 83-92
xsp = np.zeros(2*n+1); xsp[0] = xsp_1 ; xsp[-1] = xsp_end ; xsp[3] = r_wavepeak
xvflt = [g_wavepeak, r_wavepeak, i_wavepeak]

for i in [1, 2*n-1, 1]:
    xsp[i] = xvflt[int(i/2)]
    
for i in [2, 2*n-2, 1]:
    xsp[i] = (xsp[i-1] + xsp[i+1])/2


# In[ ]:


#idx of target spectrum at wavelength xsp(2) == first midpoint (g)
idx1 = (np.abs(waves - g_wavepeak)).argmin()

#idx of target spectrum at wavelength xsp(2n) == 2nd midpoint (r)
idx2 = (np.abs(waves - r_wavepeak)).argmin()

#divide idx1 and idx2 by 2, assuming everything spaced evenly
k1 = int(round(idx1/2))
k2 = idx2 + int(round((len(waves)-idx2)/2))

# sum up flux between wave_1 and k1, and then k1 and sp1, and take the ratio

fy1 = sum(fluxes[0:k1])/sum(fluxes[k1:idx1])

# sum up flux between k2 and sp1, and then sp2 and k2, and take the ratio

fy2 = sum(fluxes[k2:])/sum(fluxes[idx2:k2])


# In[ ]:


mags_ztfg = model.bandmag('ztfg', 'ab', [0])[0]
mags_ztfr = model.bandmag('ztfr', 'ab', [0])[0]
mags_ztfi = model.bandmag('ztfi', 'ab', [0])[0]

band_mags = [mags_ztfg, mags_ztfr, mags_ztfi]
mags_obs = [app_mag_g, app_mag_r, app_mag_i] #g, r, i
alam = [0.099, 0.068, 0.051] #Milky Way extinctions
mags_input = [mags_obs - alam for mags_obs, alam in zip(mags_obs, alam)]
mgdiff_r = mags_input[1] - mags_ztfr
fct = 10**(-0.4*mgdiff_r)

warpflux = fct * fluxes
warpwave = waves

spectrum_warp = sncosmo.Spectrum(warpwave, warpflux, time=50)

filt_names = ['ztfg', 'ztfr', 'ztfi']


# In[ ]:


i = 0
n = 3
ysp = np.zeros(2*n+1)

ym = 1*np.exp(31)
y2n = np.zeros(2*n+1)
met = False

# create a new spectrum 
while i < 20:
    band_mags = spectrum_warp.bandmag(filt_names, 'ab') #Magnitude through the given bandpass(es), and for the given magnitude system(s).
    mgdf = mags_input - band_mags
    # print(i, mgdf)
    # print(band_mags)
    
    for j in range(len(filt_names)):
        if np.abs(mgdf[j]) < 0.005:
            met = True
    if met == True: break
        
    mgdf = mags_input - band_mags
    i+=1
    for k in range(1, 2*n, 2): #1 and 3 and 5
        ysp[k] = (10**(-0.4 * mgdf[int(k/2)]))
    for k in range(2, 2*n-1, 2): # 2 and 4
        ysp[k] = (ysp[k-1] + ysp[k+1])/2 
    ysp[0] = fy1 * ysp[1] #0
    ysp[-1] = fy2 * ysp[2*n-1] #6
    tck = splrep(x = xsp, y = ysp)
    splint = splev(warpwave, tck)
    
    warpflux = splint * warpflux
    warpwave = waves
    spectrum_warp = sncosmo.Spectrum(warpwave, warpflux, time=50)
    
    fy1 = 1.0
    fy2 = 1.0
    
#     print(xsp)


# In[ ]:


xs = warpwave
# plt.plot(xs, splint, 'g', lw=3)
plt.figure(figsize = (10, 6))
plt.plot(xs, warpflux/np.max(warpflux), label ='spline fit')
plt.plot(spectrum.wave, spectrum.flux/np.max(spectrum.flux), label = 'original')
plt.title('spline fit to spectrum \n vs. \n original spectrum')
plt.legend()
plt.xlabel('observed wavelength [$\AA$]')
plt.ylabel('flux')


# In[ ]:


def calc_kcorr(z, spectrum, flt1, flt2):
    
    spdered = sncosmo.Spectrum(spectrum.wave/(1+z), spectrum.flux)
    mag1 = spdered.bandmag(flt1, 'ab') #restframe bandpass
    mag2 = spectrum.bandmag(flt2, 'ab') #measured bandpass
    
    kcorr = mag2 - mag1 + 2.5*np.log10(1 + z)
    return kcorr


# In[ ]:


kcorr = calc_kcorr(z_input, spectrum_warp, 'ztfg', 'besselli')

#de-redshift spectrum to get restframe V and I absolute mags

spdered = sncosmo.Spectrum(spectrum.wave, spectrum_warp.flux)

# calculate m_V (magnitude in V) and m_I (magnitude in I)

new_filtnames = ['bessellv', 'besselli']
V_rest, I_rest = spdered.bandmag(new_filtnames, 'vega')
print(f'V = {V_rest}\nI = {I_rest}')