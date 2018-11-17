# LWaves

Scripts:

1. car2sph: Converts Cartesian (txt) to Spherical Coordinates of 3D map 

<p float="left">
<img src="../img's/LW008_2D.png" width="265"/>
<img src="../img's/LW008_3D.png" width="265"/>
<img src="../img's/LW008.png" width="265"/>
</p>

2. LiDARRays: Transform spherical map (RTPA) to sliced map with grouped points for each LiDAR ray

<p float="left">
<img src="../img's/LW008_reduced.png" width="265"/>
<img src="../img's/LW008_rays2_0.png" width="265"/>
</p>

3.1 WaveTimeStack: Transform 3D LiDAR map (RTPAXYZ_rays) to Wave Time-stack of water surface elevation 

<p float="left">
<img src="../img's/LW008_RZt.png" width="265"/>
<img src="../img's/LW008_RZ.png" width="265"/>
<img src="../img's/LW008_Zt.png" width="265"/>
</p>

4. LiDARWaveSpectrum: Spectral wave analysis for LiDAR time-series (RZt)

<p float="left">
<img src="../img's/LW008_Zt_V2.png" width="265"/>
<img src="../img's/LW008_Pn.png" width="265"/>
<img src="../img's/LW008_Ef_1.png" width="265"/>
</p>

5. OSSI4 WaveSpectrum: Spectral wave analysis for OSSI (pressure sensor) time-series (RZt)

<p float="left">
<img src="../img's/OSSI4_Pt.png" width="265"/>
<img src="../img's/OSSI4_Pn.png" width="265"/>
<img src="../img's/OSSI4_Ef_1.png" width="265"/>
</p>

6. Validation: Validation of spectral wave characteristics of LiDAR w.r.t. OSSI's (Spectrum)

<p float="left">
<img src="../img's/Ef_2.png" width="400"/>
<img src="../img's/HT.png" width=400"/>
</p>

7. Precision: Sensitivity analysis of precision on spectral wave characteristics of LiDAR (Precision)

<p float="left">
<img src="../img's/Ef_4.png" width="265"/>
</p>

3.2 RunupTimeStack: Transform 3D LiDAR map (RTPAXYZ_rays) to Run-up Time-stack based on reflectivity threshold

<p float="left">
<img src="../img's/LR022_RZtA.png" width="400"/>
</p>

Functions:

1. zero_crossing: Zero crossing analysis of wave data

2. rms_wave_height: Calculation of H_rms

3. significant_wave_height: Compute the signficant wave height 

4. wave_spectrum: Compute variance spectral density spectrum of the time-series and its 90% confidence intervals

5. spectral_moment: Calculate the n^th-order spectral moment for a given frequency band [fmin fmax]

6. p2sse: Corrects a pressure signal by means of linear wave theory back to the sea surface elevation (see)

7. wavenumberGuo: Compute wavenumber according to Guo 2002
