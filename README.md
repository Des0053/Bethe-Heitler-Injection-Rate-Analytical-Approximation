# Bethe-Heitler-Injection-Rate-Analytical-Approximation
Python class-script providing an analytical approximation of the pair distribution produced due to Bethe-Heitler pair production.

Function $Q_{inj}$ describes the produced electron/positron distribution as a function of the electron/positron Lorentz factor, $\gamma_e$, for the interaction of a single proton of Lorentz factor $\gamma_p$
interacting with a single photon of energy $\epsilon$ (in $m_e c^2$ units). 

Class input parameters:
 * g: pair Lorentz factor grid in logarithmic scale (np.array form)

 * x: photon energy in $m_e c^2$ units and in logarithmic scale 

 * g_p: proton Lorentz factor value in logarithmic scale

 * R: radius of the source in lenear scale (cm)

 * ind: unit selection:
   
    - 0: $s^{-1}$

    - 1: $cm^{-3} s^{-1}$

    - 2: $(\tilde V \tilde t)^{-1}$ with $\tilde V= V [cm^{-3}] (\sigma_T R)^{-1}$ and $\tilde t= t[s] \frac{c}{R}$

   Class output:
    * BH_inj_rate: $q_{BH}$ variable returns an array of the BH injection rate computed on each point on the $\gamma_e$ (g in script symbols) given by the user **for single photon-single proton** interaction.
   
    * BH_inj_rate_arr: $q_{BH}$ variable returns an array of size [len(x), len(g)] that refers to the BH injection rate computed for a single photon of $\epsilon \in x$ energies with a single photon. The photon energies have to be given in an array and in a logarithmic scale. Each row of the output array describes the produced pair distribution across the $\gamma_e$ grid, by the interaction of the photon with a single photon of a specific energy.

The class also calculates the normalization, $A$, and the slope, $s$, of the injection rate function, $Q_{inj}$, given as:

$Q_{inj}(\gamma_e)= A(\gamma_p,  \epsilon) \cdot \exp \left[ -\frac{\left[\log_{10}\left(\frac{\gamma_e}{\gamma_{e,\rm pk}}\right)\right]^{ \mathbf{s(\gamma_p \epsilon)}}}{2 {a_1}^2}- a_2^2 \left(\frac{\gamma_{e, \rm pk}}{\gamma_e}-1\right)^2 -a_3 \frac{\gamma_e}{\gamma_{e, \rm cr}}  \right] $

**Comment:** The above function was modelled using the Leptohadronic, radiative, Monte-Carlo code ATHEvA (Mastichiadis et. al (2005), Dimitrakoudis et. al (2012)). The aforementioned code takes into account interactions between protons and photons with interaction energy, $\gamma_p E$, values $\leq 10^4$. As a result the function, $q_{BH}$, described above, has the same limitation and has been benchmarked for $\gamma_p E \leq 10^4$. This limitation is embeded by the AM_flag. When this flag is equal to 1 (default value), $\gamma_p E >10^4$ interactions are not taken into account. The user can change the AM_flag to 0 to include interaction of $\gamma_p E>10^4$ but with caution since the validity of the function in that regime is not known.   


The modelling, properties and benchmarking of the function are included in the paper "A closer look at the electromagnetic signatures of Bethe-Heitler pair production process in blazars" by D. Karavola and M. Petropoulou, 2024 (submitted to JCAP, arxiv:2401.05534).

