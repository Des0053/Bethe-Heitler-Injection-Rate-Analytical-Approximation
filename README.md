# Bethe-Heitler-Injection-Rate-Analytical-Approximation
Python class-script providing an analytical approximation of the pair distribution produced due to Bethe-Heitler pair production.

Function Q_inj describes the produced electron/positron distribution as a function of the electron/positron Lorentz factor, $\gamma_e$, for the interaction of a single proton of Lorentz factor $\gamma_p$
interacting with a single photon of energy $\epsilon$ (in $m_e c^2$ units). 

Class input parameters:
 * g: pair Lorentz factor grid in logarithmic scale (np.array form)

 * x: photon energy in $m_e c^2$ units and in logarithmic scale 

 * g_p: proton Lorentz factor value in logarithmic scale

 * R: radius of the source in lenear scale (cm)

 * ind: unit selection:
   
    - 0: $s^{-1}$

    - 1: $cm^{-3} s^{-1}$

    - 2: $\frac{1}{\tilde V \tilde t}$ with $\tilde V= \frac{V (cgs)}{ \sigma_T R}$ and $\tilde t= t \frac{c}{R}$

   Class output:
    * BH_inj_rate: q_BH variable returns an array of the BH injection rate computed on each point on the $\gamma_e$ (g in script symbols) given by the user **for single photon-single proton** interaction.
   
    * BH_inj_rate_arr: q_BH variable returns an array of size [len(x), len(g)] that refers to the BH injection rate computed for a single photon of $\epsilon \in x$ energies with a single photon. The photon energies have to be given in an array and in a logarithmic scale. Each row of the output array describes the produced pair distribution across the $\gamma_e$ grid, by the interaction of the photon with a single photon of a specific energy.

The class also calculates the normalization, A, and the slope, s, of the injection rate function, Q_inj, given as:

$Q_{inj}(\gamma_e)= \mathbf{A(\gamma_p,  \epsilon)} \cdot \exp \left[ -\frac{\left[\log_{10}\left(\frac{\gamma_e}{\gamma_{e,\rm pk}}\right)\right]^{ \mathbf{s(\gamma_p \epsilon)}}}{2 {a_1}^2}- a_2^2 \left(\frac{\gamma_{e, \rm pk}}{\gamma_e}-1\right)^2 -a_3 \frac{\gamma_e}{\gamma_{e, \rm cr}}  \right] $


The modelling, properties and benchmarking of the function are included in the paper "A closer look at the electromagnetic signatures of Bethe-Heitler pair production process in blazars" by D. Karavola and M. Petropoulou, 2024 (to be submitted to JCAP).

