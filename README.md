# Bethe-Heitler-Injection-Rate-Analytical-Approximation
Python class-script providing an analytical approximation of the pair distribution produced due to Bethe-Heitler pair production.

Function Q_inj describes the produced electron/positron distribution as a function of the electron/positron Lorentz factor, $\gamma_e$, for the interaction of a single proton of Lorentz factor $\gamma_p$
interacting with a single photon of energy $\epsilon$ (in $m_e c^2$ units). 

Class input parameters:
$\begin{itmerize}

\item g: pair Lorentz factor grid in logarithmic scale (np.array form)

\item x: photon energy in $m_e c^2$ units and in logarithmic scale 

\item g_p: proton Lorentz factor value in logarithmic scale

\item R: radius of the source in lenear scale (cm)

\end{itemize}
$

The modelling and properties of the function are desrcibed in the paper "A closer look at the electromagnetic signatures of Bethe-Heitler pair production process in blazars" by D. Karavola and M. Petropoulou (to be submitted).

