G-tools
=======

Short description
-----------------

The *interactive* scripts \`aw2giw\` and \`aw2gtau\` calculate the Green's function (either on the imaginary time axis or on the Matsubara axis) from a given spectral function. The latter is read from a two-column (first column: frequency, second column: spectral function) file.

`hk2gwfw` is an *interactive* script which processes a Hamiltonian in the `*.hk` format. Output is one of the following:

-   Green's function
-   Spectral function
-   Hybridization function

Calculations can be done either on the real frequency axis, or on the Matsubara axis. What exactly to print is decided by the user. For instance, you can print only diagonal or only nondiagonal elements. For complex quantities, you can print the real part, the imaginary part, both (two columns) or the absolute value.

How to install and run
----------------------

You need:

-   a C++ compiler supporting the C++11 standard
-   Eigen headers (for `hk2gwawfw`)
-   an MKL library (recommended for `hk2gwawfw`)

Adjust the first three lines of `Makefile` accordingly and run

    make all

This will create three executables: `aw2giw.o`, `aw2gtau.o`, and `hk2gwawfw.o`. All three codes are interactive, simply execute them. Sample spectral functions and H(k) Hamiltonians can be found in the folder `examples`.

Physics behind
--------------

*(needs MathJax, or check in `readme.pdf`)* From a spectral function to the corresponding non-interacting Green's function:

\[
G(\tau) = \int_{\omega_{min}}^{\omega_{max}}{d\omega \frac{e^{(\mu-\omega)\tau}}{1 + e^{(\mu-\omega)\beta}}A(\omega)}
\]

\[
G(i\omega_n) = \int_{\omega_{min}}^{\omega_{max}}{d\omega \frac{1}{i\omega_n - \omega + \mu}A(\omega)}
\]

From a matrix-valued non-interacting Hamiltonian to the corresponding matrix-valued Green's function (real and Matsubara axes):

\[
{\bf G}(\omega) = \frac{1}{N_k}\sum_{k=1}^{N_k}\left[\left(\omega+i\delta+\mu\right){\bf 1}-{\bf H}(k)\right]^{-1}
\]

\[
{\bf G}(i\omega_n) = \frac{1}{N_k}\sum_{k=1}^{N_k}\left[\left(i\omega_n+\mu\right){\bf 1}-{\bf H}(k)\right]^{-1}
\]

From a matrix-valued non-interacting Hamiltonian to the corresponding matrix-valued hybridization function (real and Matsubara axes):

\[
{\bf \Delta}(\omega) = \left(\omega+\mu\right){\bf 1}-{\bf G}(\omega)^{-1}
\]

\[
{\bf \Delta}(i\omega_n) = \left(i\omega_n+\mu\right){\bf 1}-{\bf G}(i\omega_n)^{-1}
\]
