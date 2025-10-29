import numpy as np
import subprocess
import matplotlib.pyplot as plt

def numerical_value(kappa, mu, theta, phi):
    cmd = ['./spinor_spherical_harmonic_valid', str(kappa), str(mu), str(theta), str(phi)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.strip().split('\t')
    values = [float(line) for line in lines]
    return np.array([values[0] + 1j * values[1] ,values[2] + 1j * values[3]])

def verify_sum(kappa, theta, phi):
    mu_list = np.arange(-abs(kappa)+0.5, abs(kappa), 1.0)
    total = np.zeros((2,2), dtype=complex)
    for mu in mu_list:
        val = numerical_value(kappa, mu, theta, phi)
        total += np.outer(val, np.conj(val))
    return total

def theoretical_value(kappa):
    j = abs(kappa) - 0.5
    theo = (2*j + 1)/(8*np.pi)
    return np.array([[theo, 0], [0, theo]], dtype=complex)


fig, axs = plt.subplots(1, 4, figsize=(16,4))
kappa_values = [-2, 1, 2, 3]
theta_val = np.linspace(0, np.pi, 50)
phi_val   = np.linspace(0, 2*np.pi, 50)
Theta, Phi = np.meshgrid(theta_val, phi_val)
for ax, kappa in zip(axs, kappa_values):
    relerror = np.zeros(Theta.shape)
    ther_val = theoretical_value(kappa)[0,0]
    for i in range(Theta.shape[0]):
        for j in range(Theta.shape[1]):
            sum_val = verify_sum(kappa, Theta[i,j], Phi[i,j])
            relerror[i,j] = np.abs(sum_val[0,0] - ther_val)/ abs(ther_val)
    c = ax.pcolormesh(Theta, Phi, relerror, shading='auto', cmap='viridis')
    ax.set_title(f'kappa={kappa}')
    ax.set_xlabel('Theta')
    ax.set_ylabel('Phi')
fig.colorbar(c, ax=axs[-1], orientation='vertical', label='Relative Error')
plt.tight_layout()
plt.savefig('validation_error.png', dpi=300)
plt.close()

phi=np.pi/4
fig, ax=plt.subplots(figsize=(8,8))
theta=np.linspace(0, np.pi, 100)
kappa=3
relerror=np.zeros(theta.shape)
for i in range(theta.shape[0]):
    sum_val = verify_sum(kappa, theta[i], phi)
    ther_val = theoretical_value(kappa)[0,0]
    relerror[i] = np.abs(sum_val[0,0] - ther_val)/ abs(ther_val)
ax.scatter(theta, relerror, label=f'kappa={kappa}, phi={phi} rad', marker='+')
ax.set_xlabel(r'$\theta$ [rad]')
ax.set_ylabel('Relative Error')
ax.set_ylim(0, 1e-12)
ax.set_title(f'kappa={kappa}, phi={phi} rad')
ax.legend()
plt.tight_layout()
plt.savefig('validation_error_line.png', dpi=300)