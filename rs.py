# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
num_frames = 100
timestep = 0.05
# timestep = 1
r_vals = np.linspace(0.0, 10.0, 501)
r_vals = r_vals[:-1]
# r_vals[0] = 1e-10
# M_vals = np.linspace(0, 25.0, len(r_vals)+1)
# M_vals = M_vals[:-1]
# M_vals[0] = 1e-10
t_vals = np.linspace(0.0, (num_frames)*timestep, num_frames+1)
t_vals = t_vals[1:]
# t_vals[0] = 1e-10
print(f"shape of r_vals : {r_vals.shape}")
# print(f"shape of M_vals : {M_vals.shape}")
print(f"shape of t_vals : {t_vals.shape}")
# print(r_vals)
# print(M_vals)
print(t_vals)

# %%
# s_vals = np.linspace(0.0, 2*np.pi, 101)
s_vals = np.linspace(0.2, 20.0, 101)
s_vals = s_vals[:-1]

# %%
num_alphas = 10
alphlowlim = 0.0
alphhighlim = 0.9
num_ds = 50
dlowlim = 0.0
dhighlim = 5.0

# %%
alpha_vals = np.linspace(alphlowlim, alphhighlim, num_alphas)
# alpha_vals[0] = 1e-10 # set alpha=0 to a small value
# alpha_vals = np.ones(num_alphas)*0.5
print(alpha_vals)

d_vals = np.linspace(dlowlim, dhighlim, num_ds+1)
d_vals = d_vals[:-1]
# d_vals[0] = 1e-10 # set d=0 to a small value
# d_vals = np.ones(num_ds)*2.0
print(d_vals)

# print the shape of the data
print(f"shape of alpha_vals : {alpha_vals.shape}")
print(f"shape of d_vals : {d_vals.shape}")

# %%
# load Z
Z = np.load("Z.npy")
print(f"shape of Z : {Z.shape}")
# print(Z)

# %%
# load u_pred
u = np.loadtxt("u.txt")
print(f"shape of u : {u.shape}")
# print(u_pred)

# %%
# # plot the u_pred(s, r) for some s values
# plt.figure(figsize=(10, 8))
# plt.title("u(s, r) for some s values")
# plt.xlabel("r")
# plt.ylabel("u(s, r)")
# plt.grid()
# # plt.ylim(0)
# for i in range(0, s_vals.shape[0], s_vals.shape[0]//10):
#     if i == 0:
#         continue
#     plt.plot(r_vals, u[i, :], label=f"s = {s_vals[i]:.2f}")

# plt.legend()
# plt.show()

# %% [markdown]
# $P_s^{i+1}(\alpha, D) = P_s^i(\alpha, D) \int \frac{u(s,r)}{u_{pred}(s,r)}Z(s, r, \alpha, D)dr$
# 
# where 
# 
# $u_{pred}(s,r) = \int \int P_s^{i}(\alpha, D) Z(s, r, \alpha, D) d\alpha dD$
# 

# %%
u_pred = np.zeros(u.shape)
print(f"shape of u_pred : {u_pred.shape}")

# %%
# initialize P_alpha_d
# P_alpha_d = np.ones((u.shape[0], num_alphas, num_ds))

# P_alpha_d = np.zeros((u.shape[0], num_alphas, num_ds))


# P_alpha_d[:, 2, 20] = 1

# initialize P_alpha_d to be a gaussian distribution with peak at (alpha, d) = (0.5, 2.5)
alpha0 = 0.5
d0 = 3.01
P_alpha_d = np.zeros((u.shape[0], num_alphas, num_ds))
for i in range(num_alphas):
    for j in range(num_ds):
        P_alpha_d[:, i, j] = np.exp(-((alpha_vals[i] - alpha0)**2 + (d_vals[j] - d0)**2) / 0.1)

# print(P_alpha_d)
# normalize P_alpha_d
P_alpha_d /= np.trapz(np.trapz(P_alpha_d, x=d_vals, axis=2), x=alpha_vals, axis=1)[:, None, None]

###### fixed alpha
# norm = (np.trapz(P_alpha_d, x=d_vals, axis=2))[:, :, None]
# P_alpha_d = P_alpha_d / norm
###### fixed alpha

# P_alpha_d = np.ones((num_alphas, num_ds))
print(f"shape of P_alpha_d : {P_alpha_d.shape}")

P_alpha_d_old = P_alpha_d.copy()

# %%
def contour_plt(P_alpha_d, i=0, save=False, name=""):
    """
    Visualise P_alpha_d as a countour plot using contourf
    """
    # plot
    # rescale the y axis to be in the range of 0 to 2
    # rescale the x axis to be in the range of 0 to 1000
    plt.figure(figsize=(4, 4))
    x = np.linspace(dlowlim, dhighlim, num_ds)
    y = np.linspace(alphlowlim, alphhighlim, num_alphas)
    plt.contourf(x, y, P_alpha_d[i], 20, cmap='RdGy')
    plt.xlabel('d')
    plt.ylabel('alpha')
    plt.colorbar()
    if save:
        plt.savefig(name)
        plt.close()
    else:
        plt.show()

# %%
contour_plt(P_alpha_d, (s_vals.shape[0]-1)//2, save=True, name=f"./P_alpha_d_img/init_contour.png")

# %% [markdown]
# $P_s^{i+1}(\alpha, D) = P_s^i(\alpha, D) \int \frac{u(s,r)}{u_{pred}(s,r)}Z(s, r, \alpha, D)dr$
# 
# where 
# 
# $u_{pred}(s,r) = \int \int P_s^{i}(\alpha, D) Z(s, r, \alpha, D) d\alpha dD$
# 

# %%
def integ1(u_pred, P_alpha_d, alpha_vals, d_vals, r_vals, s_vals):
    # implement $u_{pred}(s,k) = \int \int P_s^{i}(\alpha, D) Z(\alpha, D, s, k) d\alpha dD$
    
    ############# correct but slow code with Z ind diff #############
    # for i in range(u_pred.shape[0]):
    #     for j in range(u_pred.shape[1]):
    #         u_pred[i, j] = np.trapz(np.trapz(P_alpha_d[i] * Z[:, :, i, j], d_vals, axis=1), alpha_vals, axis=0)
    ############# correct but slow code #############
    
    ############# vectorized code #############
    u_pred = np.trapz(np.trapz(P_alpha_d[:, None, :, :] * Z, d_vals, axis=3), alpha_vals, axis=2)
    ############# vectorized code #############
    
    ###### fixed alpha
    # u_pred = np.mean(np.trapz(P_alpha_d[:, None, :, :] * Z, d_vals, axis=3), axis=2)
    ###### fixed alpha

    # normalize u_pred for each s
    u_pred = u_pred / np.trapz(u_pred, r_vals, axis=1).reshape(-1, 1)

    return u_pred

def integ2(u, u_pred, P_alpha_d, alpha_vals, d_vals, r_vals, s_vals):
    # implement $P_s^{i+1}(\alpha, D) = P_s^i(\alpha, D) \int \frac{u(s,k)}{u_{pred}(s,k)}Z(\alpha, D, s, k)dk$

    # calculate ratio = u/u_pred
    u_pred[u_pred < 1e-7] = 0
    ratio = np.divide(u, u_pred, out=np.ones_like(u), where=u_pred != 0)
    # print(f"mean ratio : {np.mean(ratio)}")
    # ratio = abs(ratio)

    ############# correct but slow code with Z ind diff #############
    # for i in range( P_alpha_d.shape[0]):
    #     for j in range(P_alpha_d.shape[1]):
    #         for k in range(P_alpha_d.shape[2]):
    #             P_alpha_d[i, j, k] *= np.trapz(ratio[i] * Z[j, k, i, :], r_vals)
    ############# correct but slow code #############
    
    ############# vectorized code #############
    P_alpha_d *= np.trapz(ratio[:, :, None, None] * Z, r_vals, axis=1)
    ############# vectorized code #############
                
    # normalize P_alpha_d for each s
    norm = np.trapz(np.trapz(P_alpha_d, x=d_vals, axis=2), x=alpha_vals, axis=1)[:, None, None]
    P_alpha_d = P_alpha_d / norm

    ###### fixed alpha
    # norm = (np.trapz(P_alpha_d, x=d_vals, axis=2))[:, :, None]
    # P_alpha_d = P_alpha_d / norm
    ###### fixed alpha

    return P_alpha_d, np.mean(ratio)

# %%
# iterate until delta is small enough
delta = 10000
num_iter = 0
ratio = 1
# while delta > 1e-3 and num_iter < 1:
while delta > 1e-7:
    u_pred = integ1(u_pred, P_alpha_d, alpha_vals, d_vals, r_vals, s_vals)
    P_alpha_d, ratio = integ2(u, u_pred, P_alpha_d, alpha_vals, d_vals, r_vals, s_vals)
    delta = np.sum((P_alpha_d - P_alpha_d_old)**2)
    # print(f"delta : {delta}")
    P_alpha_d_old = P_alpha_d.copy()
    num_iter += 1
    if num_iter % 100 == 0:
        contour_plt(P_alpha_d, (s_vals.shape[0]-1)//2, save=True, name=f"./P_alpha_d_img/contour_{num_iter}.png")
        print(f"num_iter : {num_iter}")
        print(f"mean ratio : {ratio}")
        print(f"delta : {delta}")
print(f"num_iter : {num_iter}")

# %%
contour_plt(P_alpha_d, (s_vals.shape[0]-1)//2, save=True, name=f"./P_alpha_d_img/final_contour_middle.png")
print(f"ratio = {np.mean(u/u_pred)}")

# %%
for i in range(0, len(P_alpha_d), s_vals.shape[0]//10):
    contour_plt(P_alpha_d, i, save=True, name=f"./P_alpha_d_img/final_contour_{i}.png")

# %%
# save P_alpha_d
np.save("P_alpha_d.npy", P_alpha_d)


