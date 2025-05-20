import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"]   = "28"

def MN_theo(T):
    T_c = 2/(np.log(1+np.sqrt(2)))
    return np.heaviside(T_c-T , 0)*np.abs((1-1/(np.sinh(2/T))**4))**(1/8)

def setup(title, xlabel, ylabel):
    plt.figure(figsize=(16, 9))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks()
    plt.yticks()
    plt.locator_params(nbins=6)
    plt.grid()

def dataset(N):
    def load_data(a, b, c, boundary):
        text = "2D_Data/2D_x_qN_zNS_K.npy"
        text = text.replace("x", str(a))
        text = text.replace("q", str(b))
        text = text.replace("z", str(c))
        text = text.replace("K", str(boundary))
        return np.load(text)
    

    def create(name):
        arr = []
        Ns_arr  = [1000, 10000, 1000, 10000]
        bounds  = ["per", "per", "free", "free"]
        for i in range(4):
            arr.append(load_data(name, N, Ns_arr[i], bounds[i]))
        return np.array(arr).reshape(-1, 20)
    
    U = create("U")
    C = create("C")
    M = create("M")
    return [U, C, M]
    



def plot_data(data, title, x_axis, y_axis, magn, heat, Cmax, name):
    text = "2D_Plots/2D_Z_X.pdf"
    text = text.replace("X", str(N))
    text = text.replace("Z", name)
    
    T   = np.linspace(0.2, 4, 20)
    T_c = 2/(np.log(1+np.sqrt(2)))
    
    Ns_arr = ["10^3", "10^4", "10^3", "10^4"]
    bounds = [r"periodic:", r"periodic:", r"free: \qquad         ", r"free:    \qquad       "]

    
    setup(title, x_axis, y_axis)
    if magn:
        T   = np.linspace(0.2, 4, 1000)
        plt.plot(T,  MN_theo(T), color="black", ls="--", linewidth=5, label="Theory")
        T   = np.linspace(0.2, 4, 20)
    if heat:
        plt.vlines(T_c, 0, Cmax, color="black", ls="--", linewidth=5, label="$T_c$")
    for i in range(4):
        plt.plot(T, data[i], label=r"X".replace("X", bounds[i]) + r" $N_S=Y$".replace("Y", Ns_arr[i]),linewidth=4,marker="o", markerfacecolor='white')
    plt.legend()
    plt.savefig(text)
    plt.close()


N_arr       = [10, 50, 100]
title_arr   = ["Internal Energy per Spin, N=X", "Specific Heat per Spin, N=X", "Absolute Magnetization per Spin, N=X"]
y_arr       = ["U/N$^2$", r"C/N$^2$", r"$|M|$ /N$^2$"]
magn        = [False, False, True]
heat        = [False, True, False]
Cmax        = [2, 2, 6]
name        = ["U", "C", "M"]

k = 0
for N in N_arr:
    data = dataset(N)
    for i in range(3):
        plot_data(data[i], title_arr[i].replace("X", str(N)), "T", y_arr[i], magn[i], heat[i], Cmax[k], name[i])
    k += 1







