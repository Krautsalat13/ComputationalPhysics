import numpy as np
import matplotlib.pyplot as plt
import copy
import time
import scipy as sc
import scipy.integrate as int_

start = time.time()

plt.rcParams["text.usetex"] = False
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"]   = "18"

rand = np.random
rand.seed(1556)

per_bounds = False

def setup(title, xlabel, ylabel):
    plt.figure(figsize=(10,45/8))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks()
    plt.yticks()
    plt.locator_params(nbins=6)
    plt.grid()

def UN_theo(N, T):
    def k(T):
        return 1/np.sinh(2/T)**2

    def integrand(k,x):
        return 1/np.sqrt(1 - 4*k/(1+k)**2 * np.sin(x)**2)

    def U(T):
        I = int_.quad(lambda x: integrand(k(T), x), 0, np.pi/2)[0]
        U = -1/np.tanh(2/T)*( 1 + 2/np.pi *(2*np.tanh(2/T)**2 -1 )*I )
        return U
    return np.vectorize(U)(T)

def CN_theo(N, T):
    return (N-1)/N *(1/(T*np.cosh(1/T)))**2

def MN_theo(N, T):
    T_c = 2/(np.log(1+np.sqrt(2)))
    return np.heaviside(T_c-T , 0)*np.abs((1-1/(np.sinh(2/T))**4))**(1/8)


def energy(S):
    N = len(S[0])
    S1 = S[1:,]                # exclude the first row
    S2 = S[:-1,]               # exclude the last row
    
    S3 = S[:,1:]               # exclude the first coloumn
    S4 = S[:,:-1]              # exclude the last coloumn
    
    per_bound_contr = 0
    if per_bounds:
        S5 = S[0,]
        S6 = S[N-1,]
        
        S7 = S[:,0]
        S8 = S[:,N-1]
        
        per_bound_contr = np.sum(S5*S6) + np.sum(S7*S8)
    
    # Now we can express the sums via
    E = - np.sum(S1*S2) - np.sum(S3*S4) - per_bound_contr
    return E


def energy_diff(S, N, i, j):
    if not per_bounds:
        first_term  = S[i-1,j] if i > 0 else 0
        second_term = S[i+1,j] if i < N-1 else 0 
        third_term  = S[i,j-1] if j > 0 else 0 
        fourth_term = S[i,j+1] if j < N-1 else 0  
    else:
        first_term  = S[i-1,j]
        second_term = S[i+1,j] if i < N-1 else S[0,j] 
        third_term  = S[i,j-1]
        fourth_term = S[i,j+1] if j < N-1 else S[i,0] 
        
    dE = 2*S[i,j]*(first_term + second_term + third_term + fourth_term)
    return dE

def spin_flip(S, N, T):
    i, j = (rand.random(size=2)*N).astype(int)
    r = rand.random()

    dE = energy_diff(S, N, i, j)
    q  = np.exp(-dE/T)
    if q > r:
        S[i, j] *= -1

    return S

    

def find_groundstate(T, m, N):
    S = rand.choice([-1, 1], size=(N, N))
    for i in range(m):
        if i%int(m/10) == 0:
            print(np.sum(S))
            # plt.imshow(S)
            # plt.show()
            print(i/m)
        S       = spin_flip(S, N, T)
    
    plt.imshow(S)


#find_groundstate(3, 100000, 50)





def calc_U_C(N_samples, N):
    K       = 20
    N_wait  = int(1e5)              # N_samples
    N_run   = int(N_samples*N**2)        # int(N_samples*N)
    
    T   = np.linspace(0.2, 4, K)
    UN  = np.zeros(K)
    CN  = np.zeros(K)
    MN  = np.zeros(K)

    # initial field
    S = rand.choice([-1, 1], size=(N,N))
    for j in range(N_wait):
        S = spin_flip(S, N, T[K-1])         # reach thermal eq. in N_wait steps
    print("passed N_wait-Loop")
    for k in range(K):
        #print(T[K-1-k])
        E       = energy(S)
        M_temp  = np.sum(S)/N**2

        dE      = 0
        dM      = 0
        dM_temp    = 0
        U_temp  = 0
        C_temp  = 0
        for i in range(N_run):
            E           += dE
            U_temp      += E/N_run
            C_temp      += E**2/N_run
            dM_temp        += dM/N**2
            #M_temp      += np.abs(np.sum(S)/N**2)/N_run

            i, j    = (rand.random(size=2)*N).astype(int)
            r       = rand.random()
            dE      = energy_diff(S, N, i, j)
            q       = np.exp(-dE/T[K-1-k])
            if q > r:
                S[i, j] *= -1
                dM      = 2*S[i, j]
            else:
                dE = 0
                dM = 0
        

        UN[K-1-k] = U_temp/N**2
        CN[K-1-k] = (C_temp - U_temp**2)/T[K-1-k]**2/N**2
        MN[K-1-k] = np.abs(M_temp + dM_temp)
        print("T, M", T[K-1-k], np.abs(M_temp + dM_t/N_run))
    return UN, CN, MN













    




def gen_data():
    N_arr           = np.array([10])
    #N_arr           = np.array([10, 50, 100])
    N_samples_arr   = np.array([1000, 10000])
    
    def text_rep(a, b, c):
        text = "2D_Data/2D_x_qN_zNS.npy"
        text = text.replace("x", str(a))
        text = text.replace("q", str(b))
        text = text.replace("z", str(c))
        return text

    
    for N in N_arr:
        print("Next step -------------------------------------------")
        for N_s in N_samples_arr:
            UN2, CN2, MN2 = calc_U_C(N_s, N)
            np.save(text_rep("U", N, N_s), UN2)
            np.save(text_rep("C", N, N_s), CN2)
            np.save(text_rep("M", N, N_s), MN2)

def gen_plots():
    show_plots  = True
    T = np.linspace(0.2, 4, 20)
    def load_data(a, b, c):
        text = "2D_Data/2D_x_qN_zNS.npy"
        text = text.replace("x", str(a))
        text = text.replace("q", str(b))
        text = text.replace("z", str(c))
        return np.load(text)

    U_10N_1000NS    = load_data("U", 10, 1000)
    U_10N_10000NS   = load_data("U", 10, 10000)
    U_50N_1000NS    = load_data("U", 50, 1000)
    U_50N_10000NS   = load_data("U", 50, 10000)
    U_100N_1000NS   = load_data("U", 100, 1000)
    U_100N_10000NS  = load_data("U", 100, 10000)
    
    C_10N_1000NS    = load_data("C", 10, 1000)
    C_10N_10000NS   = load_data("C", 10, 10000)
    C_50N_1000NS    = load_data("C", 50, 1000)
    C_50N_10000NS   = load_data("C", 50, 10000)
    C_100N_1000NS   = load_data("C", 100, 1000)
    C_100N_10000NS  = load_data("C", 100, 10000)

    M_10N_1000NS    = load_data("M", 10, 1000)
    M_10N_10000NS   = load_data("M", 10, 10000)
    M_50N_1000NS    = load_data("M", 50, 1000)
    M_50N_10000NS   = load_data("M", 50, 10000)
    M_100N_1000NS   = load_data("M", 100, 1000)
    M_100N_10000NS  = load_data("M", 100, 10000)
    
    U_list = [U_10N_1000NS, U_10N_10000NS, U_50N_1000NS, U_50N_10000NS, U_100N_1000NS, U_100N_10000NS]
    C_list = [C_10N_1000NS, C_10N_10000NS, C_50N_1000NS, C_50N_10000NS, C_100N_1000NS, C_100N_10000NS]
    M_list = np.abs(np.array([M_10N_1000NS, M_10N_10000NS, M_50N_1000NS, M_50N_10000NS, M_100N_1000NS, M_100N_10000NS]))
    N_list     = [10, 10, 50, 50, 100, 100]
    Ns_list    = [1000, 10000, 1000, 10000, 1000, 10000]
    
    # setup("Inner Energy per Spin", "T", "U/N")
    # plt.plot(T,  UN_theo(10, T), color="black", ls=":", linewidth=3, label="Theory, N=10")
    # plt.plot(T, UN_theo(100, T), color="black", ls="-.",linewidth=3, label="Theory, N=100")
    # plt.plot(T, UN_theo(1000, T), color="black", ls="--",linewidth=3, label="Theory, N=1000")
    # for i in range(6):
    #     plt.plot(T, U_list[i], label=r"N={}, $N_S={}$".format(N_list[i], Ns_list[i]))
    # plt.legend()
    # plt.savefig("2D_Plots/1D_UN.pdf")
    
    # setup("Specific Heat per Spin", "T", "C/N")
    # plt.plot(T,  CN_theo(10, T), color="black", ls=":", linewidth=3, label="Theory, N=10")
    # plt.plot(T, CN_theo(100, T), color="black", ls="-.",linewidth=3, label="Theory, N=100")
    # plt.plot(T, CN_theo(1000, T), color="black", ls="--",linewidth=3, label="Theory, N=1000")
    # for i in range(6):
    #     plt.plot(T, C_list[i], label=r"N={}, $N_S={}$".format(N_list[i], Ns_list[i]))
    # plt.legend()
    # plt.savefig("2D_Plots/1D_CN.pdf")

    setup("Inner Energy per Spin, N=10", "T", "$U/N^2$")
    plt.plot(T,  UN_theo(10, T), color="black", ls="--", linewidth=3, label="Theory, N=10")
    plt.plot(T, U_list[0], label=r"$N_S={}$".format(1000))
    plt.plot(T, U_list[1], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_UN_10.pdf")
    if not show_plots:
        plt.close()

    setup("Inner Energy per Spin, N=50", "T", "$U/N^2$")
    plt.plot(T,  UN_theo(50, T), color="black", ls="--", linewidth=3, label="Theory, N=50")
    plt.plot(T, U_list[2], label=r"$N_S={}$".format(1000))
    plt.plot(T, U_list[3], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_UN_50.pdf")
    if not show_plots:
        plt.close()


    setup("Inner Energy per Spin, N=100", "T", "$U/N^2$")
    plt.plot(T,  UN_theo(100, T), color="black", ls="--", linewidth=3, label="Theory, N=100")
    plt.plot(T, U_list[4], label=r"$N_S={}$".format(1000))
    plt.plot(T, U_list[5], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_UN_100.pdf")
    if not show_plots:
        plt.close()

   
    
    
    setup("Specific Heat per Spin, N=10", "T", "$C/N^2$")
    plt.plot(T,  CN_theo(10, T), color="black", ls="--", linewidth=3, label="Theory, N=10")
    plt.plot(T, C_list[0], label=r"$N_S={}$".format(1000))
    plt.plot(T, C_list[1], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_CN_10.pdf")
    if not show_plots:
        plt.close()


    setup("Specific Heat per Spin, N=50", "T", "$C/N^2$")
    plt.plot(T,  CN_theo(50, T), color="black", ls="--", linewidth=3, label="Theory, N=50")
    plt.plot(T, C_list[2], label=r"$N_S={}$".format(1000))
    plt.plot(T, C_list[3], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_CN_50.pdf")
    if not show_plots:
        plt.close()


    setup("Specific Heat per Spin, N=100", "T", "$C/N^2$")
    plt.plot(T,  CN_theo(100, T), color="black", ls="--", linewidth=3, label="Theory, N=100")
    plt.plot(T, C_list[4], label=r"$N_S={}$".format(1000))
    plt.plot(T, C_list[5], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_CN_100.pdf")
    if not show_plots:
        plt.close()


    
    
    setup("absolute Magnetization per Spin, N=10", "T", "$M/N^2$")
    #plt.plot(T,  MN_theo(10, T), color="black", ls="--", linewidth=3, label="Theory, N=10")
    plt.plot(T, M_list[0], label=r"$N_S={}$".format(1000))
    plt.plot(T, M_list[1], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_MN_10.pdf")
    if not show_plots:
        plt.close()


    setup("absolute Magnetization per Spin, N=50", "T", "$M/N^2$")
    plt.plot(T,  MN_theo(50, T), color="black", ls="--", linewidth=3, label="Theory, N=50")
    plt.plot(T, M_list[2], label=r"$N_S={}$".format(1000))
    plt.plot(T, M_list[3], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_MN_50.pdf")
    if not show_plots:
        plt.close()


    setup("absolute Magnetization per Spin, N=100", "T", "$M/N^2$")
    plt.plot(T,  MN_theo(100, T), color="black", ls="--", linewidth=3, label="Theory, N=100")
    plt.plot(T, M_list[4], label=r"$N_S={}$".format(1000))
    plt.plot(T, M_list[5], label=r"$N_S={}$".format(10000))
    plt.legend()
    plt.savefig("2D_Plots/2D_MN_100.pdf")
    if not show_plots:
        plt.close()



gen_data()
# gen_plots()


end = time.time()
print(end-start)