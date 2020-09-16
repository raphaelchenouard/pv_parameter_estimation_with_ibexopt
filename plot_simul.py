import numpy as np
from scipy.optimize import minimize,least_squares,newton
import matplotlib.pyplot as plt

I_exp = [0.764, 0.762, 0.7605, 0.7605, 0.76, 0.759, 0.757, 0.757, 0.7555, 0.754, 0.7505, 0.7465, 0.7385, 0.728, 0.7065, 0.6755, 0.632, 0.573, 0.499, 0.413, 0.3165, 0.212, 0.1035, -0.01, -0.123, -0.21]
V_exp = [-0.2057, -0.1291, -0.0588, 0.0057, 0.0646, 0.1185, 0.1678, 0.2132, 0.2545, 0.2924, 0.3269, 0.3585, 0.3873, 0.4137, 0.4373, 0.459, 0.4784, 0.496, 0.5119, 0.5265, 0.5398, 0.5521, 0.5633, 0.5736, 0.5833, 0.59]



#
#def rmse_sdm(Iph,Isd,Rs,Rsh,n):
#    #return (sum([(Iph - Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T_sdm))-1)-(VL+Rs*IL)/Rsh-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])/N)**0.5
#    return (sum([(model(Iph,Isd,Rs,Rsh,n,VL,IL)-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])/N)**0.5
#
#def sse_sdm_norm(Iph,Isd,Rs,Rsh,n):
#    #return (sum([(Iph - 1e-6*Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T_sdm))-1)-(VL+Rs*IL)/(1e2*Rsh)-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])/N)**0.5
#    return sum([(model_norm(Iph,Isd,Rs,Rsh,n,VL,IL)-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])
#
#def rmse_sdm_norm(Iph,Isd,Rs,Rsh,n):
#    #return (sum([(Iph - 1e-6*Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T_sdm))-1)-(VL+Rs*IL)/(1e2*Rsh)-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])/N)**0.5
#    return (sum([(model_norm(Iph,Isd,Rs,Rsh,n,VL,IL)-IL)**2 for (VL,IL) in zip(V_exp,I_exp)])/N)**0.5
#def error_max_sdm(Iph,Isd,Rs,Rsh,n):
#    #return max([np.abs(Iph - Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T_sdm))-1)-(VL+Rs*IL)/Rsh-IL) for (VL,IL) in zip(V_exp,I_exp)])
#    return max([abs(model(Iph,Isd,Rs,Rsh,n,VL,IL)-IL) for (VL,IL) in zip(V_exp,I_exp)])
#
#def error_max_sdm_norm(Iph,Isd,Rs,Rsh,n):
#    #return max([np.abs(Iph - 1e-6*Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T_sdm))-1)-(VL+Rs*IL)/(1e2*Rsh)-IL) for (VL,IL) in zip(V_exp,I_exp)])
#    return max([abs(model_norm(Iph,Isd,Rs,Rsh,n,VL,IL)-IL) for (VL,IL) in zip(V_exp,I_exp)])



N=len(I_exp) # 26
T = 33+273.0
k = 1.3806503*1e-23
q = 1.60217646*1e-19

Temps = [15+273,25+273,33+273,45+273,55+273]

V_plot = np.linspace(min(V_exp),max(V_exp),50)

I_plots = []
P_plots = []

P_sdm=[]
for (Vi,Ii) in zip(V_exp,I_exp):
    P_sdm.append(Vi*Ii)

# SD
def model_SD(Iph,Isd,Rs,Rsh,n,VL,IL):
    return Iph - 1e-6*Isd*(np.exp(q*(VL+Rs*IL)/(n*k*T))-1)-(VL+Rs*IL)/(1e2*Rsh)


Iph,Isd,Rs,Rsh,n = (0.760771090719 , 0.323459337407 , 0.0363725593061 , 0.537981707126 , 1.48204583019)


for T in Temps:

    I_plot = []
    for VL in V_plot:
        I_plot.append(newton((lambda IL: model_SD(Iph,Isd,Rs,Rsh,n,VL,IL)-IL),(0)))

    I_plots.append(I_plot)
    
    P_plot = []
    for (Vi,Ii) in zip(V_plot,I_plot):
        P_plot.append(Vi*Ii)
    P_plots.append(P_plot)


fig1, ax1_sdm = plt.subplots()
ax1_sdm.plot(V_exp,I_exp,'ro',label='Experimental data with T=33')

#for i,I_p in enumerate(I_plots):
#ax1_sdm.plot(V_plot,I_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax1_sdm.plot(V_plot,I_plots[2],'b-',label='Simulated data with T=33')
ax1_sdm.set_xlabel('Voltage(V)')
ax1_sdm.set_ylabel('Current(A)')
ax1_sdm.legend()
plt.title('I-V curve for SDM')
plt.savefig('IV-SDM.pdf',format='pdf')

fig2, ax2_sdm = plt.subplots()
ax2_sdm.plot(V_exp,P_sdm,'ro',label='Experimental data')
#for i,P_p in enumerate(P_plots):
#ax2_sdm.plot(V_plot,P_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax2_sdm.plot(V_plot,P_plots[2],'b-',label='Simulated data with T=33')

ax2_sdm.set_xlabel('Voltage(V)')
ax2_sdm.set_ylabel('Power(W)')
ax2_sdm.yaxis.set_label_position("right")
ax2_sdm.legend()
plt.title('P-V curve for SDM')
plt.savefig('PV-SDM.pdf',format='pdf')

# DD
def model_DD(Iph,Isd1,Isd2,Rs,Rsh,n1,n2,VL,IL):
    return Iph - 1e-6*Isd1*(np.exp(q*(VL+Rs*IL)/(n1*k*T))-1)\
        - 1e-6*Isd2*(np.exp(q*(VL+Rs*IL)/(n2*k*T))-1)\
        - (VL+Rs*IL)/(1e2*Rsh)

Iph = 0.760734104017
Isd1 = 0.875561685655
Isd2 = 0.189346589121
Rs = 0.0369244615197
Rsh = 0.564497768391
n1 = 1.9530786771
n2 = 1.4375


for T in Temps:
    
    I_plot = []
    for VL in V_plot:
        I_plot.append(newton((lambda IL: model_DD(Iph,Isd1,Isd2,Rs,Rsh,n1,n2,VL,IL)-IL),(0)))
    
    I_plots.append(I_plot)
    
    P_plot = []
    for (Vi,Ii) in zip(V_plot,I_plot):
        P_plot.append(Vi*Ii)
    P_plots.append(P_plot)
fig3, ax3_ddm = plt.subplots()
ax3_ddm.plot(V_exp,I_exp,'ro',label='Experimental data with T=33')

#for i,I_p in enumerate(I_plots):
#ax3_ddm.plot(V_plot,I_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax3_ddm.plot(V_plot,I_plots[2],'b-',label='Simulated data with T=33')
ax3_ddm.set_xlabel('Voltage(V)')
ax3_ddm.set_ylabel('Current(A)')
ax3_ddm.legend()
plt.title('I-V curve for DDM')
plt.savefig('IV-DDM.pdf',format='pdf')

fig4, ax4_ddm = plt.subplots()
ax4_ddm.plot(V_exp,P_sdm,'ro',label='Experimental data')
#for i,P_p in enumerate(P_plots):
#ax4_ddm.plot(V_plot,P_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax4_ddm.plot(V_plot,P_plots[2],'b-',label='Simulated data with T=33')

ax4_ddm.set_xlabel('Voltage(V)')
ax4_ddm.set_ylabel('Power(W)')
ax4_ddm.yaxis.set_label_position("right")
ax4_ddm.legend()
plt.title('P-V curve for DDM')
plt.savefig('PV-DDM.pdf',format='pdf')


# TD
def model_TD(Iph,Isd1,Isd2,Isd3,Rs,Rsh,n1,n2,n3,VL,IL):
    return Iph - 1e-6*Isd1*(np.exp(q*(VL+Rs*IL)/(n1*k*T))-1)\
        - 1e-6*Isd2*(np.exp(q*(VL+Rs*IL)/(n2*k*T))-1)\
        - 1e-6*Isd3*(np.exp(q*(VL+Rs*IL)/(n3*k*T))-1)\
        - (VL+Rs*IL)/(1e2*Rsh)

Iph = 0.760741890645
Isd1 = 0.202577311254
Isd2 = 0.0402347502252
Isd3 = 0.595509723623
Rs = 0.0366887919685
Rsh = 0.557170524344
n1 = 1.45269266801
n2 = 1.49017942153
n3 = 1.9906189751

for T in Temps:
    
    I_plot = []
    for VL in V_plot:
        I_plot.append(newton((lambda IL: model_TD(Iph,Isd1,Isd2,Isd3,Rs,Rsh,n1,n2,n3,VL,IL)-IL),(0)))
    
    I_plots.append(I_plot)
    
    P_plot = []
    for (Vi,Ii) in zip(V_plot,I_plot):
        P_plot.append(Vi*Ii)
    P_plots.append(P_plot)

fig5, ax5_tdm = plt.subplots()
ax5_tdm.plot(V_exp,I_exp,'ro',label='Experimental data with T=33')

#for i,I_p in enumerate(I_plots):
#ax5_tdm.plot(V_plot,I_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax5_tdm.plot(V_plot,I_plots[2],'b-',label='Simulated data with T=33')
ax5_tdm.set_xlabel('Voltage(V)')
ax5_tdm.set_ylabel('Current(A)')
ax5_tdm.legend()
plt.title('I-V curve for TDM')
plt.savefig('IV-TDM.pdf',format='pdf')

fig6, ax6_tdm = plt.subplots()
ax6_tdm.plot(V_exp,P_sdm,'ro',label='Experimental data')
#for i,P_p in enumerate(P_plots):
#ax6_tdm.plot(V_plot,P_p,'b-',label='Simulated data with T='+str(Temps[i]))
ax6_tdm.plot(V_plot,P_plots[2],'b-',label='Simulated data with T=33')

ax6_tdm.set_xlabel('Voltage(V)')
ax6_tdm.set_ylabel('Power(W)')
ax6_tdm.yaxis.set_label_position("right")
ax6_tdm.legend()
plt.title('P-V curve for TDM')
plt.savefig('PV-TDM.pdf',format='pdf')


plt.show()
