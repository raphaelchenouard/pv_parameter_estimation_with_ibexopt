import matplotlib.pyplot as plt
import numpy as np

deltas = [0.1,1,10,100,1000,2000]
iters = [10,9,9,9,9,9]

SD_loup = np.array([7.49129827085,0.0688909434095,0.053910116311,0.04518650551130,0.000536542997341,8.23551351815e-05,7.67702303276e-05,7.67702303276e-05,5.51702697187e-05,3.69774524955e-05,3.52921423861e-05,2.73231853554e-05,2.66778670569e-05,2.65046406118e-05,2.61368858309e-05,2.61368858309e-05,2.61368858309e-05,2.57687368363e-05,2.57687368363e-05,2.57023931359e-05,2.53504066716e-05,2.53504066716e-05,2.53504066716e-05,2.53504066716e-05,2.53504066716e-05,2.53296457074e-05,2.53296457074e-05,2.53296457074e-05,2.53146561963e-05,2.53146561963e-05,2.53146561963e-05,2.53146561963e-05,2.53146561963e-05,2.53146561963e-05,2.53146561963e-05,2.53121089848e-05,2.53018531357e-05,2.53018531357e-05,2.53018531357e-05,2.53018531357e-05,2.53018531357e-05,2.52884613948e-05,2.52839149171e-05,2.52839149171e-05,2.52790113993e-05,2.52790113993e-05,2.52790113993e-05,2.52789893119e-05,2.52787227298e-05,2.52785550593e-05])#,2.52785550593e-05])
SD_uplo = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.60159146222e-08,6.02091383024e-07,1.51176707396e-06,2.00477208777e-06,2.318776961e-06,2.56023931369e-06,2.56687368373e-06,2.60368858319e-06,2.66207230263e-06,4.31181902286e-06,5.15885715587e-06,5.9696790191e-06,6.78453646778e-06,7.23298903419e-06,7.75556526765e-06,8.29050055439e-06,8.75607588865e-06,9.1387132271e-06,1.19723043038e-05,1.42349590767e-05,1.61244327932e-05,1.76956590401e-05,1.91312228016e-05,2.04740300593e-05,2.16165768112e-05,2.25662200665e-05,2.32687316902e-05,2.43478714234e-05,2.48524762804e-05])#,2.51785550593e-05])

DD_loup = np.array([2.0460347509,0.287728191826,0.287728191826,0.287728191826,0.287728191826,0.287728191826,0.287728191826,0.088526092054,0.088526092054,0.088526092054,0.088526092054,0.00575339357769,8.37848995448e-05,4.52317355276e-05,2.75446635931e-05,2.75446635931e-05,2.75446635931e-05,2.75446635931e-05,2.7029319583e-05,2.7029319583e-05,2.7029319583e-05,2.7029319583e-05,2.60290409271e-05,2.60290409271e-05,2.60290409271e-05,2.60290409271e-05,2.60290409271e-05,2.57403337259e-05,2.57403337259e-05,2.57403337259e-05,2.57403337259e-05,2.57403337259e-05,2.57403337259e-05,2.57403337259e-05,2.53248262719e-05,2.53163141702e-05,2.53163141702e-05,2.53163141702e-05,2.53163141702e-05,2.53163141702e-05,2.53163141702e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05,2.52038080903e-05])
DD_uplo = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

TD_loup = np.array([7.49129827085,0.309810003575,0.0498281923014,0.0498281923014,0.0362774025583,0.0362774025583,0.0362774025583,0.0353688780919,0.0353688780919,0.0345984641854,0.0345984641854,0.0345640532744,0.0289599345759,0.000224064725986,7.66611763393e-05,7.53214728001e-05,7.53214728001e-05,5.46176545474e-05,4.32107601887e-05,4.32107601887e-05,4.2107898703e-05,4.2107898703e-05,3.16213057393e-05,3.02111210549e-05,3.02111210549e-05,2.90048715513e-05,2.90048715513e-05,2.90048715513e-05,2.90048715513e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.5353861476e-05,2.52466925623e-05,2.51463274713e-05,2.51463274713e-05,2.51463274713e-05,2.51463274713e-05,2.51463274713e-05,2.51463274713e-05,2.51463274713e-05])
TD_uplo = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])


total_t=0
cpt=0

XX=[0]

for d,k in zip(deltas,iters):
    for i in range(k):
        total_t += d
        XX.append(total_t)


print(XX)

N=26

fig1, ax1 = plt.subplots()

ax1.plot(XX[:len(SD_loup)], np.sqrt(SD_loup/N), label="SD_LOUP")

#ax.plot(XX, np.sqrt([2.52783241285e-05/N for i in range(len(XX))]), 'r:', label="Optimal SD solution")
#ax.plot(2.51783241285e-05,2.52783241285e-05,'r.',label="Optimal SD solution")


print('RMSE_SD = '+str(np.sqrt(2.52783241285e-05/N)))
#ax.text(0.1,1.1e-3,'RMSE_SD = ['+str(round(np.sqrt(2.51783241285e-05/N),10))+','+str(round(np.sqrt(2.52783241285e-5/N),10))+']')
RMSE_SD_L = np.sqrt(2.51785550593e-05/N)
RMSE_SD_U = np.sqrt(2.52785550593e-05/N)
ax1.plot(15749.0387911,RMSE_SD_U,'ro',label='SD optimal LOUP')
#ax1.text(800,RMSE_SD_U-4e-4,'['+str(round(RMSE_SD_L,7))+','+str(round(RMSE_SD_U,8))+']')
print('RMSE_SD in [',RMSE_SD_L,',',RMSE_SD_U,'] ~ ',RMSE_SD_U-RMSE_SD_L)

ax1.plot(XX[:len(DD_loup)], np.sqrt(DD_loup/N), label="DD_LOUP")
ax1.plot(XX[:len(TD_loup)], np.sqrt(TD_loup/N), label="TD_LOUP")

ax1.plot(XX[:len(SD_uplo)], np.sqrt(SD_uplo/N), "--", label="SD_UPLO")
ax1.plot(XX[:len(TD_uplo)], np.sqrt(TD_uplo/N), "--", label="TD_UPLO")
ax1.plot(XX[:len(DD_uplo)], np.sqrt(DD_uplo/N), "--", label="DD_UPLO")


ax1.set_yscale("log")
ax1.set_xscale("symlog")

ax1.set(xlabel='time (s)', ylabel='RMSE value')
ax1.grid()
ax1.legend()
#plt.title('SSE convergence curve')
plt.savefig('SSE_conv.pdf',format='pdf')


#SD_box=[80,76,84,74,68,64,84,84,72,84,952,962,892,926,812,966,778,936,940,9736,10506,11156,11128,7918,9414,10570,10282,10126,111610,107368,111482,112692,107474,109026,110718,110138,110278,1149802,1191294,1242660,1315048,1417636,1570792,1683048,1725940,1737196,3449632,3534670,4146330,0]
#DD_box=[46,38,38,38,28,40,48,42,42,38,418,464,406,402,432,460,442,476,364,5510,4736,5328,4386,4846,4382,5014,6192,6032,60892,58194,60908,63404,63240,66558,62596,62852,63114,622804,687014,583162,680140,464042,638838,691646,692020,671688,1340072,1265208,1218100,1390972,1346102]
#TD_box=[32,30,32,32,28,28,32,30,32,32,426,440,264,294,266,284,332,336,300,2832,3294,3690,3490,3950,3106,3030,3202,3284,38816,39242,44148,40204,41964,39792,39764,43282,40734,455552,429776,426142,458486,387364,482478,415578,483748,446726,899028,839372]
#ax2.plot(XX[1:len(SD_box)+1], SD_box, label="SD")
#ax2.plot(XX[1:len(DD_box)+1], DD_box, label="DD")
#ax2.plot(XX[1:len(TD_box)+1], TD_box, label="TD")

SD_loup = [0.0869310826398,0.0869310826398,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.07925358208,0.00432321234519,0.00198208946697,0.0019182659076,0.00190240188769,0.00190240188769,0.00190240188769,0.00190240188769,0.00190240188769,0.00190240188769,0.00190008730711]
SD_uplo = [0,0,0,0,0,0.000343389914367,0.00034526216672,0.000376293793676,0.000377937184892,0.000378590491925,0.000385200624438,0.000413629011101,0.000452021819298,0.00069792938874,0.000804709782787,0.000976675840751,0.00101910174918,0.00109433581968,0.00118186222047,0.00126165033065,0.0018981872198]

# starts at 0.3
DD_loup = [0.999994167365,0.998988351206,0.113078927255,0.113078927255,0.113078927255,0.113078927255,0.0851363309275,0.0851363309275,0.0851363309275,0.0851363309275,0.0797767741004,0.0797767741004,0.0797767741004,0.0797767741004,0.0797767741004,0.0797767741004,0.0797767741004,0.00202958166117,0.00194243569215,0.00191926775785,0.00191926775785,0.00191926775785,0.00191209589417,0.00191209589417,0.00191209589417,0.00191209589417,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190585687592,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190352110007,0.00190151377302,0.00190151377302,0.00190151377302,0.00190151377302]
DD_uplo = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000420226275351,0.000433819295267,0.000447807839937,0.000457026070009,0.000467220592561,0.000470254678901,0.000477165227466,0.000497323693811,0.000534614253723,0.0006151261417,0.000618664477267,0.000624851616253,0.000639791750181,0.000649814632471,0.00065782336784,0.00066697761993,0.000676191554308,0.000687442499833,0.000777605272039,0.000804538275396,0.000828635962457,0.000855293705919,0.000873191915995,0.000887565460266,0.00090532290981,0.000918674118158,0.000930589635184,0.000952452021469,0.000969645875433,0.000985438774667,0.000994011649315,0.00100068579473]

# starts at 0.2
TD_loup = [0.47928104325,0.47928104325,0.47928104325,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.137158218898,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00467296208204,0.00216076238125,0.00203141860379,0.00203141860379,0.00203141860379,0.00203141860379,0.00191957122767,0.00191957122767,0.00191957122767,0.00191957122767,0.00191423587299,0.00191423587299,0.0019088982987,0.00190677451761,0.00190677451761,0.00190677451761,0.00190677451761,0.00190677451761,0.00190677451761,0.00190677451761,0.00190442030934,0.00190442030934,0.00190442030934,0.00190442030934]
TD_uplo = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000222837634007,0.000252716330587,0.000346081869924,0.000349794455529,0.000350320052998,0.00035410139487,0.000364712526365,0.000377953444029,0.000420013634102,0.000424616265709,0.0004393698633,0.000440835300469,0.000444438830635,0.00044579586543,0.000446649170458,0.000446844666049,0.00044833898459,0.000603387109368,0.000616993510558,0.000617866040973,0.000618254805235,0.000618525794015,0.000618734029864,0.000618902431682,0.000619061441043,0.000619204077023,0.000619494614123,0.000619829831053,0.000620164799543,0.000620598456337,0.000621293378784]

fig2, ax2 = plt.subplots()

ax2.plot(XX[:len(SD_loup)], SD_loup, label="SD_LOUP")
ax2.plot(20,0.00190008730711,'ro',label='SD optimal LOUP')
ax2.plot(XX[4:len(DD_loup)+4], DD_loup, label="DD_LOUP")
ax2.plot(XX[3:len(TD_loup)+3], TD_loup, label="TD_LOUP")
ax2.plot(XX[:len(SD_uplo)], SD_uplo, "--", label="SD_UPLO")
ax2.plot(XX[4:len(DD_uplo)+4], DD_uplo, "--", label="DD_UPLO")
ax2.plot(XX[3:len(TD_uplo)+3], TD_uplo, "--", label="TD_UPLO")


ax2.set_xscale("symlog")
ax2.set_yscale("log")
ax2.set(xlabel='time (s)', ylabel='MAE value')
ax2.grid()
ax2.legend()
#plt.title('MAE convergence curve')
plt.savefig('MAE_conv.pdf',format='pdf')

plt.show()