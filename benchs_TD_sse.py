import os


    #if __main__:

deltas = [0.1,1,10,100,1000]
iters = [10,9,9,9,9]


total_t=0
cpt=0

for d,k in zip(deltas,iters):
    for i in range(k):
        total_t += d
        cmd=""
        if cpt==0:
            cmd = "ibexopt photovoltaic_TD_sse_.mbx --trace -t "+str(d)+" > res_TD_sse/res_"+str(total_t)+".txt"
        elif cpt<30:
            cmd = "ibexopt photovoltaic_TD_sse_.mbx --trace -t "+str(d)+" -i photovoltaic_TD_sse_.cov > res_TD_sse/res_"+str(total_t)+".txt"
        else:
            cmd = "ibexopt photovoltaic_TD_sse_.mbx -t "+str(d)+" -i photovoltaic_TD_sse_.cov > res_TD_sse/res_"+str(total_t)+".txt"
        print("running: ",cmd)
        os.system(cmd)
        cpt+=1

print(cpt, "runs executed!")
