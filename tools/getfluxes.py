import harm_script2

ts = []
fs = []
md = []
#jem = []
#jtot = []

starti,endi=0,1000
for dumpno in range(starti,endi+1):
    rd("dump%03d" % dumpno);

    rhor = 1+(1-a**2)**0.5
    ihor = iofr(rhor)
    #cvel()
    Tcalcud()
    ts.append(int(t)) #time
    fs.append(horfluxcalc(ihor)) #horizon flux
    md.append(mdotcalc(ihor)) #mass accretion rate
    #jem.append(jetpowcalc(0)[ihor]) #jet power (EM)
    #jtot.append(jetpowcalc(2)[ihor]) #jet power (total)

#write to file
fluxes_file = open("fluxes.dat","w+")
for line in range(len(ts)):
    fluxes_file.write(str(ts[line]))
    fluxes_file.write("\t")
    fluxes_file.write(str(fs[line]))
    fluxes_file.write("\t")
    fluxes_file.write(str(md[line]))
    fluxes_file.write("\n")
fluxes_file.close()
