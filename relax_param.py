from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.constants import c
from lmfit import Model
#
def KINGDEF(rr, r_c):#A, r_c):
    A=I0
    return  A / (1 + (rr/r_c)**2)
def KING2D(x,y,theta,e,r_c):#A,r_c):
    A=I0
    return A / (1 + (((x*np.cos(theta)+y*np.sin(theta))**2)+((e**2)*((-x*np.cos(theta)+y*np.cos(theta))**2)))/((r_c)**2))
#
cosmo = FlatLambdaCDM(H0=68.3,Om0=0.299) #Bocquet et al. 2015 cosmology
c = c.to(u.km/u.s).value
#

cluster = '2344-4224'    #cluster ID
zcl = 0.282384    #cluster redshift
ra, dec = 356.1481, -42.41    #SZ-center in deg
mag = 18.5887    #m* i-band

rel_lst = []
crt = Table.read("SPT-CLJ"+str(cluster)+"_redsequence.des",format="ascii")
magtxt = "MOF_BDF_MAG_I_CORRECTED"
###############CLIP##################
d = SkyCoord(ra*u.degree,dec*u.degree) 
catalog = SkyCoord(ra=crt["RA"]*u.degree,dec=crt["DEC"]*u.degree)
lst = d.separation(catalog).to(u.arcsec)
d_A = cosmo.angular_diameter_distance(z=zcl)
distance_Mpc = (lst*d_A).to(u.Mpc,u.dimensionless_angles())
crt["PRJ_SEP"] = np.array(distance_Mpc)
crt["POS_ANG"] = np.array(d.position_angle(SkyCoord(crt["RA"]*u.degree,crt["DEC"]*u.degree)))
clip = crt[crt["PRJ_SEP"]<=2]      
clip = clip[abs(clip[magtxt]-mag)<=2]      #instead of photo-z cut, just cut in magnitudes
################ABSMAGS###############
D = (c*zcl)/cosmo.H0.value 
clip["ABSMAG"] = clip[magtxt]-(5*np.log10(D*(10**5)))
clip["flux"] = (10**(clip["ABSMAG"]/(-2.5)))*(1.8646e-8) 
absmag = mag - (5*np.log10(D*(10**5))) 
clflux = (10**(absmag/(-2.5)))*(1.8646e-8)
####################################### 
clipp = clip[clip["PRJ_SEP"]<=1]
l1mpc = np.sum(clipp["flux"] - np.min(clipp["flux"]))    
r200 = 10**(-0.57 + 0.44*np.log10(l1mpc/clflux)) 
siga = 0.03                                   #values from wen & han 2013
sigb = 0.15                              
sigtot = (siga + sigb*(clip["PRJ_SEP"]/r200))*r200 
clip["sigtot"] = sigtot
clip["x"] = clip["PRJ_SEP"]*np.sin(clip["POS_ANG"])
clip["y"] = clip["PRJ_SEP"]*np.cos(clip["POS_ANG"])
h = plt.hist2d(clip["x"],clip["y"],bins=(200,200),range=[(-2,2),(-2,2)])
smomap = h[0]        #smooth optical map
for i in range(len(h[0][:,0])):
    for j in range(len(h[0])):
        smomap[i,j] = np.sum(clip["flux"]*(((2*np.pi*(clip["sigtot"]**2))**(-1))*np.exp((-1)*(((h[1][j]+0.01-clip["x"])**2)+((h[2][-(i+1)]-0.01-clip["y"])**2))/(2*((clip["sigtot"])**2)))))
smomap = smomap - np.min(smomap)                           
lst = []
lst1 = []
for i in range(len(h[0][:,0])):
    for j in range(len(h[0])):
        R = np.sqrt(((h[1][j]+0.01)**2) + ((h[2][-(i+1)]-0.01)**2))  
        if R <= (2/3)*r200:
            lst += [smomap[i,j]**2]                       #smomap have elements of order e-17, is this okey????
            lst1 += [(smomap[i,j] - smomap[-(i+1),-(j+1)])**2]
SS = np.sum(lst) 
alpha = np.sum(lst1)/(2*SS)                             #ASYMMETRY FACTOR!!!!!!
I0 = smomap[99,99]                  #amplitude of king model 
rho = np.linspace(0,199,200)        
kmodel = Model(KINGDEF)             #1D King model (no r_tide or r_tide = np.inf)
c200_lst = []
for phi in np.arange(-175,185,5):       #72 directions, 5 each section
    print("c200_lst completed at "+str(int(np.around(((phi+175)/355)*100)))+"%")   #useful to see progress
    jr = rho*np.cos((np.pi/2)-(phi*np.pi/180)) + 99              #cluster centric radial index component X
    ir = rho*(-np.sin((np.pi/2)-(phi*np.pi/180))) + 99           #cluster centric radial index component Y
    lst = [] 
    for k in range(len(rho)):                                    #match precise [i,j] index without repetition
        for i in range(len(smomap[:,0])):
            for j in range(len(smomap[0])):
                if i>(ir[k]-0.5) and i<(ir[k]+0.5):
                    if j>(jr[k]-0.5) and j<(jr[k]+0.5):
                        if [i,j,smomap[i,j]] in lst:
                            continue 
                        lst += [[i,j,smomap[i,j]]]              #also add bin value
    lprof = np.array(lst)
    rr = []
    for i in range(len(lprof)):               #convert [i,j] to cluster centric distance in mpc
        rr += [np.sqrt((h[1][int(np.around(lprof[i][1]))] +0.01)**2 + (h[2][-(int(np.around(lprof[i][0]))+1)] - 0.01)**2)]
    rr = np.array(rr)/r200             #normalized by r200
    r0_lst = []
    calisto = (I0 - lprof[:,2][1:])>0
    for i in range(len(rr)):
        try:
            r0_lst += [np.sqrt((lprof[i][2]*(rr[i]**2))/(I0-lprof[i][2]))]         #calculate king's r0 parameter for every element in (rr,lprof[:,2])
        except:
            r0_lst += [np.inf] 
    r0_lst = list(np.array(r0_lst)[np.array(r0_lst)>0])
    r0 = np.median(r0_lst[1:])
    if len(calisto[calisto==True])==0:
        r0 = r200+0
    params = kmodel.make_params(A=I0, r_c=r0)                      #initial params for king model
    result = kmodel.fit(lprof[:,2], params, rr=rr)               #fit to data
    r0 = np.array(result.params)[0]
    if r0<0:                       #just to be sure
        print("negative r0!")
    c200_lst += [abs(r200/r0)]             #extract best fit r0 and calculate c200
c200R = np.min(c200_lst)                      #minimum steepness factor
c200m = np.mean(c200_lst)                     #mean steepness factor
angdir = np.arange(-175,185,5)[list(c200_lst).index(np.min(c200_lst))]    #ridge direction
beta = c200R/c200m                                      #RIDGE FLATNESS!!!!!!!!
lst = []
for i in range(len(smomap[:,0])):                    #convert data to (x,y,z) list for data within (2/3)*r200
    for j in range(len(smomap[0])):
        if np.sqrt(np.around(j*0.02 - 1.99,2)**2 + np.around(1.99 - i*0.02,2)**2) > ((2/3)*r200):
            continue 
        lst += [[np.around(j*0.02 - 1.99,2), np.around(1.99 - i*0.02,2),smomap[i,j]]]
lst = np.array(lst)
kmodel2 = Model(KING2D,independent_vars=['x','y'])           #fit 2D king model
params2 = kmodel2.make_params(theta=0, e=1 , A=I0, r_c=r0)    #initial guess params (r0 is the r_c for the 1D king model at degree 180)
result2 = kmodel2.fit(lst[:,2], params2, x=lst[:,0], y=lst[:,1])
z = KING2D(lst[:,0],lst[:,1],np.array(result2.params)[0],np.array(result2.params)[1],abs(np.array(result2.params)[2]))#,abs(np.array(result2.params)[3]))
delta = np.sum((lst[:,2] - z)**2)/SS                    #NORMALIZED DEVIATION!!!!!!!!!!!!
gamma = (beta -(1.9*alpha + 3.58*delta + 0.1))/np.sqrt(1 + 1.9**2 + 3.58**2)      #wen&han2013: A=1.9; B=3.58; C=0.1; k = np.sqrt(1 + A**2 + B**2)
rel_lst += [["J"+cluster,alpha,beta,delta,gamma]]
rel_lst = np.array(rel_lst)
reltab = Table([rel_lst[:,0],rel_lst[:,1],rel_lst[:,2],rel_lst[:,3],rel_lst[:,4]],names=["SPT-CL","alpha","beta","delta","gamma"],dtype=('<U22','float64','float64','float64','float64'))
reltab.write("reltab.csv",format="csv")
