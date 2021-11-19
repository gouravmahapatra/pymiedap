print('''
 ## Testing disk integration

 Let's define a model that corresponds to a Lambertian surface with surface
 albedo 1.0. \n
 We will compare this with the analytical solution.
 ''')

modelA = pmd.Model()
modelA.wvl_list = [0.7]
del modelA.layers.gastop
del modelA.layers.haze
del modelA.layers.cloud
modelA.layers.gasbelow.rayscat=False
modelA.layers.gasbelow.tau=[0.0]
modelA.surface[0,0] = 1.0

alphas = np.linspace(0,np.pi,80) # phase angles
alphas_deg = np.degrees(alphas)  # in degrees
theta = np.pi - alphas  # scattering angle
P = 2*(np.sin(theta) - theta*np.cos(theta))/(3.*np.pi)  # analytical Lambertian phase function

pmd.planet_integrated([modelA],npix=100, alpha=alphas_deg, force=True)

mpl.plot(modelA.phase, modelA.I[0,:]-P)  #comparing the two models
mpl.xlabel('Phase angle')
mpl.ylabel('Delta I')
mpl.title('Delta between PyMieDAP and lambertian analytic solution')
mpl.show()

