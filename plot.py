# T. C.-G.	Feb 22nd 2018
# start with:	python plot.py
# to produce figures for all OCs, or with:	python plot.py  Alessi_2 ASCC_113 IC_4756
# to produce figures for chosen OCs only.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.cm.get_cmap('viridis_r')
import sys
from astropy.io import fits

hdul = fits.open('CG18_OC_members_light.fits')
ra = hdul[1].data['ra']
dec = hdul[1].data['dec']
pmra = hdul[1].data['pmra']
epmra = hdul[1].data['pmra_error']
pmdec = hdul[1].data['pmdec']
epmdec = hdul[1].data['pmdec_error']
parallax = hdul[1].data['parallax']
epar = hdul[1].data['parallax_error']
pmrapmdeccorr =  hdul[1].data['pmra_pmdec_corr']
g = hdul[1].data['G']
jmag = hdul[1].data['Jmag']
kmag = hdul[1].data['Kmag']
proba = hdul[1].data['finalProba']
ocname = hdul[1].data['ocname']
pmused = hdul[1].data['PMused']



for oc in sorted(set(ocname)):

	if len(sys.argv)>1 and (oc not in sys.argv):
		continue

	print 'Now generating figure for',oc

	# Subset data for each individual OC:
	RA = ra[ocname==oc]
	DEC = dec[ocname==oc]
	PMRA = pmra[ocname==oc]
	ePMRA = epmra[ocname==oc]
	PMDEC = pmdec[ocname==oc]
	ePMDEC = epmdec[ocname==oc]
	PAR = parallax[ocname==oc]
	ePAR = epar[ocname==oc]
	PMRAPMDECCORR = pmrapmdeccorr[ocname==oc]
	J = jmag[ocname==oc]
	K = kmag[ocname==oc]
	G = g[ocname==oc]

	PROBA = proba[ocname==oc]

	fig=plt.figure(0,figsize=(14,9)); plt.clf()

	#	RA DEC
	plt.subplot(231); plt.title(oc)
	plt.scatter( RA , DEC , c=PROBA , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap)
	plt.scatter( RA[PROBA>0.1] , DEC[PROBA>0.1] , c=PROBA[PROBA>0.1] , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap , zorder=10)
	plt.xlabel('ra [deg]'); plt.ylabel('dec [deg]')
	plt.xlim( min(RA) , max(RA) ); plt.ylim( min(DEC) , max(DEC) )
	plt.minorticks_on()

	#	PMRA PMDEC
	plt.subplot(232)
	if pmused[ocname==oc][0]=='T':
		plt.title('TGAS proper motions\nwere used for membership')
	else:
		plt.title('UCAC4 proper motions\nwere used for membership')


	plt.scatter( PMRA , PMDEC , c=PROBA , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap)
	for i in range(len(PMRA)):
		if PROBA[i]>0.1:
			plt.errorbar( PMRA[i] , PMDEC[i] , xerr=ePMRA[i] , yerr=ePMDEC[i] , fmt=',', color=cmap(PROBA[i]) , zorder=1 )
	plt.xlabel('TGAS pmra [mas/yr]'); plt.ylabel('TGAS pmdec [mas/yr]')
	plt.xlim( min(PMRA[PROBA>0]) - 5 , max(PMRA[PROBA>0]) + 5); plt.ylim( min(PMDEC[PROBA>0]) - 5 , max(PMDEC[PROBA>0]) + 5)
	plt.minorticks_on()




	#	PMRA PMDEC with ellipses showing correlations
	def a_b_theta_from_cov(covariancematrix):
		"""
		Parameters of the ellipse representing a covariance matrix.
		Input: covariance matrix 2x2
		Returns a,b,theta
			a = semimajor axis
			b = semiminor axis
			theta = rotation of major axis with respect to x (in degrees)
		"""
		eigenvalues,eigenvectors = np.linalg.eigh(covariancematrix)
		a,b = sorted([np.sqrt(v) for v in eigenvalues])[::-1]
		maxvector = [v for v,val in zip(eigenvectors,eigenvalues) if val==max(eigenvalues)][0]
		theta=np.degrees( np.arctan( maxvector[1] / maxvector[0] ) )
		return a,b,theta
	ax=fig.add_subplot(233)
	plt.xlabel('TGAS pmra [mas/yr]'); plt.ylabel('TGAS pmdec [mas/yr]')
	plt.xlim( min(PMRA[PROBA>0]) - 5 , max(PMRA[PROBA>0]) + 5); plt.ylim( min(PMDEC[PROBA>0]) - 5 , max(PMDEC[PROBA>0]) + 5)
	plt.minorticks_on()
	#build covariance matrix and plot corresponding ellipse for each point
	for i in range(len(PMRA)):
		#build ellipse:
		covmat=np.array( [[ePMRA[i]**2, PMRAPMDECCORR[i]*ePMRA[i]*ePMDEC[i]],
			      [PMRAPMDECCORR[i]*ePMRA[i]*ePMDEC[i], ePMDEC[i]**2]] )
		a,b,tilt = a_b_theta_from_cov(covmat)
		ell = matplotlib.patches.Ellipse(xy=[PMRA[i],PMDEC[i]],
			                         width=2*a,
			                         height=2*b,
			                         angle = tilt,
			                         facecolor='none',
			                         edgecolor=cmap(PROBA[i]))
		if PROBA[i]>0.1:
			ax.add_patch(ell)



	#
	plt.subplot(234)
	plt.scatter( G , PAR , c=PROBA , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap)
	for i in range(len(PMRA)):
		if PROBA[i]>0.1:
			plt.errorbar( G[i] , PAR[i] , yerr=ePAR[i] , fmt=',', color=cmap(PROBA[i]) , zorder=1 )
	plt.xlabel('G'); plt.ylabel('parallax [mas]')
	plt.ylim( min(PAR[PROBA>0]) - 1 , max(PAR[PROBA>0]) + 1)
	plt.minorticks_on()






	#	J vs JK
	plt.subplot(235)
	plt.scatter( J-K , J , c=PROBA , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap , zorder=1)
	plt.scatter( J[PROBA>0.1]-K[PROBA>0.1] , J[PROBA>0.1] , c=PROBA[PROBA>0.1] , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap , zorder=10)
	plt.xlabel('J - K'); plt.ylabel('J')
	plt.ylim( max(J) , min(J) )
	plt.minorticks_on()

	#	G vs GK
	plt.subplot(236)
	plt.scatter( G-K , G , c=PROBA , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap)
	plt.scatter( G[PROBA>0.1]-K[PROBA>0.1] , G[PROBA>0.1] , c=PROBA[PROBA>0.1] , edgecolor='none' , vmin=0 , vmax=1 , cmap=cmap , zorder=10)
	plt.xlabel('G - K'); plt.ylabel('G')
	plt.ylim( max(G) , min(G) )
	plt.minorticks_on()



	plt.tight_layout()




	plt.savefig( oc+'_6panels.png' )
	print '   figure saved as',oc+'_6panels.png'

