import numpy as np
import astroquery.vizier
import astroquery.simbad
import astroquery.gaia
import astroquery.irsa
import astropy
import astropy.units as u

from pathlib import Path
import datetime

astroquery.vizier.Vizier.ROW_LIMIT = -1

# IRSA gives annoying warnings, this will disable all warnings (so should be commented out for debugging)
#import warnings
#warnings.filterwarnings("ignore")

h=6.63e-34
k=1.38e-23
c=3e8

def getInfo(ra=-1,dec=-1,name=-1):
    simbad=astroquery.simbad.Simbad()
    simbad.add_votable_fields('sptype')
    simbad.add_votable_fields('ids')
    simbad.add_votable_fields('otype')
    simbad.add_votable_fields('plx')
    simbad.add_votable_fields('plx_bibcode')
    if (name!=-1):
        star=simbad.query_object(name)[0]
    elif (ra!=-1) & (dec!=-1):
        coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.deg))
        star=simbad.query_region(coord,radius=100*u.arcsecond)
        #finds nearest object with a 2MASS name (in case something spurious is closer to given coords)
        whichEntry=0
        while 1-('2MASS' in star[whichEntry]['IDS'].decode('utf-8')):
            ids=star[whichEntry]['IDS'].decode('utf-8')
            whichEntry+=1
        star=star[whichEntry]
    else:
        print('must enter name or coords (ra, dec) of star')
    
    ra=star['RA']
    dec=star['DEC']
    names=star['IDS'].decode('utf-8')
    startPoint = names.find('2MASS ')+len('2MASS ')
    startName = names[startPoint:]
    endPoint = startName.find('|')
    if endPoint==-1:
        tmName=startName
    else:
        tmName = startName[:endPoint]
    if (tmName[0]!='J') & (name!=-1): #if this named star doesn't have a 2mass name, looks for the nearest one that does
        print('No 2MASS entry for star: ',mainName)
        print('Looking for nearest object with a 2MASS name')
        return getInfo(ra=ra,dec=dec)
    
    starDict={}
    starDict['SimbadName']=star['MAIN_ID'].decode('utf-8') #getting rid of small b...
    starDict['2MASSID']=tmName
    starDict['ObjectType']=star['OTYPE'].decode('utf-8')
    starDict['StellarType']=star['SP_TYPE'].decode('utf-8')
    starDict['StellarTypeSource']=getRef(star['SP_BIBCODE'])
    starDict['RA']=ra
    starDict['DEC']=dec
    starDict['CoordRef']=getRef(star['COO_BIBCODE'])
    starDict['Distance']=1000/star['PLX_VALUE'] #1000 as Simbad gives value in mas
    starDict['DistanceSource']=getRef(star['PLX_BIBCODE'])
    
    Teff,rad,lum=gaiaData(ra,dec)
    TeffSource,radSource,lumSource="Gaia Collaboration 2018","Gaia Collaboration 2018","Gaia Collaboration 2018"
    if Teff<=0:
        TeffSource=""
    starDict['Teff']=Teff
    starDict['TeffSoruce']=TeffSource
    if rad<=0:
        radSource=""
    starDict['Radius']=rad
    starDict['RadiusSoruce']=radSource
    if lum<=0:
        lumSource=""
    starDict['Luminosity']=lum
    starDict['LuminositySoruce']=lumSource
    
    return starDict

def gaiaData(ra,dec):
    coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.degree))
    width = u.Quantity(5, u.arcsecond)
    data=astroquery.gaia.Gaia.query_object(coordinate=coord, width=width, height=width)
    if len(data)==0:
        return -1,-1,-1
    Teff=data['teff_val'][0]
    if np.ma.is_masked(Teff):
        Teff=-1
    rad=data['radius_val'][0]
    if np.ma.is_masked(rad):
        rad=-1
    lum=data['lum_val'][0]
    if np.ma.is_masked(lum):
        lum=-1
    return Teff,rad,lum

def getCoords(name):
    simbad=astroquery.simbad.Simbad()
    star=simbad.query_object(name)[0]
    return star['RA'],star['DEC']

def getRef(ref): #finds the reference (in the form "Herczeg+ 2014") for a given bibcode (in form "2014ApJ...786...97H" or "J/ApJ/786/97")
    catalogs=astroquery.vizier.Vizier.find_catalogs(ref)
    description=catalogs[list(catalogs.items())[0][0]].description
    reference=description[description.rfind('(')+1:description.rfind(')')]
    return reference.replace(",","") #removes any commas
    
    
def getVizierSpectra(ra,dec,windowSize=2):
    #star=astroquery.simbad.Simbad.query_object(mainName)
    #starDict=getInfo(name=name,ra=ra,dec=dec)
    raString=str(ra).replace(' ','+').strip()
    decString=str(dec).replace(' ','+').strip()
    
    url="http://vizier.u-strasbg.fr/viz-bin/sed?-c="+raString+decString+"&-c.r="+str(windowSize)+"&-c.u=arcsec"
    data=astropy.table.Table.read(url)
    
    lambdas=c/np.array(1e9*data['sed_freq']) # converting to wavelength in m
    fluxs=1e-26*c*np.array(data['sed_flux'])/lambdas**2 # converting to F_lambda in W m
    sigmas=np.nan_to_num(1e-26*c*np.array(data['sed_eflux'])/lambdas**2) # converting to error in F_lambda in W m

    source=data['sed_filter']
    tables=data['_tabname']

    order=np.argsort(lambdas)

    lambdas=lambdas[order]
    fluxs=fluxs[order]
    sigmas=sigmas[order]
    source=source[order]
    tables=tables[order]
    
    sources=np.empty_like(lambdas,dtype=object)
    telescopes=np.empty_like(lambdas,dtype=object)
    
    for i,table in enumerate(tables):
        tableName = table[:table.rfind('/')]
        catalogs=astroquery.vizier.Vizier.find_catalogs(tableName)
        if len(catalogs)==0:
            continue #for some reason can't find table of data...
        description=catalogs[tableName].description
        reference=description[description.rfind('(')+1:description.rfind(')')]
        #print('ref: ',reference)
        #sources[i]='"'+reference.replace(',','')+'"'
        sources[i]=reference.replace(',','')
        
        # Manually maps each data point to the telescope it comes from (probably should tabulate this elsewhere and work from that...)
        if ('Herschel' in source[i]) | (reference == 'Marsh+ 2016'):
            telescopes[i]='Herschel'
        elif 'Spitzer' in source[i]:
            telescopes[i]='Spitzer'
        elif 'WISE' in source[i]:
            telescopes[i]='WISE'
        elif ('2MASS' in source[i]) | ('Johnson' in source[i]):
            telescopes[i]='2MASS'
        elif ('GAIA' in source[i]) | ('Gaia' in source[i]):
            telescopes[i]='Gaia'
        elif 'PAN' in source[i]:
            telescopes[i]='Pan-STARRS'
        elif ('IRAS' in source[i]) | (reference == 'Saunders+ 2000'):
            telescopes[i]='IRAS'
        elif 'POSS' in source[i]:
            telescopes[i]='POSS'
        elif 'AKARI' in source[i]:
            telescopes[i]='AKARI'
        elif 'SDSS' in source[i]:
            telescopes[i]='SDSS'
        elif ('SMA' in source[i]) | (reference == 'Andrews+ 2013') | (reference == 'Harris+ 2012'):
            telescopes[i]='SMA'
        elif 'Cousins' in source[i]:
            telescopes[i]='Cousins'
        elif ('ALMA' in source[i]) | (reference == 'Pascucci+ 2016'):
            telescopes[i]='ALMA'
        elif 'ISO' in source[i]:
            telescopes[i]='ISO'
        elif 'DENIS' in source[i]:
            telescopes[i]='DENIS'
        elif 'HIP' in source[i]:
            telescopes[i]='Hipparcos'
        elif 'GALEX' in source[i]:
            telescopes[i]='GALEX'
        elif 'XMM' in source[i]:
            telescopes[i]='XMM'
        elif ('SCUBA' in source[i]) | (reference == 'Mohanty+ 2013'):
            telescopes[i]='SCUBA'
        elif 'UKIDSS' in source[i]:
            telescopes[i]='UKIDSS'
        elif 'MKO' in source[i]:
            telescopes[i]='MKO'
        elif 'MSX' in source[i]:
            telescopes[i]='MSX'
        elif (reference == 'Dzib+ 2015'):
            telescopes[i]='VLA'
        elif (reference == 'Belloche+ 2011'):
            telescopes[i]='APEX'
        else:
            telescopes[i]='Unspecified'
            if source[i][0]!=':':
                print('___did you know about the ',source[i],' telescope?')
                
        if sources[i]=='Meng+ 2017': #points from this paper seem completely wrong (https://ui.adsabs.harvard.edu/#abs/2017ApJ...836...34M/abstract)
            telescopes[i]='Unspecified'
        #telescopes[i]='"'+telescopes[i]+'"'
    return lambdas,fluxs,sigmas,sources,telescopes

def queryIrsa(ra,dec,windowSize=2): #at the moment only queries the Herchel point source catalogs, but open to suggestions
    wavelengths=[70,100,160,250,350,500]
    coord=astropy.coordinates.SkyCoord(ra,dec,frame='icrs',unit=(u.hourangle, u.deg))
    ls=[]
    fs=[]
    es=[]
    for wl in wavelengths:
        if wl<200:
            catalog='ppsc_'+str(wl)
        else:
            catalog='spsc'+str(wl)
        table=astroquery.irsa.Irsa.query_region(coord,catalog=catalog,radius=windowSize*u.arcsecond,verbose=False)
        if len(table)==0:
            continue
        else:
            flux=table['flux'][0]
            error=flux/table['snr'][0]
            if 3*error > flux: # definition of upper limit
                flux=3*error
                error=0
            ls.append(wl*1.0e-6)
            fs.append(flux*c*1e-17/wl**2) #converting to F_lambda in SI
            es.append(error*c*1e-17/wl**2)
    return np.array(ls),np.array(fs),np.array(es)

def getIrsaSpectra(tmName):
    mainName,tmName,ra,dec,spType,objectType=getInfo(name='2MASS '+tmName)
    irsaLs,irsaFs,irsaEs=queryIrsa(ra,dec)
    if irsaLs[irsaLs<2e-4].size>0:
        addToSpectra(tmName,irsaLs[irsaLs<2e-4],irsaFs[irsaLs<2e-4],irsaEs[irsaLs<2e-4],2017,'Marton+ 2017','Herschel')
    if irsaLs[irsaLs>2e-4].size>0:
        addToSpectra(tmName,irsaLs[irsaLs>2e-4],irsaFs[irsaLs>2e-4],irsaEs[irsaLs>2e-4],2017,'Schulz+ 2017','Herschel')
        
def getSpectra(name=-1,ra=-1,dec=-1,saveDir=-1,windowSize=2,overwrite=0):
    starDict=getInfo(name=name,ra=ra,dec=dec)
    ra=starDict['RA']
    dec=starDict['DEC']
    
    if (saveDir!=-1) & (overwrite!=1):
        if saveDir[-1]!='/':
            saveDir=saveDir+'/'
        fName=saveDir+starDict['2MASSID']+'.ecsv'
        if Path(fName).is_file():
            if overwrite==0:
                print('Data file for this star already exists')
                print('You can find it at:')
                print(fName)
                print('Remember, all stars must have a 2MASS ID to be read, otherwise will default to nearest star in 2MASS catalog')
                print('The 2MASS ID of this star is: ',starDict['2MASSID'])
                print('Also known as ',starDict['SimbadName'],' at coordinates ',starDict['RA'],', ',starDict['DEC'])
                print('Returning the data from that file')
                print('To overwrite the file function use argument overwrite = 1 (default = 0)')
                print('Or to supress this warning set overwrite = -1')
            #table=astropy.io.ascii.read(fName)
            return getSpectraFromFile(starDict['2MASSID'],saveDir)
    
    columnNames=('lambda', 'flux', 'error','source','telescope')
    #HOW CAN I SAVE AS SCIENTIFIC NOTATION FOR A FIXED NUMBER OF SIGNIFICANT FIGURES???
    dataTypes=('f4','f4','f4','object','object') # note: strings must be treated as objects not strings to allow variable length 
    vLs,vFs,vEs,vSs,vTs = getVizierSpectra(ra,dec,windowSize=windowSize) #spectra from vizier
    table=astropy.table.Table([vLs,vFs,vEs,vSs,vTs],names=columnNames,dtype=dataTypes)
    
    iLs,iFs,iEs = queryIrsa(ra,dec,windowSize=windowSize)
    for i in range(iLs.size):
        telescope='Herschel'
        if iLs[i]>2e-4:
            source='Schulz+ 2017'
        else:
            source='Marton+ 2017'
        table.add_row([iLs[i],iFs[i],iEs[i],source,telescope])
    
    for key in starDict.keys():
        table.meta[key]=str(starDict[key])
    table.meta['FileCreated']=str(datetime.datetime.today()).split()[0]
    comments=['All units are SI, I leave it to the user to convert/add astropy units',
    'Meta data about the star is stored under indivdual fields',
    'e.g. if data stored in variable called "dataTable" the R.A. of the star can be found via "dataTable.meta["RA"]".',
    'Everything intended to be read into astropy tables - either directly or via the getSpectraFromFile() function.',
    'See GITHUB-REPO for more details.',
    'Please cite SOME-PAPER if you make use of this tool or data.']
    table.meta['comments']=comments
    
    if saveDir!=-1:
        if saveDir[-1]!='/':
            saveDir=saveDir+'/'
        fName=saveDir+starDict['2MASSID']+'.ecsv'
        astropy.io.ascii.write(table,output=fName, format='ecsv')
    return table
    
def getSpectraFromFile(twoMassID,saveDir):
    if saveDir[-1]!='/':
        saveDir=saveDir+'/'
    fName=saveDir+twoMassID+'.ecsv'
    if Path(fName).is_file():
        table=astropy.io.ascii.read(fName)
        #unfortunately all the meta data is converted to strings on saving, translating appropriate data back here
        table.meta['Distance']=float(table.meta['Distance'])
        table.meta['Teff']=float(table.meta['Teff'])
        table.meta['Radius']=float(table.meta['Radius'])
        table.meta['Luminosity']=float(table.meta['Luminosity'])
        return table
    else:
        print('No file found at: ')
        print(fName)
        return -1