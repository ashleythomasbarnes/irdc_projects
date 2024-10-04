from imports import *

#
def get_Reff(A): 
    """Get effective radius from geom. area
    A = props['area_exact']
    """
    R_eff = (A/np.pi)**0.5
    return(R_eff/au.arcsec)

#
def get_bgsubflux(leaf, bpix):
    """Get background subtracted flux of leaf"""
    
    Iv = leaf.values()
    min_Iv = leaf.vmin
    Iv_bgsub = Iv - min_Iv
    sum_Iv = np.nansum(Iv_bgsub)
    Sv = sum_Iv / bpix 
    
    return(Sv)

#
def get_bgsubfluxes(leaves, bpix):
    """Get background subtracted flux of leaves
    Uses get_bgsubflux module"""
    
    nleaves = len(leaves)
    bgsubflux = np.empty(nleaves) *np.nan
    
    for leaf, i in zip(leaves, range(nleaves)):
        bgsubflux[i] = get_bgsubflux(leaf, bpix)
    
    bgsubflux = bgsubflux *au.Jy
        
    return(bgsubflux)

#
def get_statsflux(props, data, index_map, unit=au.Jy):
    """Return max flux within each leaf"""
    maxflux = np.empty(len(props['_idx'])) *np.nan
    meanflux = np.empty(len(props['_idx'])) *np.nan
    
    for idx in props['_idx']:
        mask = np.where(index_map==idx)
        data_masked = data[mask]
        idx_ = np.where(props['_idx'] == idx)[0]
        
        maxflux[idx_] = np.nanmax(data_masked) 
        meanflux[idx_] = np.nanmean(data_masked) 
        
    return(maxflux*unit, meanflux*unit)

#
def pruneleaves(index_map, props):
    """Return index array containing only leaves"""
    for x in range(index_map.shape[0]):
        for y in range(index_map.shape[1]):
            if index_map[x,y] not in list(props['_idx']):
                index_map[x,y] = -1
            
    return(index_map)

def get_cropmap(hdu_, region, verbose=True):
    
    hdu_tmp = hdu_.copy() #Creating dummy file
    hdu_tmp.data = [hdu_tmp.data, hdu_tmp.data]

    hdu_tmp.header['CDELT3'] = 1
    hdu_tmp.header['CRPIX3'] = 1
    hdu_tmp.header['CRVAL3'] = 1
    hdu_tmp.header['CUNIT3'] = 'km/s'
    hdu_tmp.header['CTYPE3'] = 'VELO'

    cube = SpectralCube.read(hdu_tmp)
    shape_old = cube.hdu.shape
    minax = np.nanmin(shape_old[1:])
    
    cube = cube.subcube_from_ds9region(region)
    cube = cube[:,0:minax,0:minax]
    
    hdu_tmp = cube.hdu 
    hdu_tmp.data = hdu_tmp.data[0]
    del hdu_tmp.header['*3']
    hdu_tmp.header['WCSAXES'] = 2

    shape_new = hdu_tmp.shape
    
    if verbose: 
        print('[INFO] Using region - %s' %region) 
        print('[INFO] Shape: [%s,%s] ---> [%s,%s]' %(shape_old[1], shape_old[2], shape_new[0], shape_new[1]))

    hdu_= hdu_tmp
    
    return(hdu_)