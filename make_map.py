import exp_map as em

bands,colormaps,coverage=em.bands()

date='20190331'

exp_file='exposures_'+date+'.csv'

exposures=em.get_exp_info(outfile=exp_file)

#title='DECam '+b+'-band sky coverage - Mar 31, 2019',

for b in bands:
    coverage[b] = em.make_exp_map(data_frame=exposures,
                                  title=None, 
                                  color=colormaps[b],
                                  band=b,
                                  out_file='decam_footprint_'+b+'_'+date+'.pdf',
                                  show_ra_labels_in_hours=True,
                                  width_of_mw_band=0.,
                                  minteff=em.teffcut(b))


