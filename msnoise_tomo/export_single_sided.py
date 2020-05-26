from msnoise.api import *


def main():
    """

    :return:
    """
    db = connect()
    components_to_compute = get_components_to_compute(db)
    maxlag = float(get_config(db, "maxlag"))
    if not os.path.isdir("TOMO_SAC"):
        os.makedirs("TOMO_SAC")
    for station1, station2 in get_station_pairs(db, used=True):
        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)
        pair = "%s:%s" % (sta1, sta2)
        dist = get_interstation_distance(station1, station2,
                                         station1.coordinates)
        if dist == 0:
            logging.error("Interstation distance is 0.0 km, skipping this"
                          "pair: %s" % pair)
            continue
        if is_dtt_next_job(db, jobtype='TOMO_SAC', ref=pair):
            logging.info(
                "We will export the one-sided CCFs for REF: %s" % pair)

            for f in get_filters(db, all=False):
                filterid = int(f.ref)
                for components in components_to_compute:
                    ref_name = pair.replace('.', '_').replace(':', '_')
                    rf = os.path.join("STACKS", "%02i" %
                                      filterid, "REF", components, ref_name + ".MSEED")
                    ref = read(rf)[0]

                    st1 = ref.copy()
                    st1.data = st1.data[ref.stats.npts//2:]
                    st1.stats.starttime = 0
                    st1.stats.sac = AttribDict()
                    st1.stats.sac.b = 0
                    st1.stats.sac.e = maxlag
                    st1.stats.sac.stla = station1.Y
                    st1.stats.sac.stlo = station1.X
                    st1.stats.sac.evla = station2.Y
                    st1.stats.sac.evlo = station2.X
                    st1.stats.sac.scale = 1
                    st1.stats.sac.lcalda = 1
                    st1.stats.sac.dist = dist
                    st1.stats.sac.npts = st1.stats.npts
                    st1.stats.sac.kevnm = sta1
                    st1.stats.sac.kstnm = sta2
                    # fn1 = "%s_%s_%s_%s_REAL.sac"%(station1.net, station1.sta, station2.net, station2.sta)
                    # fn1 = os.path.join("TOMO_SAC", "%02i" % filterid, components, fn1)
                    # st1.write(fn1, format="SAC")

                    st2 = ref.copy()
                    st2.data = st2.data[:ref.stats.npts//2+1][::-1]
                    st2.stats.starttime = 0
                    st2.stats.sac = AttribDict()
                    st2.stats.sac.b = 0
                    st2.stats.sac.e = maxlag
                    st2.stats.sac.stla = station2.Y
                    st2.stats.sac.stlo = station2.X
                    st2.stats.sac.evla = station1.Y
                    st2.stats.sac.evlo = station1.X
                    st2.stats.sac.scale = 1
                    st2.stats.sac.lcalda = 1
                    st2.stats.sac.dist = get_interstation_distance(station1, station2, station1.coordinates)
                    st2.stats.sac.npts = st2.stats.npts
                    st2.stats.sac.kevnm = sta2
                    st2.stats.sac.kstnm = sta1
                    # fn2 = "%s_%s_%s_%s_REAL.sac"%(station2.net, station2.sta, station1.net, station1.sta)
                    # fn2 = os.path.join("TOMO_SAC", "%02i" % filterid, components, fn2)
                    # st2.write(fn2, format="SAC")

                    st3 = st1
                    st3.data += st2.data
                    st3.data /= 2.0
                    fn3 = "%s_%s_%s_%s_MEAN.sac"%(station1.net, station1.sta, station2.net, station2.sta)
                    outpath = os.path.join("TOMO_SAC", "%02i" % filterid, components)
                    if not os.path.isdir(outpath):
                        os.makedirs(outpath)
                    fn3 = os.path.join(outpath, fn3)
                    st3.write(fn3, format="SAC")

            update_job(db, "REF", pair, jobtype='TOMO_SAC', flag='D')
            update_job(db, "REF", pair, jobtype='TOMO_FTAN', flag='T')


if __name__ == "__main__":
    main()